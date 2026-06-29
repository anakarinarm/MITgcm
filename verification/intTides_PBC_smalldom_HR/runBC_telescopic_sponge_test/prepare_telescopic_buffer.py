#!/usr/bin/env python3
"""Build a telescopic numerical buffer without adding physical structure.

The original 600 x 600 domain is embedded unchanged at [j=60:660, i=60:660].
Sixty cells are added west, south, and north.  Their widths grow by 2.5%
per cell toward the outer open boundary.  Bathymetry and hydrography are
extended normally from the old boundary. Existing complex M2 harmonics are
retained over the original coordinate spans and tapered to zero over the new
corner segments. A minimum, area-weighted barotropic correction based on the
initialized stretched-grid metrics restores zero complex boundary transport.
"""

from pathlib import Path
import argparse

import numpy as np


HERE = Path(__file__).resolve().parent
EXPERIMENT = HERE.parent
INPUT = EXPERIMENT / "input"

OLD_NX = OLD_NY = 600
NR = 60
BUFFER = 60
NX = OLD_NX + BUFFER
NY = OLD_NY + 2 * BUFFER
STRETCH = 1.025
M2_PERIOD = 44712.0

DX_OLD_FILE = INPUT / "smalldom_600x600_BTS_dx.bin"
DY_OLD_FILE = INPUT / "smalldom_600x600_BTS_dy.bin"
BATHY_OLD_FILE = INPUT / "smalldom_600x600_BTS_bat_finite.bin"
TEMP_OLD_FILE = INPUT / "IniTemp60x600x600.bin"
SALT_OLD_FILE = INPUT / "IniSalt60x600x600.bin"

TAG = "telescopic60"


def read_be64(path: Path, shape: tuple[int, ...]) -> np.ndarray:
    values = np.fromfile(path, dtype=">f8")
    expected = int(np.prod(shape))
    if values.size != expected:
        raise ValueError(f"{path}: expected {expected} values, got {values.size}")
    return values.reshape(shape).astype(np.float64)


def write_be64(path: Path, values: np.ndarray) -> None:
    np.asarray(values, dtype=">f8").tofile(path)
    print(f"wrote {path.name}: shape={values.shape}")


def stretched_outer_cells(interface_width: float) -> np.ndarray:
    """Widths ordered from the outer boundary toward the old domain."""
    return interface_width * STRETCH ** np.arange(BUFFER - 1, -1, -1)


def build_grid_and_state() -> tuple[np.ndarray, np.ndarray]:
    dx_old = read_be64(DX_OLD_FILE, (OLD_NX,))
    dy_old = read_be64(DY_OLD_FILE, (OLD_NY,))
    dx_west = stretched_outer_cells(dx_old[0])
    dy_south = stretched_outer_cells(dy_old[0])
    dy_north = stretched_outer_cells(dy_old[-1])[::-1]
    dx = np.concatenate((dx_west, dx_old))
    dy = np.concatenate((dy_south, dy_old, dy_north))

    write_be64(INPUT / f"smalldom_{NX}x{NY}_{TAG}_dx.bin", dx)
    write_be64(INPUT / f"smalldom_{NX}x{NY}_{TAG}_dy.bin", dy)

    bathy_old = read_be64(BATHY_OLD_FILE, (OLD_NY, OLD_NX))
    bathy = np.pad(bathy_old, ((BUFFER, BUFFER), (BUFFER, 0)), mode="edge")
    write_be64(INPUT / f"smalldom_{NX}x{NY}_{TAG}_bat.bin", bathy)

    # These initial fields are horizontally uniform at every level.  Extend
    # them explicitly so the buffer starts with exactly the same stratification.
    for old_path, stem in ((TEMP_OLD_FILE, "IniTemp"),
                           (SALT_OLD_FILE, "IniSalt")):
        old = np.memmap(old_path, dtype=">f8", mode="r",
                        shape=(NR, OLD_NY, OLD_NX))
        profile = np.asarray(old[:, 0, 0], dtype=np.float64)
        if not all(np.all(old[k] == profile[k]) for k in range(NR)):
            raise ValueError(f"{old_path} is not horizontally uniform")
        output = INPUT / f"{stem}{NR}x{NX}x{NY}_{TAG}.bin"
        mm = np.memmap(output, dtype=">f8", mode="w+", shape=(NR, NY, NX))
        for k, value in enumerate(profile):
            mm[k, :, :] = value
        mm.flush()
        del mm
        print(f"wrote {output.name}: shape={(NR, NY, NX)}")

    return dx_old, dy_old


def extend_boundary_tracers() -> None:
    # Files contain two (time, depth, horizontal) records with identical
    # background profiles.  Edge extension therefore introduces no anomaly.
    for side, old_n, pads in (
        ("W", OLD_NY, (BUFFER, BUFFER)),
        ("S", OLD_NX, (BUFFER, 0)),
        ("N", OLD_NX, (BUFFER, 0)),
    ):
        for tracer in ("T", "S"):
            source = INPUT / f"lin{tracer}{side}.bin"
            values = read_be64(source, (2, NR, old_n))
            extended = np.pad(values, ((0, 0), (0, 0), pads), mode="edge")
            write_be64(INPUT / f"lin{tracer}{side}_{TAG}.bin", extended)


def read_new_boundary_areas() -> dict[str, np.ndarray]:
    """Read exact wet face areas produced by MITgcm on the stretched grid."""
    required = [HERE / f"{name}.data" for name in
                ("DRF", "DXG", "DYG", "hFacW", "hFacS")]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Exact tide balancing requires a zero-step initialization first; "
            f"missing {missing}"
        )
    drf = np.fromfile(HERE / "DRF.data", dtype=">f4").astype(np.float64)
    dxg = np.fromfile(HERE / "DXG.data", dtype=">f4").reshape(NY, NX)
    dyg = np.fromfile(HERE / "DYG.data", dtype=">f4").reshape(NY, NX)
    hfw = np.fromfile(HERE / "hFacW.data", dtype=">f4").reshape(NR, NY, NX)
    hfs = np.fromfile(HERE / "hFacS.data", dtype=">f4").reshape(NR, NY, NX)
    return {
        "W": np.sum(hfw[:, :, 1] * drf[:, None], axis=0) * dyg[:, 0],
        "S": np.sum(hfs[:, 1, :] * drf[:, None], axis=0) * dxg[1, :],
        "N": np.sum(hfs[:, -1, :] * drf[:, None], axis=0) * dxg[-1, :],
    }


def read_harmonic(side: str) -> np.ndarray:
    component = "u" if side == "W" else "v"
    amp = read_be64(INPUT / f"OB{side}{component}am_M2_zeroNet.bin", (OLD_NX,))
    phase = read_be64(INPUT / f"OB{side}{component}ph_M2_zeroNet.bin", (OLD_NX,))
    return amp * np.exp(-2j * np.pi * phase / M2_PERIOD)


def outward_transport(coeff: dict[str, np.ndarray],
                      areas: dict[str, np.ndarray]) -> complex:
    return (-np.sum(coeff["W"] * areas["W"])
            - np.sum(coeff["S"] * areas["S"])
            + np.sum(coeff["N"] * areas["N"]))


def corner_taper() -> np.ndarray:
    """Raised-cosine weight: zero at the outer corner, one at old-domain edge."""
    phase = np.linspace(0.0, 0.5 * np.pi, BUFFER)
    return np.sin(phase) ** 2


def extend_tides() -> None:
    old = {side: read_harmonic(side) for side in ("W", "S", "N")}
    taper = corner_taper()
    coeff = {
        "W": np.concatenate((old["W"][0] * taper,
                              old["W"],
                              old["W"][-1] * taper[::-1])),
        "S": np.concatenate((old["S"][0] * taper, old["S"])),
        "N": np.concatenate((old["N"][0] * taper, old["N"])),
    }
    areas = read_new_boundary_areas()
    before = outward_transport(coeff, areas)
    total_area = sum(np.sum(area) for area in areas.values())
    outward_correction = -before / total_area
    wet = {side: areas[side] > 0 for side in areas}
    coeff["W"][wet["W"]] -= outward_correction
    coeff["S"][wet["S"]] -= outward_correction
    coeff["N"][wet["N"]] += outward_correction
    for side in coeff:
        coeff[side][~wet[side]] = 0.0
    after = outward_transport(coeff, areas)

    print(f"extended-boundary transport before correction: {abs(before):.9g} m3/s")
    print(f"uniform outward correction: {abs(outward_correction):.9g} m/s")
    print(f"extended-boundary transport after correction: {abs(after):.9g} m3/s")

    for side in ("W", "S", "N"):
        component = "u" if side == "W" else "v"
        amplitude = np.abs(coeff[side])
        phase = np.mod(-np.angle(coeff[side]), 2 * np.pi) * M2_PERIOD / (2 * np.pi)
        phase[amplitude == 0] = 0.0
        write_be64(INPUT / f"OB{side}{component}am_M2_zeroNet_{TAG}.bin", amplitude)
        write_be64(INPUT / f"OB{side}{component}ph_M2_zeroNet_{TAG}.bin", phase)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--grid-only", action="store_true",
        help="write grid/state inputs before the required zero-step initialization",
    )
    args = parser.parse_args()
    dx_old, dy_old = build_grid_and_state()
    extend_boundary_tracers()
    if not args.grid_only:
        extend_tides()
    else:
        print("grid-only stage complete; initialize MITgcm, then rerun without "
              "--grid-only to balance and write tidal forcing")
    west_width = np.sum(stretched_outer_cells(dx_old[0]))
    south_width = np.sum(stretched_outer_cells(dy_old[0]))
    print(f"new xgOrigin = {242.5 - west_width:.12f} degrees east")
    print(f"new ygOrigin = {31.25 - south_width:.12f} degrees north")
    print(f"original domain indices (Python): j={BUFFER}:{BUFFER+OLD_NY}, "
          f"i={BUFFER}:{BUFFER+OLD_NX}")


if __name__ == "__main__":
    main()
