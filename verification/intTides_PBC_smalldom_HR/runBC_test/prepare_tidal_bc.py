#!/usr/bin/env python3
"""Prepare transport-consistent MITgcm tidal open-boundary velocities.

The TPXO velocity harmonics are interpolated as complex numbers directly onto
the native MITgcm U/V boundary faces.  For this idealized experiment the
outward complex transport is then constrained to zero using the minimum
area-weighted, depth-uniform velocity correction shared by all open faces.

This deliberately removes the domain-mean elevation tide implied by TPXO while
retaining its spatially varying incident barotropic current as closely as
possible.  It is not a reconstruction of the complete TPXO elevation solution.
"""

from pathlib import Path

import numpy as np
import xarray as xr
from scipy.interpolate import RegularGridInterpolator


# User configuration ---------------------------------------------------------
TPXO_VELOCITY_FILE = Path(
    "/Users/Karina/Research/Canyons/data/tides/u_tpxo9.v5a.nc"
)
CONSTITUENTS = {"M2": 44712.0}  # exact periods used by data.obcs [s]
BALANCE_MODE = "zero_net"       # appropriate for the intended idealized run


HERE = Path(__file__).resolve().parent
EXPERIMENT = HERE.parent
GRID_DIR = EXPERIMENT / "runTSini_mmtm"
OUTPUT_DIR = EXPERIMENT / "input"
NX = NY = 600
NR = 60
DEL_R = np.array([
    1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5,
    5, 6, 6, 7, 8, 8, 9, 10, 11, 12, 13, 14, 16, 17, 19, 21, 22, 25,
    27, 29, 32, 35, 38, 42, 46, 50, 54, 59, 65, 71, 77, 84, 92, 101,
    110, 120, 131, 143, 153, 168,
], dtype=np.float64)


def read_mds_grid(name: str, levels: int = 1) -> np.ndarray:
    """Read the existing float32 MITgcm grid output in (k,j,i) order."""
    values = np.fromfile(GRID_DIR / f"{name}.data", dtype=">f4")
    expected = levels * NY * NX
    if values.size != expected:
        raise ValueError(f"{name}: expected {expected} values, got {values.size}")
    if levels == 1:
        return values.reshape(NY, NX).astype(np.float64)
    return values.reshape(levels, NY, NX).astype(np.float64)


def boundary_geometry() -> dict[str, dict[str, np.ndarray]]:
    """Return exact face coordinates and wet cross-sectional areas."""
    xc, yc = read_mds_grid("XC"), read_mds_grid("YC")
    xg, yg = read_mds_grid("XG"), read_mds_grid("YG")
    dxg, dyg = read_mds_grid("DXG"), read_mds_grid("DYG")
    hfac_w = read_mds_grid("hFacW", NR)
    hfac_s = read_mds_grid("hFacS", NR)

    # These indices follow data.obcs and MITgcm's OBCS application stencil:
    # OB_Iwest=1 -> active U face i=2; OB_Jsouth=1 -> V face j=2;
    # OB_Jnorth=Ny -> V face j=Ny (Fortran one-based indices).
    west_area = np.sum(hfac_w[:, :, 1] * DEL_R[:, None], axis=0) * dyg[:, 0]
    south_area = np.sum(hfac_s[:, 1, :] * DEL_R[:, None], axis=0) * dxg[1, :]
    north_area = np.sum(hfac_s[:, -1, :] * DEL_R[:, None], axis=0) * dxg[-1, :]

    return {
        "W": {"lon": xg[:, 1], "lat": yc[:, 1], "area": west_area},
        "S": {"lon": xc[1, :], "lat": yg[1, :], "area": south_area},
        "N": {"lon": xc[-1, :], "lat": yg[-1, :], "area": north_area},
    }


def constituent_index(dataset: xr.Dataset, name: str) -> int:
    names = [value.decode().strip().upper() for value in dataset.con.values]
    try:
        return names.index(name.upper())
    except ValueError as exc:
        raise ValueError(f"TPXO constituent {name!r} not found in {names}") from exc


def complex_interpolator(
    dataset: xr.Dataset, component: str, constituent: str, target_lon: np.ndarray,
    target_lat: np.ndarray,
) -> np.ndarray:
    """Linearly interpolate A*exp(-i*phase) rather than wrapped phase."""
    index = constituent_index(dataset, constituent)
    suffix = "u" if component == "u" else "v"
    x = np.asarray(dataset[f"lon_{suffix}"][:, 0])
    y = np.asarray(dataset[f"lat_{suffix}"][0, :])

    target_x = np.mod(target_lon, 360.0)
    buffer = 0.5
    ix = np.flatnonzero((x >= target_x.min() - buffer) &
                        (x <= target_x.max() + buffer))
    iy = np.flatnonzero((y >= target_lat.min() - buffer) &
                        (y <= target_lat.max() + buffer))
    if ix.size < 2 or iy.size < 2:
        raise ValueError("TPXO subset does not bracket the MITgcm boundary")

    amp_name, phase_name = ("ua", "up") if component == "u" else ("va", "vp")
    amplitude = np.asarray(dataset[amp_name].isel(
        nc=index, nx=slice(ix[0], ix[-1] + 1), ny=slice(iy[0], iy[-1] + 1)
    ), dtype=np.float64) * 0.01  # cm/s -> m/s
    phase = np.asarray(dataset[phase_name].isel(
        nc=index, nx=slice(ix[0], ix[-1] + 1), ny=slice(iy[0], iy[-1] + 1)
    ), dtype=np.float64)
    harmonic = amplitude * np.exp(-1j * np.deg2rad(phase))

    interpolator = RegularGridInterpolator(
        (x[ix], y[iy]), harmonic, method="linear", bounds_error=True
    )
    points = np.column_stack((target_x, target_lat))
    return interpolator(points)


def outward_transport(coefficients: dict[str, np.ndarray], geometry: dict) -> complex:
    """Complex volume transport, positive outward from the model domain."""
    return (
        np.sum(-coefficients["W"] * geometry["W"]["area"])
        + np.sum(-coefficients["S"] * geometry["S"]["area"])
        + np.sum(+coefficients["N"] * geometry["N"]["area"])
    )


def impose_zero_transport(coefficients: dict[str, np.ndarray], geometry: dict) -> complex:
    """Apply the least-squares uniform outward correction over all wet faces."""
    initial = outward_transport(coefficients, geometry)
    total_area = sum(np.sum(item["area"]) for item in geometry.values())
    outward_correction = -initial / total_area

    wet_w = geometry["W"]["area"] > 0
    wet_s = geometry["S"]["area"] > 0
    wet_n = geometry["N"]["area"] > 0
    coefficients["W"][wet_w] -= outward_correction
    coefficients["S"][wet_s] -= outward_correction
    coefficients["N"][wet_n] += outward_correction
    coefficients["W"][~wet_w] = 0
    coefficients["S"][~wet_s] = 0
    coefficients["N"][~wet_n] = 0
    return outward_correction


def amplitude_and_phase_seconds(harmonic: np.ndarray, period: float) -> tuple:
    amplitude = np.abs(harmonic)
    phase_radians = np.mod(-np.angle(harmonic), 2 * np.pi)
    phase_seconds = phase_radians * period / (2 * np.pi)
    phase_seconds[amplitude == 0] = 0
    return amplitude, phase_seconds


def write_records(path: Path, records: list[np.ndarray]) -> None:
    """Write one full boundary vector per constituent in MDS record order."""
    np.stack(records, axis=0).astype(">f8").tofile(path)


def main() -> None:
    geometry = boundary_geometry()
    output: dict[str, list[np.ndarray]] = {
        "W_amp": [], "W_phase": [], "S_amp": [], "S_phase": [],
        "N_amp": [], "N_phase": [],
    }

    with xr.open_dataset(TPXO_VELOCITY_FILE) as dataset:
        for name, period in CONSTITUENTS.items():
            coefficients = {
                "W": complex_interpolator(dataset, "u", name,
                                            geometry["W"]["lon"], geometry["W"]["lat"]),
                "S": complex_interpolator(dataset, "v", name,
                                            geometry["S"]["lon"], geometry["S"]["lat"]),
                "N": complex_interpolator(dataset, "v", name,
                                            geometry["N"]["lon"], geometry["N"]["lat"]),
            }
            for side in coefficients:
                coefficients[side][geometry[side]["area"] <= 0] = 0

            before = outward_transport(coefficients, geometry)
            correction = 0j
            if BALANCE_MODE == "zero_net":
                correction = impose_zero_transport(coefficients, geometry)
            after = outward_transport(coefficients, geometry)

            print(f"{name}: net transport before = {abs(before):.6g} m3/s "
                  f"at phase {np.angle(before, deg=True):.3f} deg")
            print(f"{name}: uniform outward correction = {abs(correction):.6g} m/s "
                  f"at phase {np.angle(correction, deg=True):.3f} deg")
            print(f"{name}: net transport after  = {abs(after):.6g} m3/s")

            for side in ("W", "S", "N"):
                amp, phase = amplitude_and_phase_seconds(coefficients[side], period)
                output[f"{side}_amp"].append(amp)
                output[f"{side}_phase"].append(phase)

    tag = "_".join(CONSTITUENTS) + "_zeroNet"
    names = {
        "W_amp": f"OBWuam_{tag}.bin", "W_phase": f"OBWuph_{tag}.bin",
        "S_amp": f"OBSvam_{tag}.bin", "S_phase": f"OBSvph_{tag}.bin",
        "N_amp": f"OBNvam_{tag}.bin", "N_phase": f"OBNvph_{tag}.bin",
    }
    for key, filename in names.items():
        path = OUTPUT_DIR / filename
        write_records(path, output[key])
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
