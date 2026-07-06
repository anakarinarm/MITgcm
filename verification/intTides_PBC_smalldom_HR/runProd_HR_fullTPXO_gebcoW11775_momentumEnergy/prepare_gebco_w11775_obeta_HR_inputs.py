#!/usr/bin/env python3
"""Prepare HR TPXO SSH+velocity OBCS inputs with GEBCO bathymetry.

This run tests the raw full-TPXO SSH+velocity boundary architecture on an HR
grid whose inner physical-domain southwest corner is near
117.75 W, 31.25 N.  A 60-cell stretched sponge is added outside the physical
domain on W/S/N.  The physical-domain cells are kept close to square.

Bathymetry is interpolated from the broad local GEBCO file instead of
edge-extending the previous small-domain bathymetry.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr
from scipy.interpolate import RegularGridInterpolator


GEBCO_FILE = Path(
    "/Users/karina/Research/Canyons/data/bathy/"
    "GEBCO_07_Jun_2024_180ba79203c1/gebco_2023_n34.5_s29.5_w-120.0_e-115.0.nc"
)
TPXO_VELOCITY_FILE = Path("/Users/karina/Research/Canyons/data/tides/u_tpxo9.v5a.nc")
TPXO_ELEVATION_FILE = Path("/Users/karina/Research/Canyons/data/tides/h_tpxo9.v5a.nc")

CONSTITUENT = "M2"
PERIOD = 44712.0
N_RECORDS = 16
NR = 60
SPONGE = 60
PHYSICAL_NX = 690
PHYSICAL_NY = 600
NX = PHYSICAL_NX + SPONGE
NY = PHYSICAL_NY + 2 * SPONGE
MIN_WET_DEPTH = 10.0
CLOSE_EASTERN_CELLS = 1

OLD_HR_TAG = "660x720_telescopic60"
TAG = "HR750x720_gebcoW11775_physW11775_tel60"
OBCS_TAG = f"M2_fullTPXO_{TAG}_obeta16"

PHYSICAL_XG_WEST = 242.25
PHYSICAL_XG_EAST = 243.601840000000
PHYSICAL_YG_SOUTH = 31.25

HERE = Path(__file__).resolve().parent
EXPERIMENT = HERE.parent
INPUT = EXPERIMENT / "input"


DRF = np.array(
    [
        1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5,
        5, 6, 6, 7, 8, 8, 9, 10, 11, 12, 13, 14, 16, 17, 19, 21, 22,
        25, 27, 29, 32, 35, 38, 42, 46, 50, 54, 59, 65, 71, 77, 84,
        92, 101, 110, 120, 131, 143, 153, 168,
    ],
    dtype=np.float64,
)


def read_be64(path: Path, shape: tuple[int, ...]) -> np.ndarray:
    values = np.fromfile(path, dtype=">f8")
    expected = int(np.prod(shape))
    if values.size != expected:
        raise ValueError(f"{path}: expected {expected} values, got {values.size}")
    return values.reshape(shape).astype(np.float64)


def write_be64(path: Path, array: np.ndarray) -> None:
    np.asarray(array, dtype=">f8").tofile(path)
    print(f"wrote {path.name}: shape={array.shape}, bytes={path.stat().st_size}")


def build_widths() -> tuple[np.ndarray, np.ndarray, float, float]:
    dy = read_be64(INPUT / f"smalldom_{OLD_HR_TAG}_dy.bin", (NY,))
    if NY != PHYSICAL_NY + 2 * SPONGE:
        raise ValueError("NY must equal PHYSICAL_NY + 2*SPONGE for this grid design")

    lat_ref = PHYSICAL_YG_SOUTH + 0.5 * dy[SPONGE:-SPONGE].sum()
    core_dx = (PHYSICAL_XG_EAST - PHYSICAL_XG_WEST) / PHYSICAL_NX
    # Match the physical width of the S/N sponge cells in the W sponge.  The
    # south-sponge sequence is ordered boundary -> interior, which is exactly
    # what we need for a west sponge ordered west boundary -> physical domain.
    west_sponge_dx = dy[:SPONGE] / np.cos(np.deg2rad(lat_ref))
    dx = np.r_[west_sponge_dx, np.full(PHYSICAL_NX, core_dx, dtype=np.float64)]
    xg_origin = PHYSICAL_XG_WEST - west_sponge_dx.sum()
    yg_origin = PHYSICAL_YG_SOUTH - dy[:SPONGE].sum()

    write_be64(INPUT / f"smalldom_{TAG}_dx.bin", dx)
    write_be64(INPUT / f"smalldom_{TAG}_dy.bin", dy)
    print(f"model grid west/east: {xg_origin:.12f} E to {xg_origin + dx.sum():.12f} E")
    print(f"model grid south/north: {yg_origin:.12f} N to {yg_origin + dy.sum():.12f} N")
    print(f"physical west/east: {PHYSICAL_XG_WEST:.12f} E to {PHYSICAL_XG_EAST:.12f} E")
    print(f"physical south/north: {PHYSICAL_YG_SOUTH:.12f} N to {PHYSICAL_YG_SOUTH + dy[SPONGE:-SPONGE].sum():.12f} N")
    print(f"physical core dx/dy ratio in meters: {core_dx * np.cos(np.deg2rad(lat_ref)) / np.median(dy[SPONGE:-SPONGE]):.4f}")
    return dx, dy, xg_origin, yg_origin


def centers_edges(origin: float, widths: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    edges = origin + np.r_[0.0, np.cumsum(widths)]
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, edges[:-1]


def horizontal_grid(
    dx: np.ndarray,
    dy: np.ndarray,
    xg_origin: float,
    yg_origin: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    xc_1d, xg_1d = centers_edges(xg_origin, dx)
    yc_1d, yg_1d = centers_edges(yg_origin, dy)
    xc = np.repeat(xc_1d[None, :], NY, axis=0)
    xg = np.repeat(xg_1d[None, :], NY, axis=0)
    yc = np.repeat(yc_1d[:, None], NX, axis=1)
    yg = np.repeat(yg_1d[:, None], NX, axis=1)
    return xc, yc, xg, yg


def flatten_sponge_bathymetry(bathy: np.ndarray) -> None:
    """Remove bathymetric gradients normal to W/S/N sponge buffers."""
    if SPONGE <= 0:
        return
    if bathy.shape != (NY, NX):
        raise ValueError(f"Expected bathy shape {(NY, NX)}, got {bathy.shape}")
    bathy[:, :SPONGE] = bathy[:, [SPONGE]]
    bathy[:SPONGE, :] = bathy[[SPONGE], :]
    bathy[NY - SPONGE :, :] = bathy[[NY - SPONGE - 1], :]
    print("flattened bathymetry normal to the W/S/N sponge buffers")


def interpolate_gebco_bathymetry(xc: np.ndarray, yc: np.ndarray) -> np.ndarray:
    with xr.open_dataset(GEBCO_FILE) as ds:
        lon = np.asarray(ds["lon"], dtype=np.float64)
        lat = np.asarray(ds["lat"], dtype=np.float64)
        elev = np.asarray(ds["elevation"], dtype=np.float64)
    interp = RegularGridInterpolator(
        (lat, lon),
        elev,
        method="linear",
        bounds_error=False,
        fill_value=0.0,
    )
    lon_w = np.where(xc > 180.0, xc - 360.0, xc)
    elev_model = interp(np.column_stack((yc.ravel(), lon_w.ravel()))).reshape(xc.shape)
    bathy = np.where(elev_model < 0.0, elev_model, 0.0)
    wet = bathy < 0.0
    shallow = wet & (-bathy < MIN_WET_DEPTH)
    bathy[shallow] = -MIN_WET_DEPTH
    flatten_sponge_bathymetry(bathy)
    if CLOSE_EASTERN_CELLS > 0:
        bathy[:, -CLOSE_EASTERN_CELLS:] = 0.0
        print(
            "closed eastern boundary by setting the last "
            f"{CLOSE_EASTERN_CELLS} C-grid column(s) to land"
        )
    wet = bathy < 0.0
    print(f"GEBCO wet cells: {int(wet.sum())}/{bathy.size}")
    print(f"deepened wet cells shallower than {MIN_WET_DEPTH:g} m: {int(shallow.sum())}")
    if wet.any():
        print(f"wet depth min/max: {float((-bathy[wet]).min()):.3f}/{float((-bathy[wet]).max()):.3f} m")
    write_be64(INPUT / f"smalldom_{TAG}_bat_min10m_obeta.bin", bathy)
    return bathy


def prepare_hydrography() -> None:
    for stem in ("IniTemp", "IniSalt"):
        old = np.memmap(INPUT / f"{stem}{NR}x660x720_telescopic60.bin", dtype=">f8", mode="r", shape=(NR, NY, 660))
        profile = np.asarray(old[:, 0, 0], dtype=np.float64)
        out = np.memmap(INPUT / f"{stem}{NR}x{NX}x{NY}_{TAG}.bin", dtype=">f8", mode="w+", shape=(NR, NY, NX))
        for k, value in enumerate(profile):
            out[k, :, :] = value
        out.flush()
        del out
        del old
        print(f"wrote {stem}{NR}x{NX}x{NY}_{TAG}.bin")


def duplicate_boundary_profiles() -> None:
    for source_name, source_n, target_n in [
        ("linTN_telescopic60.bin", 660, NX),
        ("linTS_telescopic60.bin", 660, NX),
        ("linSN_telescopic60.bin", 660, NX),
        ("linSS_telescopic60.bin", 660, NX),
        ("linTW_telescopic60.bin", NY, NY),
        ("linSW_telescopic60.bin", NY, NY),
    ]:
        source = INPUT / source_name
        values = read_be64(source, (2, NR, source_n))
        profile = values[:1, :, :1]
        extended = np.repeat(profile, target_n, axis=2)
        records = np.repeat(extended, N_RECORDS, axis=0)
        out_name = source_name.replace(".bin", f"_{TAG}_obeta16.bin")
        write_be64(INPUT / out_name, records)


def constituent_index(dataset: xr.Dataset, name: str) -> int:
    names = [value.decode().strip().upper() for value in dataset.con.values]
    try:
        return names.index(name.upper())
    except ValueError as exc:
        raise ValueError(f"TPXO constituent {name!r} not found in {names}") from exc


def tpxo_complex_interpolator(
    dataset: xr.Dataset,
    amp_name: str,
    phase_name: str,
    lon_name: str,
    lat_name: str,
    target_lon: np.ndarray,
    target_lat: np.ndarray,
    amp_scale: float = 1.0,
) -> np.ndarray:
    index = constituent_index(dataset, CONSTITUENT)
    lon = np.asarray(dataset[lon_name][:, 0], dtype=np.float64)
    lat = np.asarray(dataset[lat_name][0, :], dtype=np.float64)
    target_x = np.mod(target_lon, 360.0)
    target_y = target_lat
    buffer = 0.5
    ix = np.flatnonzero((lon >= target_x.min() - buffer) & (lon <= target_x.max() + buffer))
    iy = np.flatnonzero((lat >= target_y.min() - buffer) & (lat <= target_y.max() + buffer))
    if ix.size < 2 or iy.size < 2:
        raise ValueError(f"{amp_name}: TPXO subset does not bracket target boundary")
    amp = np.asarray(dataset[amp_name].isel(nc=index, nx=slice(ix[0], ix[-1] + 1), ny=slice(iy[0], iy[-1] + 1)), dtype=np.float64)
    phase = np.asarray(dataset[phase_name].isel(nc=index, nx=slice(ix[0], ix[-1] + 1), ny=slice(iy[0], iy[-1] + 1)), dtype=np.float64)
    harmonic = amp_scale * amp * np.exp(-1j * np.deg2rad(phase))
    interp = RegularGridInterpolator(
        (lon[ix], lat[iy]),
        harmonic,
        method="linear",
        bounds_error=False,
        fill_value=0.0 + 0.0j,
    )
    return interp(np.column_stack((target_x, target_y)))


def make_records_1d(harmonic: np.ndarray) -> np.ndarray:
    omega_t = 2 * np.pi * np.arange(N_RECORDS, dtype=np.float64) / N_RECORDS
    return np.real(harmonic[None, :] * np.exp(1j * omega_t[:, None]))


def vertical_wet_mask(section_depth: np.ndarray) -> np.ndarray:
    z_top = np.r_[0.0, np.cumsum(DRF[:-1])]
    depth = np.maximum(section_depth, 0.0)
    return (depth[None, :] > z_top[:, None]).astype(np.float64)


def make_velocity_records(harmonic: np.ndarray, section_depth: np.ndarray) -> np.ndarray:
    surface = make_records_1d(harmonic)
    records = np.repeat(surface[:, None, :], NR, axis=1)
    records *= vertical_wet_mask(section_depth)[None, :, :]
    return records


def prepare_tpxo_forcing(xc: np.ndarray, yc: np.ndarray, xg: np.ndarray, yg: np.ndarray, bathy: np.ndarray) -> None:
    depth = np.where(bathy < 0, -bathy, 0.0)
    with xr.open_dataset(TPXO_ELEVATION_FILE) as elev:
        eta_w = tpxo_complex_interpolator(elev, "ha", "hp", "lon_z", "lat_z", xc[:, 0], yc[:, 0])
        eta_s = tpxo_complex_interpolator(elev, "ha", "hp", "lon_z", "lat_z", xc[0, :], yc[0, :])
        eta_n = tpxo_complex_interpolator(elev, "ha", "hp", "lon_z", "lat_z", xc[-1, :], yc[-1, :])

    with xr.open_dataset(TPXO_VELOCITY_FILE) as vel:
        # TPXO current amplitudes are cm/s.
        u_w = tpxo_complex_interpolator(vel, "ua", "up", "lon_u", "lat_u", xg[:, 1], yc[:, 1], 0.01)
        v_w = tpxo_complex_interpolator(vel, "va", "vp", "lon_v", "lat_v", xc[:, 0], yg[:, 0], 0.01)
        u_s = tpxo_complex_interpolator(vel, "ua", "up", "lon_u", "lat_u", xg[0, :], yc[0, :], 0.01)
        v_s = tpxo_complex_interpolator(vel, "va", "vp", "lon_v", "lat_v", xc[1, :], yg[1, :], 0.01)
        u_n = tpxo_complex_interpolator(vel, "ua", "up", "lon_u", "lat_u", xg[-1, :], yc[-1, :], 0.01)
        v_n = tpxo_complex_interpolator(vel, "va", "vp", "lon_v", "lat_v", xc[-1, :], yg[-1, :], 0.01)

    write_be64(INPUT / f"OBWeta_{OBCS_TAG}.bin", make_records_1d(eta_w))
    write_be64(INPUT / f"OBSeta_{OBCS_TAG}.bin", make_records_1d(eta_s))
    write_be64(INPUT / f"OBNeta_{OBCS_TAG}.bin", make_records_1d(eta_n))

    write_be64(INPUT / f"OBWu_{OBCS_TAG}.bin", make_velocity_records(u_w, depth[:, 1]))
    write_be64(INPUT / f"OBWv_{OBCS_TAG}.bin", make_velocity_records(v_w, depth[:, 0]))
    write_be64(INPUT / f"OBSu_{OBCS_TAG}.bin", make_velocity_records(u_s, depth[0, :]))
    write_be64(INPUT / f"OBSv_{OBCS_TAG}.bin", make_velocity_records(v_s, depth[1, :]))
    write_be64(INPUT / f"OBNu_{OBCS_TAG}.bin", make_velocity_records(u_n, depth[-1, :]))
    write_be64(INPUT / f"OBNv_{OBCS_TAG}.bin", make_velocity_records(v_n, depth[-1, :]))


def main() -> None:
    dx, dy, xg_origin, yg_origin = build_widths()
    xc, yc, xg, yg = horizontal_grid(dx, dy, xg_origin, yg_origin)
    bathy = interpolate_gebco_bathymetry(xc, yc)
    prepare_hydrography()
    duplicate_boundary_profiles()
    prepare_tpxo_forcing(xc, yc, xg, yg, bathy)
    print("Done.")


if __name__ == "__main__":
    main()
