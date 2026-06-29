#!/usr/bin/env python3
"""Prepare TPXO M2 SSH + barotropic velocity OBCS records for the LR test.

This implements the Ponte & Cornuelle style boundary experiment: sea-surface
height and barotropic velocity are put on the same M2 clock and read through
standard OBCS prescribed fields.  The existing `useOBCStides` amplitude/phase
path is deliberately not used here because it has no matching SSH harmonic
interface in this configuration.

The output files contain 16 evenly spaced records over one M2 period.  T/S
boundary profiles are duplicated to the same record count so all prescribed
OBCS fields share `externForcingPeriod=2794.5` and
`externForcingCycle=44712.`.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr
from scipy.interpolate import RegularGridInterpolator


TPXO_VELOCITY_FILE = Path("/Users/Karina/Research/Canyons/data/tides/u_tpxo9.v5a.nc")
TPXO_ELEVATION_FILE = Path("/Users/Karina/Research/Canyons/data/tides/h_tpxo9.v5a.nc")

CONSTITUENT = "M2"
PERIOD = 44712.0
N_RECORDS = 16
TAG = "M2_fullTPXO_LR220x240_telescopic20_obeta16"

HERE = Path(__file__).resolve().parent
EXPERIMENT = HERE.parent
GRID_DIR = EXPERIMENT / "runBC_LR_fullTPXO_balanced"
OUTPUT_DIR = EXPERIMENT / "input"

NX = 220
NY = 240
NR = 60


def read_mds_grid(name: str, levels: int = 1) -> np.ndarray:
    values = np.fromfile(GRID_DIR / f"{name}.data", dtype=">f4")
    expected = levels * NY * NX
    if values.size != expected:
        raise ValueError(f"{name}: expected {expected} values, got {values.size}")
    if levels == 1:
        return values.reshape(NY, NX).astype(np.float64)
    return values.reshape(levels, NY, NX).astype(np.float64)


def constituent_index(dataset: xr.Dataset, name: str) -> int:
    names = [value.decode().strip().upper() for value in dataset.con.values]
    try:
        return names.index(name.upper())
    except ValueError as exc:
        raise ValueError(f"TPXO constituent {name!r} not found in {names}") from exc


def complex_interpolator(
    dataset: xr.Dataset,
    amp_name: str,
    phase_name: str,
    lon_name: str,
    lat_name: str,
    constituent: str,
    target_lon: np.ndarray,
    target_lat: np.ndarray,
    amp_scale: float = 1.0,
) -> np.ndarray:
    """Interpolate A*exp(-i*phase) onto target points."""
    index = constituent_index(dataset, constituent)
    lon = np.asarray(dataset[lon_name][:, 0])
    lat = np.asarray(dataset[lat_name][0, :])

    target_x = np.mod(target_lon, 360.0)
    target_y = target_lat
    buffer = 0.5
    ix = np.flatnonzero((lon >= target_x.min() - buffer) & (lon <= target_x.max() + buffer))
    iy = np.flatnonzero((lat >= target_y.min() - buffer) & (lat <= target_y.max() + buffer))
    if ix.size < 2 or iy.size < 2:
        raise ValueError(f"{amp_name}: TPXO subset does not bracket target boundary")

    amplitude = np.asarray(
        dataset[amp_name].isel(nc=index, nx=slice(ix[0], ix[-1] + 1), ny=slice(iy[0], iy[-1] + 1)),
        dtype=np.float64,
    ) * amp_scale
    phase = np.asarray(
        dataset[phase_name].isel(nc=index, nx=slice(ix[0], ix[-1] + 1), ny=slice(iy[0], iy[-1] + 1)),
        dtype=np.float64,
    )
    harmonic = amplitude * np.exp(-1j * np.deg2rad(phase))

    interpolator = RegularGridInterpolator(
        (lon[ix], lat[iy]), harmonic, method="linear", bounds_error=True
    )
    return interpolator(np.column_stack((target_x, target_y)))


def make_records_1d(harmonic: np.ndarray) -> np.ndarray:
    omega_t = 2 * np.pi * np.arange(N_RECORDS, dtype=np.float64) / N_RECORDS
    return np.real(harmonic[None, :] * np.exp(1j * omega_t[:, None]))


def make_velocity_records(harmonic: np.ndarray, wet_mask: np.ndarray) -> np.ndarray:
    """Return records in MDS section order: record, k, horizontal-index."""
    surface = make_records_1d(harmonic)
    records = np.repeat(surface[:, None, :], NR, axis=1)
    records *= wet_mask[None, :, :]
    return records


def write_be64(path: Path, array: np.ndarray) -> None:
    array.astype(">f8").tofile(path)
    print(f"wrote {path}  shape={array.shape}  bytes={path.stat().st_size}")


def duplicate_boundary_profile(source_name: str, horizontal_size: int) -> None:
    source = OUTPUT_DIR / source_name
    destination = OUTPUT_DIR / source_name.replace(".bin", "_obeta16.bin")
    values = np.fromfile(source, dtype=">f8")
    record_size = NR * horizontal_size
    if values.size % record_size != 0:
        raise ValueError(f"{source}: {values.size} values is not an integer number of records")
    first = values[:record_size].reshape(NR, horizontal_size)
    records = np.repeat(first[None, :, :], N_RECORDS, axis=0)
    write_be64(destination, records)


def main() -> None:
    xc = read_mds_grid("XC")
    yc = read_mds_grid("YC")
    xg = read_mds_grid("XG")
    yg = read_mds_grid("YG")
    hfac_c = read_mds_grid("hFacC", NR)
    hfac_w = read_mds_grid("hFacW", NR)
    hfac_s = read_mds_grid("hFacS", NR)

    wet_w_u = (hfac_w[:, :, 1] > 0).astype(np.float64)
    wet_w_v = (hfac_c[:, :, 0] > 0).astype(np.float64)
    wet_s_v = (hfac_s[:, 1, :] > 0).astype(np.float64)
    wet_s_u = (hfac_c[:, 0, :] > 0).astype(np.float64)
    wet_n_v = (hfac_s[:, -1, :] > 0).astype(np.float64)
    wet_n_u = (hfac_c[:, -1, :] > 0).astype(np.float64)

    with xr.open_dataset(TPXO_ELEVATION_FILE) as elev:
        eta_w = complex_interpolator(elev, "ha", "hp", "lon_z", "lat_z", CONSTITUENT, xc[:, 0], yc[:, 0])
        eta_s = complex_interpolator(elev, "ha", "hp", "lon_z", "lat_z", CONSTITUENT, xc[0, :], yc[0, :])
        eta_n = complex_interpolator(elev, "ha", "hp", "lon_z", "lat_z", CONSTITUENT, xc[-1, :], yc[-1, :])

    with xr.open_dataset(TPXO_VELOCITY_FILE) as vel:
        # TPXO velocity amplitudes are cm/s.
        u_w = complex_interpolator(vel, "ua", "up", "lon_u", "lat_u", CONSTITUENT, xg[:, 1], yc[:, 1], 0.01)
        v_w = complex_interpolator(vel, "va", "vp", "lon_v", "lat_v", CONSTITUENT, xc[:, 0], yg[:, 0], 0.01)
        u_s = complex_interpolator(vel, "ua", "up", "lon_u", "lat_u", CONSTITUENT, xg[0, :], yc[0, :], 0.01)
        v_s = complex_interpolator(vel, "va", "vp", "lon_v", "lat_v", CONSTITUENT, xc[1, :], yg[1, :], 0.01)
        u_n = complex_interpolator(vel, "ua", "up", "lon_u", "lat_u", CONSTITUENT, xg[-1, :], yc[-1, :], 0.01)
        v_n = complex_interpolator(vel, "va", "vp", "lon_v", "lat_v", CONSTITUENT, xc[-1, :], yg[-1, :], 0.01)

    write_be64(OUTPUT_DIR / f"OBWeta_{TAG}.bin", make_records_1d(eta_w))
    write_be64(OUTPUT_DIR / f"OBSeta_{TAG}.bin", make_records_1d(eta_s))
    write_be64(OUTPUT_DIR / f"OBNeta_{TAG}.bin", make_records_1d(eta_n))

    write_be64(OUTPUT_DIR / f"OBWu_{TAG}.bin", make_velocity_records(u_w, wet_w_u))
    write_be64(OUTPUT_DIR / f"OBWv_{TAG}.bin", make_velocity_records(v_w, wet_w_v))
    write_be64(OUTPUT_DIR / f"OBSu_{TAG}.bin", make_velocity_records(u_s, wet_s_u))
    write_be64(OUTPUT_DIR / f"OBSv_{TAG}.bin", make_velocity_records(v_s, wet_s_v))
    write_be64(OUTPUT_DIR / f"OBNu_{TAG}.bin", make_velocity_records(u_n, wet_n_u))
    write_be64(OUTPUT_DIR / f"OBNv_{TAG}.bin", make_velocity_records(v_n, wet_n_v))

    for name, size in [
        ("linTN_LR220x240_telescopic20.bin", NX),
        ("linTS_LR220x240_telescopic20.bin", NX),
        ("linSN_LR220x240_telescopic20.bin", NX),
        ("linSS_LR220x240_telescopic20.bin", NX),
        ("linTW_LR220x240_telescopic20.bin", NY),
        ("linSW_LR220x240_telescopic20.bin", NY),
    ]:
        duplicate_boundary_profile(name, size)

    print("Done. Configure data with externForcingPeriod=2794.5 and externForcingCycle=44712.")


if __name__ == "__main__":
    main()
