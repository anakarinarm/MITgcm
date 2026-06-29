#!/usr/bin/env python3
"""Modal characteristic split of M2 boundary-normal energy.

This is a companion to ``check_boundary_reflections.py``.  It uses the same
MITgcm diagnostics (UVEL, VVEL, PHIHYD, ETAN), fits the post-ramp M2 harmonic,
projects the baroclinic pressure potential and normal velocity onto simple
hydrostatic vertical modes, and estimates outward/inward characteristic energy
for each mode.

The decomposition is intentionally diagnostic, not a replacement for a full
modal radiation boundary condition.  The vertical modes use a WKB/cosine
approximation based on the initial T/S stratification and the model's linear
EOS.  The result is most useful for answering: "is this boundary dominated by
prescribed incoming tide, or by outward internal-wave radiation/reflection?"
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class Meta:
    nx: int
    ny: int
    nz: int
    dtype: np.dtype
    nrecords: int
    fields: list[str]
    time: float | None


def parse_meta(path: Path) -> Meta:
    text = path.read_text()
    dim_block = re.search(r"dimList\s*=\s*\[(.*?)\];", text, re.S)
    if not dim_block:
        raise ValueError(f"Cannot parse dimList from {path}")
    numbers = [int(x) for x in re.findall(r"[-+]?\d+", dim_block.group(1))]
    dims = numbers[0::3]
    if len(dims) == 2:
        nx, ny = dims
        nz = 1
    elif len(dims) == 3:
        nx, ny, nz = dims
    else:
        raise ValueError(f"Unexpected dimensions {dims} in {path}")
    precision = re.search(r"dataprec\s*=\s*\[\s*'([^']+)'", text).group(1)
    dtype = np.dtype(">f4" if precision.strip() == "float32" else ">f8")
    nrecords = int(re.search(r"nrecords\s*=\s*\[\s*(\d+)", text).group(1))
    field_block = re.search(r"fldList\s*=\s*\{(.*?)\};", text, re.S)
    fields = re.findall(r"'([^']+)'", field_block.group(1)) if field_block else []
    fields = [field.strip() for field in fields]
    time_match = re.search(r"timeInterval\s*=\s*\[\s*([^\]]+)", text)
    time = float(time_match.group(1)) if time_match else None
    return Meta(nx, ny, nz, dtype, nrecords, fields, time)


def read_mds(path: Path) -> tuple[np.ndarray, Meta]:
    meta_path = path.with_suffix(".meta") if path.suffix == ".data" else path
    data_path = meta_path.with_suffix(".data")
    meta = parse_meta(meta_path)
    values = np.fromfile(data_path, dtype=meta.dtype)
    expected = meta.nrecords * meta.nz * meta.ny * meta.nx
    if values.size != expected:
        raise ValueError(f"{data_path}: expected {expected} values, got {values.size}")
    if not np.isfinite(values).all():
        bad = values.size - np.count_nonzero(np.isfinite(values))
        raise ValueError(f"{data_path}: contains {bad} non-finite values")
    return values.reshape(meta.nrecords, meta.nz, meta.ny, meta.nx), meta


def field_records(array: np.ndarray, meta: Meta, field: str, levels: int) -> np.ndarray:
    if field not in meta.fields:
        raise ValueError(f"{field} absent from fields {meta.fields}")
    if meta.nrecords != len(meta.fields) or meta.nz != levels:
        raise ValueError(
            f"Expected {len(meta.fields)} records with {levels} levels, "
            f"got nrecords={meta.nrecords}, nz={meta.nz}"
        )
    return array[meta.fields.index(field)]


def grid_field(run_dir: Path, name: str) -> np.ndarray:
    meta = parse_meta(run_dir / f"{name}.meta")
    raw = np.fromfile(run_dir / f"{name}.data", dtype=meta.dtype)
    expected = meta.nz * meta.ny * meta.nx
    if raw.size != expected:
        raise ValueError(f"{name}: expected {expected} values, got {raw.size}")
    result = raw.reshape(meta.nz, meta.ny, meta.nx)
    return result if meta.nz > 1 else result[0]


def fit_harmonic(values: np.ndarray, times: np.ndarray, period: float) -> np.ndarray:
    centered = (times - times.mean()) / period
    omega_t = 2 * np.pi * times / period
    design = np.column_stack((np.ones(times.size), centered,
                              np.cos(omega_t), np.sin(omega_t)))
    coefficients = np.linalg.pinv(design) @ values.reshape(times.size, -1)
    harmonic = coefficients[2] - 1j * coefficients[3]
    return harmonic.reshape(values.shape[1:])


def remove_depth_mean(values: np.ndarray, vertical_weights: np.ndarray) -> np.ndarray:
    denominator = np.sum(vertical_weights, axis=0)
    numerator = np.sum(values * vertical_weights[None, ...], axis=1)
    mean = np.divide(numerator, denominator[None, :],
                     out=np.zeros_like(numerator), where=denominator[None, :] > 0)
    return values - mean[:, None, :]


def namelist_value(data_text: str, name: str, default: float | None = None) -> float:
    match = re.search(rf"\b{name}\s*=\s*([^,\n]+)", data_text)
    if not match:
        if default is None:
            raise ValueError(f"{name} not found in data")
        return default
    return float(match.group(1).replace("D", "E").replace("d", "e"))


def namelist_string(data_text: str, name: str) -> str:
    match = re.search(rf"\b{name}\s*=\s*'([^']+)'", data_text)
    if not match:
        raise ValueError(f"{name} not found in data")
    return match.group(1)


def initial_profile(run_dir: Path, filename: str, nz: int, ny: int, nx: int) -> np.ndarray:
    path = (run_dir / filename).resolve()
    nvalues = nz * ny * nx
    bytes_per_value = path.stat().st_size // nvalues
    if bytes_per_value == 4:
        dtype = ">f4"
    elif bytes_per_value == 8:
        dtype = ">f8"
    else:
        raise ValueError(f"Cannot infer precision for {path}")
    raw = np.memmap(path, dtype=dtype, mode="r", shape=(nz, ny, nx))
    # These files are horizontally uniform in this experiment; average a tiny
    # corner sample to avoid reading the whole 200+ MB array into memory.
    return np.asarray(raw[:, : min(4, ny), : min(4, nx)].mean(axis=(1, 2)))


def stratification_and_modes(run_dir: Path, nmodes: int, rho0: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data_text = (run_dir / "data").read_text()
    alpha = namelist_value(data_text, "tAlpha")
    beta = namelist_value(data_text, "sBeta")
    gravity = namelist_value(data_text, "gravity", 9.81)
    theta_file = namelist_string(data_text, "hydrogThetaFile")
    salt_file = namelist_string(data_text, "hydrogSaltFile")

    drf = np.asarray(grid_field(run_dir, "DRF")).reshape(-1)
    rc = np.asarray(grid_field(run_dir, "RC")).reshape(-1)
    nz = drf.size
    rac = grid_field(run_dir, "RAC")
    ny, nx = rac.shape
    theta = initial_profile(run_dir, theta_file, nz, ny, nx)
    salt = initial_profile(run_dir, salt_file, nz, ny, nx)

    z_down = -rc
    density_anomaly = rho0 * (-alpha * (theta - theta[0]) + beta * (salt - salt[0]))
    drho_dz = np.gradient(density_anomaly, z_down)
    n2 = gravity / rho0 * drho_dz
    n2 = np.maximum(n2, 1e-8)
    nfreq = np.sqrt(n2)
    buoyancy_integral = np.sum(nfreq * drf)
    depth = np.sum(drf)

    modes = []
    speeds = []
    stretched = np.cumsum(drf) - 0.5 * drf
    s = stretched / depth
    weights = drf
    for mode in range(1, nmodes + 1):
        shape = np.cos(mode * np.pi * s)
        shape -= np.sum(shape * weights) / np.sum(weights)
        norm = math.sqrt(np.sum(shape * shape * weights) / np.sum(weights))
        shape /= norm
        modes.append(shape)
        speeds.append(buoyancy_integral / (mode * np.pi))
    return np.asarray(modes), np.asarray(speeds), n2


def load_harmonics(run_dir: Path, offsets: list[int], period: float,
                   ramp_periods: float) -> tuple[dict, dict, np.ndarray]:
    files_3d = sorted(run_dir.glob("bc_wave_3d.[0-9]*.meta"))
    files_eta = sorted(run_dir.glob("bc_wave_eta.[0-9]*.meta"))
    if not files_3d or not files_eta:
        raise FileNotFoundError("No bc_wave_3d/bc_wave_eta diagnostics found")
    if len(files_3d) != len(files_eta):
        raise ValueError("3-D and ETAN diagnostic counts differ")

    drf = np.asarray(grid_field(run_dir, "DRF")).reshape(-1)
    hfac = np.asarray(grid_field(run_dir, "hFacC"))
    dxc, dyc = grid_field(run_dir, "DXC"), grid_field(run_dir, "DYC")
    rac = grid_field(run_dir, "RAC")
    nr, ny, nx = hfac.shape
    sections = {side: {offset: {"p": [], "u": []} for offset in offsets}
                for side in ("W", "S", "N")}
    times = []

    for file_3d, file_eta in zip(files_3d, files_eta):
        array_3d, meta_3d = read_mds(file_3d)
        array_eta, meta_eta = read_mds(file_eta)
        if abs(meta_3d.time - meta_eta.time) > 0.1:
            raise ValueError("3-D and ETAN snapshots are not synchronous")
        u = field_records(array_3d, meta_3d, "UVEL", nr)
        v = field_records(array_3d, meta_3d, "VVEL", nr)
        phi = field_records(array_3d, meta_3d, "PHIHYD", nr)
        eta = field_records(array_eta, meta_eta, "ETAN", 1)[0]
        times.append(meta_3d.time)
        for offset in offsets:
            iw = 1 + offset
            js = 1 + offset
            jn = ny - 1 - offset
            sections["W"][offset]["p"].append(phi[:, :, iw] + 9.81 * eta[:, iw])
            sections["W"][offset]["u"].append(-0.5 * (u[:, :, iw] + u[:, :, iw + 1]))
            sections["S"][offset]["p"].append(phi[:, js, :] + 9.81 * eta[js, :])
            sections["S"][offset]["u"].append(-0.5 * (v[:, js, :] + v[:, js + 1, :]))
            sections["N"][offset]["p"].append(phi[:, jn, :] + 9.81 * eta[jn, :])
            sections["N"][offset]["u"].append(+0.5 * (v[:, jn, :] + v[:, jn + 1, :]))

    times = np.asarray(times)
    keep = times >= ramp_periods * period - 0.1
    if np.count_nonzero(keep) < 8:
        raise ValueError("Fewer than eight post-ramp samples are available")

    harmonics = {side: {} for side in ("W", "S", "N")}
    for side in ("W", "S", "N"):
        for offset in offsets:
            pressure = np.asarray(sections[side][offset]["p"])[keep]
            velocity = np.asarray(sections[side][offset]["u"])[keep]
            if side == "W":
                i = 1 + offset
                hfac_line = hfac[:, :, i]
                horizontal = dyc[:, i]
                distance = np.sum(np.nanmean(dxc[:, 1:i + 1], axis=0))
            elif side == "S":
                j = 1 + offset
                hfac_line = hfac[:, j, :]
                horizontal = dxc[j, :]
                distance = np.sum(np.nanmean(dyc[1:j + 1, :], axis=1))
            else:
                j = ny - 1 - offset
                hfac_line = hfac[:, j, :]
                horizontal = dxc[j, :]
                distance = np.sum(np.nanmean(dyc[j:ny - 1, :], axis=1))
            vertical = hfac_line * drf[:, None]
            pressure = remove_depth_mean(pressure, vertical)
            velocity = remove_depth_mean(velocity, vertical)
            harmonics[side][offset] = {
                "p": fit_harmonic(pressure, times[keep], period),
                "u": fit_harmonic(velocity, times[keep], period),
                "vertical": vertical,
                "horizontal": horizontal,
                "distance_km": distance / 1000,
            }
    grid = {"drf": drf, "hfac": hfac, "rac": rac}
    return harmonics, grid, times[keep]


def modal_projection(
    profile: np.ndarray, vertical: np.ndarray, modes: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    nmodes = modes.shape[0]
    out = np.zeros((nmodes, profile.shape[1]), dtype=np.complex128)
    norms = np.zeros((nmodes, profile.shape[1]))
    for n in range(nmodes):
        mode = modes[n, :, None]
        denom = np.sum(vertical * mode * mode, axis=0)
        numer = np.sum(profile * vertical * mode, axis=0)
        out[n] = np.divide(numer, denom, out=np.zeros_like(numer), where=denom > 0)
        norms[n] = denom
    return out, norms


def analyze(args: argparse.Namespace) -> list[dict]:
    run_dir = args.run_dir.resolve()
    modes, speeds, n2 = stratification_and_modes(run_dir, args.modes, args.rho0)
    harmonics, _, _ = load_harmonics(run_dir, args.offsets, args.period, args.ramp_periods)

    print("Estimated WKB long-wave mode speeds (m/s): " +
          ", ".join(f"c{n + 1}={c:.3f}" for n, c in enumerate(speeds)))
    print(f"N2 range used for modes: {np.min(n2):.3e} to {np.max(n2):.3e} s^-2")

    rows = []
    for side in ("W", "S", "N"):
        for offset in args.offsets:
            item = harmonics[side][offset]
            p_modes, p_norms = modal_projection(item["p"], item["vertical"], modes)
            u_modes, u_norms = modal_projection(item["u"], item["vertical"], modes)
            modal_norms = 0.5 * (p_norms + u_norms)
            wet_width = np.where(np.sum(item["vertical"], axis=0) > 0,
                                 item["horizontal"], 0.0)
            for n, c in enumerate(speeds, start=1):
                p_over_c = p_modes[n - 1] / c
                u_out = 0.5 * (u_modes[n - 1] + p_over_c)
                u_in = 0.5 * (p_over_c - u_modes[n - 1])
                weights = wet_width * modal_norms[n - 1]
                flux_out = 0.5 * args.rho0 * c * np.sum(weights * np.abs(u_out) ** 2)
                flux_in = 0.5 * args.rho0 * c * np.sum(weights * np.abs(u_in) ** 2)
                net = flux_out - flux_in
                ratio = math.sqrt(flux_in / flux_out) if flux_out > 0 else math.nan
                incoming_fraction = flux_in / (flux_in + flux_out) if flux_in + flux_out > 0 else math.nan
                rows.append({
                    "boundary": side,
                    "offset_cells": offset,
                    "distance_km": item["distance_km"],
                    "mode": n,
                    "phase_speed_m_s": c,
                    "outward_characteristic_W": flux_out,
                    "inward_characteristic_W": flux_in,
                    "net_outward_W": net,
                    "inward_to_outward_amplitude": ratio,
                    "incoming_energy_fraction": incoming_fraction,
                })
    return rows


def write_results(run_dir: Path, rows: list[dict]) -> None:
    csv_path = run_dir / "modal_characteristic_summary.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {csv_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--period", type=float, default=44712.0)
    parser.add_argument("--ramp-periods", type=float, default=3.0)
    parser.add_argument("--rho0", type=float, default=1026.0)
    parser.add_argument("--modes", type=int, default=5)
    parser.add_argument("--offsets", type=lambda value: [int(x) for x in value.split(",")],
                        default=[22, 23, 27, 33])
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    output = analyze(arguments)
    write_results(arguments.run_dir.resolve(), output)
