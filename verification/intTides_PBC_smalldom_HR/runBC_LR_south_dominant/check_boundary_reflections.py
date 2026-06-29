#!/usr/bin/env python3
"""Diagnose internal-wave reflection at MITgcm open boundaries.

The script fits the complex M2 harmonic after the forcing ramp, removes the
depth-mean velocity and pressure, and evaluates outward baroclinic
pressure-work flux on several sections parallel to each open boundary.

It also reports a normalized progressive-wave index chi. Under the local,
single-mode/equal-impedance approximation, the reflected-to-outgoing amplitude
ratio is sqrt((1-chi)/(1+chi)). This is a useful reflection proxy, not an exact
multimode characteristic decomposition. Robust conclusions require the proxy
and outward flux to be consistent across several nearby sections.
"""

from __future__ import annotations

import argparse
import csv
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
        raise ValueError(
            f"{data_path}: contains {bad} non-finite values; the model run "
            "is numerically invalid and reflection metrics cannot be computed"
        )
    return values.reshape(meta.nrecords, meta.nz, meta.ny, meta.nx), meta


def diagnostic_files(run_dir: Path, prefix: str) -> list[Path]:
    return sorted(run_dir.glob(f"{prefix}.[0-9]*.meta"))


def field_records(array: np.ndarray, meta: Meta, field: str, levels: int) -> np.ndarray:
    if field not in meta.fields:
        raise ValueError(f"{field} absent from fields {meta.fields}")
    if meta.nrecords != len(meta.fields) or meta.nz != levels:
        raise ValueError(
            f"Expected {len(meta.fields)} field records with {levels} levels, "
            f"but metadata reports nrecords={meta.nrecords}, nz={meta.nz}. "
            "Check that 2-D and 3-D diagnostics are in separate output streams."
        )
    index = meta.fields.index(field)
    return array[index]


def grid_field(run_dir: Path, name: str) -> np.ndarray:
    meta = parse_meta(run_dir / f"{name}.meta")
    raw = np.fromfile(run_dir / f"{name}.data", dtype=meta.dtype)
    # Standard MITgcm grid files encode z in dimList rather than nrecords.
    expected = meta.nz * meta.ny * meta.nx
    if raw.size != expected:
        raise ValueError(f"{name}: expected {expected} grid values, got {raw.size}")
    result = raw.reshape(meta.nz, meta.ny, meta.nx)
    return result if meta.nz > 1 else result[0]


def fit_harmonic(values: np.ndarray, times: np.ndarray, period: float) -> np.ndarray:
    """Fit constant, trend, cos and sin; return C for Re(C exp(i omega t))."""
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


def load_sections(run_dir: Path, offsets: list[int], period: float,
                  ramp_periods: float) -> tuple[dict, np.ndarray, np.ndarray, dict]:
    files_3d = diagnostic_files(run_dir, "bc_wave_3d")
    files_eta = diagnostic_files(run_dir, "bc_wave_eta")
    if not files_3d or not files_eta:
        old = diagnostic_files(run_dir, "bc_wave_state")
        if old:
            _, old_meta = read_mds(old[0])
            levels = old_meta.nrecords // max(len(old_meta.fields), 1)
            raise RuntimeError(
                "The completed bc_wave_state output contains only "
                f"{levels} vertical level. ETAN shared its stream with 3-D fields, "
                "so MITgcm selected their common surface level. Rerun with the "
                "corrected data.diagnostics (bc_wave_3d + bc_wave_eta)."
            )
        raise FileNotFoundError("No bc_wave_3d/bc_wave_eta diagnostics found")
    if len(files_3d) != len(files_eta):
        raise ValueError("The 3-D and ETAN diagnostic file counts differ")

    drf = np.asarray(grid_field(run_dir, "DRF")).reshape(-1)
    nr = drf.size
    hfac = np.asarray(grid_field(run_dir, "hFacC"))
    dxc, dyc = grid_field(run_dir, "DXC"), grid_field(run_dir, "DYC")
    rac = grid_field(run_dir, "RAC")
    ny, nx = rac.shape

    sections = {side: {offset: {"p": [], "u": []} for offset in offsets}
                for side in ("W", "S", "N")}
    eta_means, times = [], []
    wet_surface = hfac[0] > 0
    area_sum = np.sum(rac[wet_surface])

    for file_3d, file_eta in zip(files_3d, files_eta):
        array_3d, meta_3d = read_mds(file_3d)
        array_eta, meta_eta = read_mds(file_eta)
        if meta_3d.time is None or meta_eta.time is None:
            raise ValueError("Diagnostic metadata has no timeInterval")
        if abs(meta_3d.time - meta_eta.time) > 0.1:
            raise ValueError("3-D and ETAN snapshots are not synchronous")
        u = field_records(array_3d, meta_3d, "UVEL", nr)
        v = field_records(array_3d, meta_3d, "VVEL", nr)
        phi = field_records(array_3d, meta_3d, "PHIHYD", nr)
        eta = field_records(array_eta, meta_eta, "ETAN", 1)[0]
        times.append(meta_3d.time)
        eta_means.append(np.sum(eta[wet_surface] * rac[wet_surface]) / area_sum)

        for offset in offsets:
            iw = 1 + offset
            js = 1 + offset
            jn = ny - 1 - offset
            if iw + 1 >= nx or js + 1 >= ny or jn < 0:
                raise ValueError(f"Offset {offset} lies outside the domain")

            # Hydrostatic pressure potential and velocity interpolated to C cells.
            sections["W"][offset]["p"].append(phi[:, :, iw] + 9.81 * eta[:, iw])
            sections["W"][offset]["u"].append(-0.5 * (u[:, :, iw] + u[:, :, iw + 1]))
            sections["S"][offset]["p"].append(phi[:, js, :] + 9.81 * eta[js, :])
            sections["S"][offset]["u"].append(-0.5 * (v[:, js, :] + v[:, js + 1, :]))
            sections["N"][offset]["p"].append(phi[:, jn, :] + 9.81 * eta[jn, :])
            sections["N"][offset]["u"].append(+0.5 * (v[:, jn, :] + v[:, jn + 1, :]))

    times = np.asarray(times)
    eta_means = np.asarray(eta_means)
    keep = times >= ramp_periods * period - 0.1
    if np.count_nonzero(keep) < 8:
        raise ValueError("Fewer than eight post-ramp samples are available")

    geometry = {"drf": drf, "hfac": hfac, "dxc": dxc, "dyc": dyc,
                "nx": nx, "ny": ny, "keep": keep}
    return sections, times, eta_means, geometry


def analyze(args: argparse.Namespace) -> list[dict]:
    run_dir = args.run_dir.resolve()
    sections, times, eta_means, grid = load_sections(
        run_dir, args.offsets, args.period, args.ramp_periods
    )
    keep = grid["keep"]
    post_times = times[keep]
    results = []

    for side in ("W", "S", "N"):
        for offset in args.offsets:
            pressure = np.asarray(sections[side][offset]["p"])[keep]
            velocity = np.asarray(sections[side][offset]["u"])[keep]
            if side == "W":
                i = 1 + offset
                hfac_line = grid["hfac"][:, :, i]
                horizontal = grid["dyc"][:, i]
                # Sum actual metric widths: offset*local-dx is wrong on a
                # telescopic grid. Index 1 is the active OBCS-adjacent cell.
                distance = np.sum(np.nanmean(grid["dxc"][:, 1:i + 1], axis=0))
            elif side == "S":
                j = 1 + offset
                hfac_line = grid["hfac"][:, j, :]
                horizontal = grid["dxc"][j, :]
                distance = np.sum(np.nanmean(grid["dyc"][1:j + 1, :], axis=1))
            else:
                j = grid["ny"] - 1 - offset
                hfac_line = grid["hfac"][:, j, :]
                horizontal = grid["dxc"][j, :]
                distance = np.sum(np.nanmean(
                    grid["dyc"][j:grid["ny"] - 1, :], axis=1))

            vertical = hfac_line * grid["drf"][:, None]
            pressure = remove_depth_mean(pressure, vertical)
            velocity = remove_depth_mean(velocity, vertical)
            p_hat = fit_harmonic(pressure, post_times, args.period)
            u_hat = fit_harmonic(velocity, post_times, args.period)
            weights = vertical * horizontal[None, :]
            valid = weights > 0

            cross = np.sum(weights[valid] * np.real(p_hat[valid] * np.conj(u_hat[valid])))
            p2 = np.sum(weights[valid] * np.abs(p_hat[valid]) ** 2)
            u2 = np.sum(weights[valid] * np.abs(u_hat[valid]) ** 2)
            chi = cross / np.sqrt(p2 * u2) if p2 > 0 and u2 > 0 else np.nan
            flux = 0.5 * args.rho0 * cross
            if np.isfinite(chi) and chi > 0:
                energy_ratio = max(0.0, (1 - min(chi, 1.0)) / (1 + min(chi, 1.0)))
                amplitude_ratio = np.sqrt(energy_ratio)
            else:
                energy_ratio = np.nan
                amplitude_ratio = np.nan
            results.append({
                "boundary": side, "offset_cells": offset,
                "distance_km": distance / 1000, "outward_flux_W": flux,
                "progressive_index": chi,
                "reflection_energy_proxy": energy_ratio,
                "reflection_amplitude_proxy": amplitude_ratio,
            })

    eta_hat = fit_harmonic(eta_means[keep, None], post_times, args.period).item()
    print(f"Post-ramp domain-mean ETAN M2 amplitude: {abs(eta_hat):.6g} m")
    print(f"Post-ramp domain-mean ETAN M2 phase: "
          f"{np.mod(-np.angle(eta_hat, deg=True), 360):.3f} deg GMT")
    return results


def write_results(run_dir: Path, results: list[dict]) -> None:
    import os
    import tempfile

    cache_root = Path(tempfile.gettempdir()) / "mitgcm-plot-cache"
    os.environ.setdefault("MPLCONFIGDIR", str(cache_root / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_root))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    csv_path = run_dir / "boundary_reflection_summary.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for side, label in (("W", "West"), ("S", "South"), ("N", "North")):
        subset = [item for item in results if item["boundary"] == side]
        distance = [item["distance_km"] for item in subset]
        flux = [item["outward_flux_W"] / 1e6 for item in subset]
        reflection = [item["reflection_amplitude_proxy"] for item in subset]
        axes[0].plot(distance, flux, "o-", label=label)
        axes[1].plot(distance, reflection, "o-", label=label)
    axes[0].axhline(0, color="0.3", linewidth=0.8)
    axes[0].set(xlabel="Distance inward from boundary (km)",
                ylabel="Outward baroclinic M2 flux (MW)")
    axes[1].axhline(0.1, color="0.5", linestyle="--", linewidth=0.8,
                    label="10% amplitude")
    axes[1].set(xlabel="Distance inward from boundary (km)",
                ylabel="Reflection-amplitude proxy", ylim=(0, 1))
    for axis in axes:
        axis.grid(alpha=0.25)
        axis.legend()
    fig.tight_layout()
    figure_path = run_dir / "boundary_reflection_summary.png"
    fig.savefig(figure_path, dpi=180)
    print(f"Wrote {csv_path}")
    print(f"Wrote {figure_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--period", type=float, default=44712.0)
    parser.add_argument("--ramp-periods", type=float, default=3.0)
    parser.add_argument("--rho0", type=float, default=1026.0)
    parser.add_argument("--offsets", type=lambda value: [int(x) for x in value.split(",")],
                        # LR sponge occupies offsets 1:20. These sample
                        # approximately the same physical distances as HR
                        # old-domain offsets 5, 10, 20 and 40 cells.
                        default=[22, 23, 27, 33])
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    output = analyze(arguments)
    write_results(arguments.run_dir.resolve(), output)
