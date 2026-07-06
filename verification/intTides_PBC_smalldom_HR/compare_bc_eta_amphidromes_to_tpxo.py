#!/usr/bin/env python3
"""Compare boundary-test M2 SSH amphidrome structure against TPXO9.

The boundary tests mostly saved the light diagnostic stream
``bc_wave_eta.*``.  This script fits an M2 harmonic to those ETAN snapshots,
interpolates the TPXO9 M2 elevation harmonic onto the same MITgcm grid, and
evaluates whether each test develops a model SSH node/amphidrome similar to
the TPXO reference.

The key metric is not only "minimum |eta|" but also phase winding around that
minimum.  A robust amphidrome should have small |eta| and phase winding close
to +/- 2*pi around a closed ring.
"""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import xarray as xr
from scipy.interpolate import RegularGridInterpolator

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import cmocean

    CMAP_AMP = cmocean.cm.amp
    CMAP_BALANCE = cmocean.cm.balance
    CMAP_PHASE = cmocean.cm.phase
    CMAP_DELTA = cmocean.cm.delta
except Exception:  # pragma: no cover
    CMAP_AMP = "viridis"
    CMAP_BALANCE = "RdBu_r"
    CMAP_PHASE = "twilight"
    CMAP_DELTA = "PuOr_r"


@dataclass(frozen=True)
class Meta:
    nx: int
    ny: int
    nz: int
    dtype: np.dtype
    nrecords: int
    fields: list[str]
    time: float | None


@dataclass(frozen=True)
class Point:
    name: str
    j: int
    i: int


def parse_meta(path: Path) -> Meta:
    text = path.read_text(errors="ignore")
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
    precision = re.search(r"dataprec\s*=\s*\[\s*'([^']+)'", text).group(1).strip()
    dtype = np.dtype(">f4" if precision == "float32" else ">f8")
    nrecords_match = re.search(r"nrecords\s*=\s*\[\s*(\d+)", text)
    nrecords = int(nrecords_match.group(1)) if nrecords_match else 1
    field_block = re.search(r"fldList\s*=\s*\{(.*?)\};", text, re.S)
    fields = [field.strip() for field in re.findall(r"'([^']+)'", field_block.group(1))] if field_block else []
    time_match = re.search(r"timeInterval\s*=\s*\[\s*([^\]]+)", text)
    time = float(time_match.group(1)) if time_match else None
    return Meta(nx, ny, nz, dtype, nrecords, fields, time)


def read_mds_field(meta_path: Path) -> tuple[np.ndarray, Meta]:
    meta = parse_meta(meta_path)
    data_path = meta_path.with_suffix(".data")
    expected = meta.nrecords * meta.nz * meta.ny * meta.nx
    values = np.fromfile(data_path, dtype=meta.dtype)
    if values.size != expected:
        raise ValueError(f"{data_path}: expected {expected} values, got {values.size}")
    array = values.reshape(meta.nrecords, meta.nz, meta.ny, meta.nx).astype(np.float64)
    return array, meta


def read_grid(run_dir: Path, name: str) -> np.ndarray:
    array, meta = read_mds_field(run_dir / f"{name}.meta")
    if meta.nrecords != 1:
        raise ValueError(f"{name}: expected one record, got {meta.nrecords}")
    result = array[0]
    return result if meta.nz > 1 else result[0]


def collect_eta_files(run_dir: Path) -> list[Path]:
    files = sorted(p for p in run_dir.glob("bc_wave_eta.[0-9]*.meta") if not p.name.startswith("._"))
    if len(files) < 8:
        raise FileNotFoundError(f"{run_dir}: found only {len(files)} bc_wave_eta frames")
    return files


def crop_slices(nx: int, ny: int, sponge: int) -> tuple[slice, slice]:
    if sponge < 0:
        raise ValueError("sponge must be non-negative")
    if sponge == 0:
        return slice(0, ny), slice(0, nx)
    if 2 * sponge >= ny or sponge >= nx:
        raise ValueError(f"sponge={sponge} incompatible with nx={nx}, ny={ny}")
    return slice(sponge, ny - sponge), slice(sponge, nx)


def harmonic_weights(times: np.ndarray, period: float) -> tuple[np.ndarray, np.ndarray]:
    centered = (times - times.mean()) / period
    omega_t = 2 * np.pi * times / period
    design = np.column_stack((np.ones(times.size), centered, np.cos(omega_t), np.sin(omega_t)))
    pinv = np.linalg.pinv(design)
    return pinv[2], pinv[3]


def fit_eta_harmonic(files: list[Path], ysl: slice, xsl: slice, period: float, spinup_periods: float) -> tuple[np.ndarray, np.ndarray]:
    metas = [parse_meta(path) for path in files]
    times_all = np.asarray([meta.time for meta in metas], dtype=np.float64)
    keep = np.isfinite(times_all) & (times_all >= spinup_periods * period - 0.1)
    files = [path for path, ok in zip(files, keep) if ok]
    times = times_all[keep]
    if times.size < 8:
        raise ValueError("Fewer than eight post-spinup ETAN samples are available")
    cos_w, sin_w = harmonic_weights(times, period)

    shape = (ysl.stop - ysl.start, xsl.stop - xsl.start)
    eta_cos = np.zeros(shape, dtype=np.float64)
    eta_sin = np.zeros(shape, dtype=np.float64)
    for n, path in enumerate(files):
        array, meta = read_mds_field(path)
        if meta.fields and "ETAN" not in meta.fields:
            raise ValueError(f"{path}: expected ETAN, found {meta.fields}")
        eta = array[0, 0, ysl, xsl]
        eta_cos += cos_w[n] * eta
        eta_sin += sin_w[n] * eta
    return eta_cos - 1j * eta_sin, times


def constituent_index(dataset: xr.Dataset, name: str) -> int:
    names = [value.decode().strip().upper() for value in dataset.con.values]
    try:
        return names.index(name.upper())
    except ValueError as exc:
        raise ValueError(f"TPXO constituent {name!r} not found in {names}") from exc


def interpolate_tpxo_eta(tpxo_elevation: Path, constituent: str, lon_target: np.ndarray, lat_target: np.ndarray) -> np.ndarray:
    with xr.open_dataset(tpxo_elevation) as dataset:
        index = constituent_index(dataset, constituent)
        lon = np.asarray(dataset["lon_z"][:, 0], dtype=np.float64)
        lat = np.asarray(dataset["lat_z"][0, :], dtype=np.float64)
        target_x = np.mod(lon_target, 360.0)
        target_y = lat_target
        buffer = 0.5
        ix = np.flatnonzero((lon >= np.nanmin(target_x) - buffer) & (lon <= np.nanmax(target_x) + buffer))
        iy = np.flatnonzero((lat >= np.nanmin(target_y) - buffer) & (lat <= np.nanmax(target_y) + buffer))
        if ix.size < 2 or iy.size < 2:
            raise ValueError("TPXO subset does not bracket the model grid")
        amp = np.asarray(dataset["ha"].isel(nc=index, nx=slice(ix[0], ix[-1] + 1), ny=slice(iy[0], iy[-1] + 1)), dtype=np.float64)
        phase = np.asarray(dataset["hp"].isel(nc=index, nx=slice(ix[0], ix[-1] + 1), ny=slice(iy[0], iy[-1] + 1)), dtype=np.float64)
    harmonic = amp * np.exp(-1j * np.deg2rad(phase))
    interp = RegularGridInterpolator(
        (lon[ix], lat[iy]),
        harmonic,
        method="linear",
        bounds_error=False,
        fill_value=np.nan + 1j * np.nan,
    )
    return interp(np.column_stack((target_x.ravel(), target_y.ravel()))).reshape(target_x.shape)


def wrap_phase(angle: np.ndarray | float) -> np.ndarray | float:
    return (np.asarray(angle) + np.pi) % (2 * np.pi) - np.pi


def ring_indices(j: int, i: int, radius: int, shape: tuple[int, int]) -> tuple[np.ndarray, np.ndarray]:
    ny, nx = shape
    js: list[int] = []
    is_: list[int] = []
    for di in range(-radius, radius + 1):
        for dj in (-radius, radius):
            jj = j + dj
            ii = i + di
            if 0 <= jj < ny and 0 <= ii < nx:
                js.append(jj)
                is_.append(ii)
    for dj in range(-radius + 1, radius):
        for di in (radius, -radius):
            jj = j + dj
            ii = i + di
            if 0 <= jj < ny and 0 <= ii < nx:
                js.append(jj)
                is_.append(ii)
    return np.asarray(js, dtype=int), np.asarray(is_, dtype=int)


def phase_winding(phase: np.ndarray, point: Point, radius: int, wet: np.ndarray) -> float:
    js, is_ = ring_indices(point.j, point.i, radius, phase.shape)
    if js.size < 8:
        return np.nan
    values = phase[js, is_]
    ok = wet[js, is_] & np.isfinite(values)
    if np.count_nonzero(ok) < max(8, int(0.75 * values.size)):
        return np.nan
    values = values[ok]
    return float(np.nansum(wrap_phase(np.diff(np.r_[values, values[0]]))))


def nearest_index(xc: np.ndarray, yc: np.ndarray, lon: float, lat: float, wet: np.ndarray) -> tuple[int, int]:
    distance2 = (xc - lon) ** 2 + (yc - lat) ** 2
    distance2 = np.where(wet, distance2, np.inf)
    return tuple(int(v) for v in np.unravel_index(np.nanargmin(distance2), distance2.shape))


def min_amp_point(name: str, eta: np.ndarray, wet: np.ndarray, depth: np.ndarray, min_depth: float) -> Point:
    score = np.where(wet & (depth >= min_depth) & np.isfinite(eta), np.abs(eta), np.inf)
    if not np.isfinite(score).any():
        raise ValueError(f"No wet finite cells deeper than {min_depth:g} m")
    j, i = np.unravel_index(np.nanargmin(score), score.shape)
    return Point(name, int(j), int(i))


def km_distance(lon1: float, lat1: float, lon2: float, lat2: float) -> float:
    lat0 = np.deg2rad(0.5 * (lat1 + lat2))
    dx = (lon2 - lon1) * 111.0 * np.cos(lat0)
    dy = (lat2 - lat1) * 111.0
    return float(np.hypot(dx, dy))


def masked(field: np.ndarray, wet: np.ndarray) -> np.ma.MaskedArray:
    return np.ma.masked_where((~wet) | (~np.isfinite(field)), field)


def draw_markers(ax, points: list[Point], xc: np.ndarray, yc: np.ndarray) -> None:
    styles = {
        "model_min_eta_amp": dict(marker="*", s=155, c="white", edgecolors="k", linewidths=0.9, label="model min |η|"),
        "tpxo_min_eta_amp": dict(marker="P", s=105, c="gold", edgecolors="k", linewidths=0.8, label="TPXO min |η|"),
        "canyon_marker": dict(marker="X", s=90, c="magenta", edgecolors="k", linewidths=0.7, label="canyon"),
    }
    used: set[str] = set()
    for point in points:
        style = styles[point.name].copy()
        if style["label"] in used:
            style.pop("label")
        else:
            used.add(style["label"])
        ax.scatter(float(xc[point.j, point.i]), float(yc[point.j, point.i]), zorder=9, **style)


def plot_run(
    save_path: Path,
    run_name: str,
    xc: np.ndarray,
    yc: np.ndarray,
    depth: np.ndarray,
    wet: np.ndarray,
    model_eta: np.ndarray,
    tpxo_eta: np.ndarray,
    points: list[Point],
    args: argparse.Namespace,
) -> None:
    model_amp = np.abs(model_eta)
    tpxo_amp = np.abs(tpxo_eta)
    amp_diff = model_amp - tpxo_amp
    phase_diff = wrap_phase(np.angle(model_eta) - np.angle(tpxo_eta))
    model_amp_lim = np.nanpercentile(model_amp[wet], 99.0)
    tpxo_amp_lim = np.nanpercentile(tpxo_amp[wet], 99.0)
    amp_diff_lim = np.nanpercentile(np.abs(amp_diff[wet & np.isfinite(amp_diff)]), 98.0)

    fig, axes = plt.subplots(2, 3, figsize=(15.8, 8.8), constrained_layout=True)
    for ax in axes.ravel():
        ax.set_facecolor(args.land_color)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("Longitude")
        ax.contour(xc, yc, masked(depth, wet), levels=args.bathy_contours, colors="k", linewidths=0.25, alpha=0.35)
    for ax in axes[:, 0]:
        ax.set_ylabel("Latitude")
    fields = [
        (model_amp, "Model |η|", CMAP_AMP, 0, model_amp_lim, "m"),
        (tpxo_amp, "TPXO |η|", CMAP_AMP, 0, tpxo_amp_lim, "m"),
        (amp_diff, "Model − TPXO |η|", CMAP_BALANCE, -amp_diff_lim, amp_diff_lim, "m"),
        (np.angle(model_eta), "Model η phase", CMAP_PHASE, -np.pi, np.pi, "rad"),
        (np.angle(tpxo_eta), "TPXO η phase", CMAP_PHASE, -np.pi, np.pi, "rad"),
        (phase_diff, "Model − TPXO η phase", CMAP_DELTA, -np.pi, np.pi, "rad"),
    ]
    for ax, (field, title, cmap, vmin, vmax, label) in zip(axes.ravel(), fields):
        m = ax.pcolormesh(xc, yc, masked(field, wet), shading="auto", cmap=cmap, vmin=vmin, vmax=vmax)
        draw_markers(ax, points, xc, yc)
        ax.set_title(title)
        fig.colorbar(m, ax=ax, shrink=0.84, label=label)
    axes[0, 2].legend(loc="upper right", frameon=True, fontsize=8)
    fig.suptitle(f"{run_name}: M2 ETAN amphidrome comparison against TPXO9")
    fig.savefig(save_path, dpi=args.dpi)
    plt.close(fig)


def analyze_run(run_dir: Path, save_dir: Path, args: argparse.Namespace) -> dict[str, float | str]:
    eta_files = collect_eta_files(run_dir)
    meta0 = parse_meta(eta_files[0])
    ysl, xsl = crop_slices(meta0.nx, meta0.ny, 0 if args.include_sponge else args.sponge_cells)
    xc = read_grid(run_dir, "XC")[ysl, xsl]
    yc = read_grid(run_dir, "YC")[ysl, xsl]
    depth = read_grid(run_dir, "Depth")[ysl, xsl]
    wet = read_grid(run_dir, "hFacC")[0, ysl, xsl] > 0

    model_eta, times = fit_eta_harmonic(eta_files, ysl, xsl, args.period, args.spinup_periods)
    tpxo_eta = interpolate_tpxo_eta(args.tpxo_elevation, args.constituent, xc, yc)
    wet = wet & np.isfinite(model_eta) & np.isfinite(tpxo_eta)

    model_min = min_amp_point("model_min_eta_amp", model_eta, wet, depth, args.min_depth_for_amphidrome)
    tpxo_min = min_amp_point("tpxo_min_eta_amp", tpxo_eta, wet, depth, args.min_depth_for_amphidrome)
    cj, ci = nearest_index(xc, yc, args.canyon_lon, args.canyon_lat, wet)
    canyon = Point("canyon_marker", cj, ci)
    points = [model_min, tpxo_min, canyon]

    phase_model = np.angle(model_eta)
    phase_tpxo = np.angle(tpxo_eta)
    finite = wet & np.isfinite(model_eta) & np.isfinite(tpxo_eta)
    row: dict[str, float | str] = {
        "run": run_dir.name,
        "samples_used": int(times.size),
        "time_start_s": float(times.min()),
        "time_end_s": float(times.max()),
        "model_eta_amp_p50_m": float(np.nanpercentile(np.abs(model_eta[finite]), 50.0)),
        "model_eta_amp_p95_m": float(np.nanpercentile(np.abs(model_eta[finite]), 95.0)),
        "model_eta_amp_max_m": float(np.nanmax(np.abs(model_eta[finite]))),
        "tpxo_eta_amp_p50_m": float(np.nanpercentile(np.abs(tpxo_eta[finite]), 50.0)),
        "model_min_lon": float(xc[model_min.j, model_min.i]),
        "model_min_lat": float(yc[model_min.j, model_min.i]),
        "model_min_depth_m": float(depth[model_min.j, model_min.i]),
        "model_min_eta_amp_m": float(np.abs(model_eta[model_min.j, model_min.i])),
        "model_min_eta_phase_rad": float(phase_model[model_min.j, model_min.i]),
        "model_min_winding_r4_rad": phase_winding(phase_model, model_min, 4, wet),
        "model_min_winding_r8_rad": phase_winding(phase_model, model_min, 8, wet),
        "tpxo_at_model_min_eta_amp_m": float(np.abs(tpxo_eta[model_min.j, model_min.i])),
        "tpxo_at_model_min_eta_phase_rad": float(phase_tpxo[model_min.j, model_min.i]),
        "tpxo_at_model_min_winding_r8_rad": phase_winding(phase_tpxo, model_min, 8, wet),
        "tpxo_min_lon": float(xc[tpxo_min.j, tpxo_min.i]),
        "tpxo_min_lat": float(yc[tpxo_min.j, tpxo_min.i]),
        "tpxo_min_depth_m": float(depth[tpxo_min.j, tpxo_min.i]),
        "tpxo_min_eta_amp_m": float(np.abs(tpxo_eta[tpxo_min.j, tpxo_min.i])),
        "model_at_tpxo_min_eta_amp_m": float(np.abs(model_eta[tpxo_min.j, tpxo_min.i])),
        "model_to_tpxo_min_distance_km": km_distance(xc[model_min.j, model_min.i], yc[model_min.j, model_min.i], xc[tpxo_min.j, tpxo_min.i], yc[tpxo_min.j, tpxo_min.i]),
        "model_min_to_canyon_km": km_distance(xc[model_min.j, model_min.i], yc[model_min.j, model_min.i], xc[canyon.j, canyon.i], yc[canyon.j, canyon.i]),
        "canyon_model_eta_amp_m": float(np.abs(model_eta[canyon.j, canyon.i])),
        "canyon_tpxo_eta_amp_m": float(np.abs(tpxo_eta[canyon.j, canyon.i])),
        "canyon_phase_diff_rad": float(wrap_phase(phase_model[canyon.j, canyon.i] - phase_tpxo[canyon.j, canyon.i])),
        "median_abs_eta_amp_diff_m": float(np.nanmedian(np.abs(np.abs(model_eta[finite]) - np.abs(tpxo_eta[finite])))),
        "median_abs_phase_diff_rad": float(np.nanmedian(np.abs(wrap_phase(phase_model[finite] - phase_tpxo[finite])))),
    }

    plot_run(save_dir / f"{run_dir.name}_m2_eta_tpxo_comparison.png", run_dir.name, xc, yc, depth, wet, model_eta, tpxo_eta, points, args)
    return row


def write_summary(path: Path, rows: list[dict[str, float | str]], args: argparse.Namespace) -> None:
    lines = [
        "# Boundary-test M2 SSH amphidrome comparison to TPXO9",
        "",
        f"TPXO elevation file: `{args.tpxo_elevation}`",
        f"Constituent: {args.constituent}",
        f"Spin-up excluded before harmonic fit: {args.spinup_periods:g} M2 periods",
        f"Sponge crop W/S/N: {0 if args.include_sponge else args.sponge_cells} cells",
        f"Minimum depth for min-|η| search: {args.min_depth_for_amphidrome:g} m",
        "",
        "## Key interpretation columns",
        "",
        "- `model_min_eta_amp_m`: how close the model gets to an SSH node.",
        "- `model_min_winding_r8_rad`: values near ±6.283 imply amphidromic phase winding.",
        "- `tpxo_at_model_min_eta_amp_m` and `tpxo_at_model_min_winding_r8_rad`: whether TPXO has a comparable feature at the model node.",
        "- `model_to_tpxo_min_distance_km`: separation between model and TPXO low-|η| candidates.",
        "",
        "## Compact result",
        "",
        "| run | model |η| p95 (m) | model min lon,lat | model |η|min (m) | model r8 winding (rad) | TPXO |η| at model min (m) | TPXO r8 winding at model min | distance to TPXO min (km) | canyon model/TPXO |η| (m) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['run']} | {row['model_eta_amp_p95_m']:.4g} | "
            f"{row['model_min_lon']:.4f}, {row['model_min_lat']:.4f} | "
            f"{row['model_min_eta_amp_m']:.4g} | {row['model_min_winding_r8_rad']:.4g} | "
            f"{row['tpxo_at_model_min_eta_amp_m']:.4g} | {row['tpxo_at_model_min_winding_r8_rad']:.4g} | "
            f"{row['model_to_tpxo_min_distance_km']:.2f} | "
            f"{row['canyon_model_eta_amp_m']:.4g}/{row['canyon_tpxo_eta_amp_m']:.4g} |"
        )
    lines.extend(
        [
            "",
            "## Reading the result",
            "",
            "If a TPXO-forced test has a near-zero model |η| with ~2π winding, but TPXO has large |η| and no winding at the same point, then the node is being produced by the model configuration rather than inherited from TPXO.",
            "",
        ]
    )
    path.write_text("\n".join(lines))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "runs",
        nargs="*",
        type=Path,
        help="Boundary-test run directories. Defaults to available LR fullTPXO bc_wave_eta runs.",
    )
    parser.add_argument("--save-dir", type=Path, default=Path("intTides_PBC_smalldom_HR/bc_eta_tpxo_amphidrome_comparison"))
    parser.add_argument("--tpxo-elevation", type=Path, default=Path("/Users/karina/Research/Canyons/data/tides/h_tpxo9.v5a.nc"))
    parser.add_argument("--constituent", default="M2")
    parser.add_argument("--period", type=float, default=44712.0)
    parser.add_argument("--spinup-periods", type=float, default=1.0)
    parser.add_argument("--sponge-cells", type=int, default=20)
    parser.add_argument("--include-sponge", action="store_true")
    parser.add_argument("--min-depth-for-amphidrome", type=float, default=50.0)
    parser.add_argument("--canyon-lon", type=float, default=243.25)
    parser.add_argument("--canyon-lat", type=float, default=31.785)
    parser.add_argument("--land-color", default="tan")
    parser.add_argument("--dpi", type=int, default=160)
    parser.add_argument("--bathy-contours", nargs="+", type=float, default=[50, 100, 200, 500, 1000, 1500, 2000])
    return parser.parse_args()


def default_runs() -> list[Path]:
    root = Path("intTides_PBC_smalldom_HR")
    names = [
        "runBC_LR_fullTPXO_balanced",
        "runBC_LR_fullTPXO_obeta",
        "runBC_LR_fullTPXO_obeta_balanced",
        "runBC_LR_fullTPXO_obeta_weakvel025",
        "runBC_LR_fullTPXO_obeta_weakvel010",
        "runBC_LR_fullTPXO_obeta_zerovel",
    ]
    return [root / name for name in names if (root / name).is_dir()]


def main() -> None:
    args = parse_args()
    runs = args.runs or default_runs()
    args.save_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, float | str]] = []
    for run in runs:
        run = run.resolve()
        print(f"Analyzing {run}")
        try:
            rows.append(analyze_run(run, args.save_dir.resolve(), args))
        except Exception as exc:
            print(f"WARNING: skipped {run}: {exc}")
    if not rows:
        raise RuntimeError("No runs were successfully analyzed")

    csv_path = args.save_dir / "bc_eta_tpxo_amphidrome_summary.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    md_path = args.save_dir / "bc_eta_tpxo_amphidrome_summary.md"
    write_summary(md_path, rows, args)
    print(f"Wrote {csv_path}")
    print(f"Wrote {md_path}")


if __name__ == "__main__":
    main()
