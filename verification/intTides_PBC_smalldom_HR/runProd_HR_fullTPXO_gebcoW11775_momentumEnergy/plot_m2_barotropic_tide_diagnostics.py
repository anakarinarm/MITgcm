#!/usr/bin/env python3
"""Plot M2 free-surface and barotropic-tide diagnostics.

This script addresses the "is there an amphidrome and does it matter for the
canyon?" question.  It fits M2 harmonics after the spin-up period and produces:

1. ETAN amplitude and phase maps;
2. barotropic velocity amplitude and phase-vector maps;
3. a bathymetry/context map marking the amphidromic candidate and a canyon
   marker;
4. a small text summary with coordinates and local amplitudes.

Land is masked, so the Matplotlib axes facecolor shows through.  Use
``--land-color tan`` or ``--land-color gray``.
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import cmocean

    CMAP_AMP = cmocean.cm.amp
    CMAP_PHASE = cmocean.cm.phase
    CMAP_SPEED = cmocean.cm.speed
    CMAP_DEEP = cmocean.cm.deep
except Exception:  # pragma: no cover
    CMAP_AMP = "viridis"
    CMAP_PHASE = "twilight"
    CMAP_SPEED = "viridis"
    CMAP_DEEP = "Blues"


@dataclass
class Meta:
    nx: int
    ny: int
    nz: int
    dtype: np.dtype
    nrecords: int
    fields: list[str]
    time: float | None


@dataclass
class FramePair:
    meta_2d: Path
    data_2d: Path
    meta_3d: Path
    data_3d: Path
    time_s: float


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
    fields = re.findall(r"'([^']+)'", field_block.group(1)) if field_block else []
    fields = [field.strip() for field in fields]
    time_match = re.search(r"timeInterval\s*=\s*\[\s*([^\]]+)", text)
    time = float(time_match.group(1)) if time_match else None
    return Meta(nx, ny, nz, dtype, nrecords, fields, time)


def memmap_mds(meta_path: Path) -> tuple[np.memmap, Meta]:
    meta = parse_meta(meta_path)
    data_path = meta_path.with_suffix(".data")
    expected = meta.nrecords * meta.nz * meta.ny * meta.nx
    actual = data_path.stat().st_size // meta.dtype.itemsize
    if actual != expected:
        raise ValueError(f"{data_path}: expected {expected} values, got {actual}")
    array = np.memmap(
        data_path,
        dtype=meta.dtype,
        mode="r",
        shape=(meta.nrecords, meta.nz, meta.ny, meta.nx),
    )
    return array, meta


def read_grid(run_dir: Path, name: str) -> np.ndarray:
    array, meta = memmap_mds(run_dir / f"{name}.meta")
    if meta.nrecords != 1:
        raise ValueError(f"{name}: expected one record, got {meta.nrecords}")
    result = np.asarray(array[0], dtype=np.float64)
    return result if meta.nz > 1 else result[0]


def field_index(meta: Meta, field: str) -> int:
    if field not in meta.fields:
        raise ValueError(f"{field} absent from {meta.fields}")
    return meta.fields.index(field)


def collect_frame_pairs(diag_dir: Path) -> list[FramePair]:
    files_2d = {
        p.name.split(".")[1]: p
        for p in diag_dir.glob("diag_state_2d.[0-9]*.meta")
        if not p.name.startswith("._")
    }
    files_3d = {
        p.name.split(".")[1]: p
        for p in diag_dir.glob("diag_state_3d.[0-9]*.meta")
        if not p.name.startswith("._")
    }
    keys = sorted(set(files_2d).intersection(files_3d))
    if not keys:
        raise FileNotFoundError(f"No matching diag_state_2d/3d files in {diag_dir}")

    pairs: list[FramePair] = []
    for key in keys:
        meta_2d = files_2d[key]
        meta_3d = files_3d[key]
        m2 = parse_meta(meta_2d)
        m3 = parse_meta(meta_3d)
        if m2.time is None or m3.time is None:
            raise ValueError(f"Missing timeInterval in {meta_2d} or {meta_3d}")
        if abs(m2.time - m3.time) > 0.1:
            raise ValueError(f"2-D/3-D times differ for iteration {key}")
        pairs.append(FramePair(meta_2d, meta_2d.with_suffix(".data"), meta_3d, meta_3d.with_suffix(".data"), m2.time))
    return pairs


def crop_slices(nx: int, ny: int, sponge: int) -> tuple[slice, slice]:
    if sponge < 0:
        raise ValueError("sponge must be non-negative")
    if 2 * sponge >= ny or sponge >= nx:
        raise ValueError(f"sponge={sponge} incompatible with nx={nx}, ny={ny}")
    return slice(sponge, ny - sponge), slice(sponge, nx)


def harmonic_weights(times: np.ndarray, period: float) -> tuple[np.ndarray, np.ndarray]:
    centered = (times - times.mean()) / period
    omega_t = 2 * np.pi * times / period
    design = np.column_stack((np.ones(times.size), centered, np.cos(omega_t), np.sin(omega_t)))
    pinv = np.linalg.pinv(design)
    return pinv[2], pinv[3]


def depth_mean_velocity_to_c(
    u_mass: np.ndarray,
    v_mass: np.ndarray,
    drf: np.ndarray,
    hfac_w: np.ndarray,
    hfac_s: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    zw = hfac_w * drf[:, None, None]
    zs = hfac_s * drf[:, None, None]
    hw = np.sum(zw, axis=0)
    hs = np.sum(zs, axis=0)
    u_face = np.divide(np.sum(u_mass * zw, axis=0), hw, out=np.zeros_like(hw), where=hw > 0)
    v_face = np.divide(np.sum(v_mass * zs, axis=0), hs, out=np.zeros_like(hs), where=hs > 0)
    u_c = np.empty_like(u_face)
    v_c = np.empty_like(v_face)
    u_c[:, :-1] = 0.5 * (u_face[:, :-1] + u_face[:, 1:])
    u_c[:, -1] = u_face[:, -1]
    v_c[:-1, :] = 0.5 * (v_face[:-1, :] + v_face[1:, :])
    v_c[-1, :] = v_face[-1, :]
    return u_c, v_c


def fit_tide(
    pairs: list[FramePair],
    run_dir: Path,
    diag_dir: Path,
    ysl: slice,
    xsl: slice,
    period: float,
    spinup_periods: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    times_all = np.asarray([p.time_s for p in pairs], dtype=np.float64)
    keep = times_all >= spinup_periods * period - 0.1
    pairs = [p for p, ok in zip(pairs, keep) if ok]
    times = np.asarray([p.time_s for p in pairs], dtype=np.float64)
    if times.size < 8:
        raise ValueError("Fewer than eight post-spinup samples are available")
    cos_w, sin_w = harmonic_weights(times, period)

    drf = read_grid(run_dir, "DRF").reshape(-1)
    hfac_w = read_grid(run_dir, "hFacW")
    hfac_s = read_grid(run_dir, "hFacS")

    meta_2d = parse_meta(pairs[0].meta_2d)
    meta_3d = parse_meta(pairs[0].meta_3d)
    i_eta = field_index(meta_2d, "ETAN")
    i_u = field_index(meta_3d, "UVELMASS")
    i_v = field_index(meta_3d, "VVELMASS")

    shape = (ysl.stop - ysl.start, xsl.stop - xsl.start)
    eta_cos = np.zeros(shape, dtype=np.float64)
    eta_sin = np.zeros(shape, dtype=np.float64)
    u_cos = np.zeros(shape, dtype=np.float64)
    u_sin = np.zeros(shape, dtype=np.float64)
    v_cos = np.zeros(shape, dtype=np.float64)
    v_sin = np.zeros(shape, dtype=np.float64)

    for n, pair in enumerate(pairs):
        arr2, meta2 = memmap_mds(pair.meta_2d)
        arr3, meta3 = memmap_mds(pair.meta_3d)
        eta = np.asarray(arr2[i_eta, 0, ysl, xsl], dtype=np.float64)
        u_mass = np.asarray(arr3[i_u], dtype=np.float64)
        v_mass = np.asarray(arr3[i_v], dtype=np.float64)
        u_bt, v_bt = depth_mean_velocity_to_c(u_mass, v_mass, drf, hfac_w, hfac_s)
        u_bt = u_bt[ysl, xsl]
        v_bt = v_bt[ysl, xsl]

        eta_cos += cos_w[n] * eta
        eta_sin += sin_w[n] * eta
        u_cos += cos_w[n] * u_bt
        u_sin += sin_w[n] * u_bt
        v_cos += cos_w[n] * v_bt
        v_sin += sin_w[n] * v_bt

    eta_hat = eta_cos - 1j * eta_sin
    u_hat = u_cos - 1j * u_sin
    v_hat = v_cos - 1j * v_sin
    return eta_hat, u_hat, v_hat, times


def central_gradient_metric(field: np.ndarray, dx: np.ndarray, dy: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    dfdx = np.full_like(field, np.nan, dtype=np.float64)
    dfdy = np.full_like(field, np.nan, dtype=np.float64)
    denom_x = dx[:, 1:-1] + 0.5 * (dx[:, :-2] + dx[:, 2:])
    denom_y = dy[1:-1, :] + 0.5 * (dy[:-2, :] + dy[2:, :])
    dfdx[:, 1:-1] = (field[:, 2:] - field[:, :-2]) / denom_x
    dfdy[1:-1, :] = (field[2:, :] - field[:-2, :]) / denom_y
    dfdx[:, 0] = (field[:, 1] - field[:, 0]) / (0.5 * (dx[:, 0] + dx[:, 1]))
    dfdx[:, -1] = (field[:, -1] - field[:, -2]) / (0.5 * (dx[:, -1] + dx[:, -2]))
    dfdy[0, :] = (field[1, :] - field[0, :]) / (0.5 * (dy[0, :] + dy[1, :]))
    dfdy[-1, :] = (field[-1, :] - field[-2, :]) / (0.5 * (dy[-1, :] + dy[-2, :]))
    return dfdx, dfdy


def nearest_index(xc: np.ndarray, yc: np.ndarray, lon: float, lat: float, wet: np.ndarray) -> tuple[int, int]:
    distance2 = (xc - lon) ** 2 + (yc - lat) ** 2
    distance2 = np.where(wet, distance2, np.inf)
    return tuple(int(i) for i in np.unravel_index(np.nanargmin(distance2), distance2.shape))


def diagnostic_points(
    eta_amp: np.ndarray,
    speed_amp: np.ndarray,
    depth: np.ndarray,
    slope: np.ndarray,
    wet: np.ndarray,
    xc: np.ndarray,
    yc: np.ndarray,
    args: argparse.Namespace,
) -> dict[str, tuple[int, int]]:
    # Amphidrome candidate: minimum ETAN amplitude in wet cells, but avoid
    # cells where the fit is numerically zero everywhere.
    amp_for_min = np.where(wet & np.isfinite(eta_amp), eta_amp, np.inf)
    amph = tuple(int(i) for i in np.unravel_index(np.nanargmin(amp_for_min), amp_for_min.shape))

    if args.canyon_lon is not None and args.canyon_lat is not None:
        canyon = nearest_index(xc, yc, args.canyon_lon, args.canyon_lat, wet)
        canyon_key = "user_canyon_marker"
    else:
        # Automatic canyon marker: steepest wet topography in the physical plot.
        slope_for_max = np.where(wet & np.isfinite(slope), slope, -np.inf)
        canyon = tuple(int(i) for i in np.unravel_index(np.nanargmax(slope_for_max), slope_for_max.shape))
        canyon_key = "auto_steepest_slope_marker"

    return {"amphidrome_candidate": amph, canyon_key: canyon}


def km_distance(lon1: float, lat1: float, lon2: float, lat2: float) -> float:
    lat0 = np.deg2rad(0.5 * (lat1 + lat2))
    dx = (lon2 - lon1) * 111.0 * np.cos(lat0)
    dy = (lat2 - lat1) * 111.0
    return float(np.hypot(dx, dy))


def mask(array: np.ndarray, wet: np.ndarray) -> np.ma.MaskedArray:
    return np.ma.masked_where((~wet) | (~np.isfinite(array)), array)


def draw_markers(ax, points: dict[str, tuple[int, int]], xc: np.ndarray, yc: np.ndarray) -> None:
    styles = {
        "amphidrome_candidate": dict(marker="*", s=190, c="white", edgecolors="k", linewidths=1.1, label="ETAN amp min"),
        "user_canyon_marker": dict(marker="X", s=110, c="magenta", edgecolors="k", linewidths=0.8, label="user canyon marker"),
        "auto_steepest_slope_marker": dict(marker="X", s=100, c="magenta", edgecolors="k", linewidths=0.8, label="auto steepest slope"),
    }
    for name, (j, i) in points.items():
        ax.scatter(float(xc[j, i]), float(yc[j, i]), **styles[name], zorder=8)


def draw_physical_edges(ax, run_dir: Path, cells: int) -> None:
    if cells <= 0:
        return
    full_xc = read_grid(run_dir, "XC")
    full_yc = read_grid(run_dir, "YC")
    ax.axvline(float(full_xc[0, cells]), color="k", lw=1.0, ls="--", alpha=0.7)
    ax.axhline(float(full_yc[cells, 0]), color="k", lw=1.0, ls="--", alpha=0.7)
    ax.axhline(float(full_yc[-cells, 0]), color="k", lw=1.0, ls="--", alpha=0.7)


def plot_diagnostics(
    run_dir: Path,
    save_dir: Path,
    xc: np.ndarray,
    yc: np.ndarray,
    wet: np.ndarray,
    depth: np.ndarray,
    eta_hat: np.ndarray,
    u_hat: np.ndarray,
    v_hat: np.ndarray,
    slope: np.ndarray,
    points: dict[str, tuple[int, int]],
    args: argparse.Namespace,
) -> list[Path]:
    save_dir.mkdir(parents=True, exist_ok=True)
    eta_amp = np.abs(eta_hat)
    eta_phase = np.angle(eta_hat)
    speed_amp = np.hypot(np.abs(u_hat), np.abs(v_hat))
    # Phase-vector at t=0: useful for seeing propagation/rotation.
    u_phase0 = np.real(u_hat)
    v_phase0 = np.real(v_hat)

    outputs: list[Path] = []

    fig, axes = plt.subplots(1, 2, figsize=(13.4, 6.2), constrained_layout=True)
    for ax in axes:
        ax.set_facecolor(args.land_color)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("Longitude")
        if args.mark_sponge_boundary:
            draw_physical_edges(ax, run_dir, args.sponge_boundary_cells)
    axes[0].set_ylabel("Latitude")

    amp_lim = np.nanpercentile(eta_amp[wet], 99.0)
    m0 = axes[0].pcolormesh(xc, yc, mask(eta_amp, wet), shading="auto", cmap=CMAP_AMP, vmin=0, vmax=amp_lim)
    axes[0].contour(xc, yc, mask(depth, wet), levels=args.bathy_contours, colors="k", linewidths=0.35, alpha=0.45)
    axes[0].set_title("M2 ETAN amplitude")
    draw_markers(axes[0], points, xc, yc)
    fig.colorbar(m0, ax=axes[0], shrink=0.86, label="m")

    m1 = axes[1].pcolormesh(xc, yc, mask(eta_phase, wet), shading="auto", cmap=CMAP_PHASE, vmin=-np.pi, vmax=np.pi)
    axes[1].contour(xc, yc, mask(depth, wet), levels=args.bathy_contours, colors="k", linewidths=0.35, alpha=0.45)
    axes[1].set_title("M2 ETAN phase")
    draw_markers(axes[1], points, xc, yc)
    axes[1].legend(loc="upper right", frameon=True, fontsize=8)
    cb = fig.colorbar(m1, ax=axes[1], shrink=0.86, label="radians")
    cb.set_ticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi])
    cb.set_ticklabels(["$-\\pi$", "$-\\pi/2$", "0", "$\\pi/2$", "$\\pi$"])
    fig.suptitle("M2 free-surface amphidrome diagnostic")
    out = save_dir / "m2_etan_amplitude_phase.png"
    fig.savefig(out, dpi=args.dpi)
    plt.close(fig)
    outputs.append(out)

    fig, axes = plt.subplots(1, 2, figsize=(13.4, 6.2), constrained_layout=True)
    for ax in axes:
        ax.set_facecolor(args.land_color)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("Longitude")
        if args.mark_sponge_boundary:
            draw_physical_edges(ax, run_dir, args.sponge_boundary_cells)
    axes[0].set_ylabel("Latitude")

    speed_lim = np.nanpercentile(speed_amp[wet], 98.5)
    m0 = axes[0].pcolormesh(xc, yc, mask(speed_amp, wet), shading="auto", cmap=CMAP_SPEED, vmin=0, vmax=speed_lim)
    axes[0].contour(xc, yc, mask(depth, wet), levels=args.bathy_contours, colors="k", linewidths=0.35, alpha=0.45)
    axes[0].set_title("M2 barotropic speed amplitude")
    draw_markers(axes[0], points, xc, yc)
    fig.colorbar(m0, ax=axes[0], shrink=0.86, label="m s$^{-1}$")

    m1 = axes[1].pcolormesh(xc, yc, mask(speed_amp, wet), shading="auto", cmap=CMAP_SPEED, vmin=0, vmax=speed_lim)
    stride = args.quiver_stride
    qmask = wet[::stride, ::stride]
    qx = xc[::stride, ::stride]
    qy = yc[::stride, ::stride]
    qu = np.where(qmask, u_phase0[::stride, ::stride], np.nan)
    qv = np.where(qmask, v_phase0[::stride, ::stride], np.nan)
    axes[1].quiver(qx, qy, qu, qv, color="k", alpha=0.65, pivot="middle", width=0.0021)
    axes[1].contour(xc, yc, mask(depth, wet), levels=args.bathy_contours, colors="k", linewidths=0.35, alpha=0.45)
    axes[1].set_title("M2 BT phase vector, Re(U,V) at t=0")
    draw_markers(axes[1], points, xc, yc)
    fig.colorbar(m1, ax=axes[1], shrink=0.86, label="m s$^{-1}$")
    fig.suptitle("M2 barotropic tide over canyon/topography")
    out = save_dir / "m2_barotropic_velocity_amplitude_phase_vectors.png"
    fig.savefig(out, dpi=args.dpi)
    plt.close(fig)
    outputs.append(out)

    fig, ax = plt.subplots(1, 1, figsize=(7.8, 7.0), constrained_layout=True)
    ax.set_facecolor(args.land_color)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    if args.mark_sponge_boundary:
        draw_physical_edges(ax, run_dir, args.sponge_boundary_cells)
    max_depth = np.nanpercentile(depth[wet], 99.5)
    m = ax.pcolormesh(xc, yc, mask(depth, wet), shading="auto", cmap=CMAP_DEEP, vmin=0, vmax=max_depth)
    ax.contour(xc, yc, mask(slope, wet), levels=8, colors="0.15", linewidths=0.35, alpha=0.45)
    draw_markers(ax, points, xc, yc)
    ax.legend(loc="upper right", frameon=True, fontsize=8)
    ax.set_title("Bathymetry and diagnostic markers")
    fig.colorbar(m, ax=ax, shrink=0.88, label="depth m")
    out = save_dir / "m2_bathymetry_amphidrome_canyon_context.png"
    fig.savefig(out, dpi=args.dpi)
    plt.close(fig)
    outputs.append(out)

    return outputs


def write_summary(
    path: Path,
    xc: np.ndarray,
    yc: np.ndarray,
    depth: np.ndarray,
    slope: np.ndarray,
    eta_hat: np.ndarray,
    u_hat: np.ndarray,
    v_hat: np.ndarray,
    points: dict[str, tuple[int, int]],
    times: np.ndarray,
    args: argparse.Namespace,
) -> None:
    eta_amp = np.abs(eta_hat)
    eta_phase = np.angle(eta_hat)
    speed_amp = np.hypot(np.abs(u_hat), np.abs(v_hat))
    lines = [
        "# M2 barotropic tide diagnostic summary",
        "",
        f"Samples used: {times.size}",
        f"Time range: {times.min():.1f} to {times.max():.1f} s",
        f"Spin-up excluded: {args.spinup_periods:g} M2",
        f"Land/axis facecolor: {args.land_color}",
        f"Sponge cells cropped: {args.sponge_cells}",
        "",
    ]
    for name, (j, i) in points.items():
        lines.extend(
            [
                f"## {name}",
                "",
                f"- lon: {xc[j, i]:.8f}",
                f"- lat: {yc[j, i]:.8f}",
                f"- depth: {depth[j, i]:.3f} m",
                f"- slope proxy |grad H|: {slope[j, i]:.6g}",
                f"- ETAN M2 amplitude: {eta_amp[j, i]:.6g} m",
                f"- ETAN M2 phase: {eta_phase[j, i]:.6g} rad",
                f"- barotropic speed amplitude proxy: {speed_amp[j, i]:.6g} m/s",
                "",
            ]
        )
    a = points["amphidrome_candidate"]
    canyon_key = "user_canyon_marker" if "user_canyon_marker" in points else "auto_steepest_slope_marker"
    c = points[canyon_key]
    lines.extend(
        [
            "## Separation",
            "",
            "- approximate great-circle/planar distance between amphidrome candidate "
            f"and {canyon_key}: {km_distance(xc[a], yc[a], xc[c], yc[c]):.3f} km",
            "",
            "## Notes",
            "",
            "- The amphidrome candidate is the minimum wet-cell fitted M2 ETAN amplitude",
            "  in the plotted/cropped domain.",
            "- The canyon marker is user-specified if --canyon-lon/--canyon-lat are given;",
            "  otherwise it is the maximum wet-cell bathymetric slope proxy in the plotted domain.",
            "- Barotropic speed amplitude here is sqrt(|U_bt_hat|^2 + |V_bt_hat|^2).",
            "- For a quantitative TPXO comparison, provide an interior TPXO M2 eta/u/v",
            "  reference on this grid; this script currently diagnoses the MITgcm response.",
            "",
        ]
    )
    path.write_text("\n".join(lines))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--diag-dir", type=Path, default=None)
    parser.add_argument("--save-dir", type=Path, default=None)
    parser.add_argument("--period", type=float, default=44712.0)
    parser.add_argument("--spinup-periods", type=float, default=3.0)
    parser.add_argument("--sponge-cells", type=int, default=60)
    parser.add_argument("--sponge-boundary-cells", type=int, default=60)
    parser.add_argument("--include-sponge", action="store_true")
    parser.add_argument("--mark-sponge-boundary", action="store_true")
    parser.add_argument("--land-color", default="tan", help="Matplotlib axes facecolor for land, e.g. tan or gray.")
    parser.add_argument("--quiver-stride", type=int, default=24)
    parser.add_argument("--dpi", type=int, default=170)
    parser.add_argument("--canyon-lon", type=float, default=None)
    parser.add_argument("--canyon-lat", type=float, default=None)
    parser.add_argument(
        "--bathy-contours",
        nargs="+",
        type=float,
        default=[50, 100, 200, 500, 1000, 1500, 2000],
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_dir = args.run_dir.resolve()
    diag_dir = (args.diag_dir or (run_dir / "output")).resolve()
    save_dir = (args.save_dir or (diag_dir / "m2_barotropic_tide_diagnostics")).resolve()
    pairs = collect_frame_pairs(diag_dir)
    meta_2d = parse_meta(pairs[0].meta_2d)
    sponge = 0 if args.include_sponge else args.sponge_cells
    ysl, xsl = crop_slices(meta_2d.nx, meta_2d.ny, sponge)

    xc = read_grid(run_dir, "XC")[ysl, xsl]
    yc = read_grid(run_dir, "YC")[ysl, xsl]
    depth = read_grid(run_dir, "Depth")[ysl, xsl]
    hfac0 = read_grid(run_dir, "hFacC")[0, ysl, xsl]
    wet = hfac0 > 0
    dxc = read_grid(run_dir, "DXC")[ysl, xsl]
    dyc = read_grid(run_dir, "DYC")[ysl, xsl]
    dHdx, dHdy = central_gradient_metric(depth, dxc, dyc)
    slope = np.hypot(dHdx, dHdy)

    print(f"Run directory: {run_dir}")
    print(f"Diagnostics:   {diag_dir}")
    print(f"Output:        {save_dir}")
    print(f"Frames found:  {len(pairs)}")
    print(f"Crop:          j={ysl.start}:{ysl.stop}, i={xsl.start}:{xsl.stop}")
    print("Fitting M2 ETAN and barotropic velocity...")
    eta_hat, u_hat, v_hat, times = fit_tide(pairs, run_dir, diag_dir, ysl, xsl, args.period, args.spinup_periods)

    points = diagnostic_points(np.abs(eta_hat), np.hypot(np.abs(u_hat), np.abs(v_hat)), depth, slope, wet, xc, yc, args)
    outputs = plot_diagnostics(run_dir, save_dir, xc, yc, wet, depth, eta_hat, u_hat, v_hat, slope, points, args)
    summary = save_dir / "m2_barotropic_tide_diagnostic_summary.md"
    write_summary(summary, xc, yc, depth, slope, eta_hat, u_hat, v_hat, points, times, args)

    for name, (j, i) in points.items():
        print(f"{name}: lon={xc[j, i]:.6f}, lat={yc[j, i]:.6f}, depth={depth[j, i]:.1f} m")
    for out in outputs:
        print(f"Wrote {out}")
    print(f"Wrote {summary}")


if __name__ == "__main__":
    main()
