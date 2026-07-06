#!/usr/bin/env python3
"""Animate free surface and M2 phase from production 2-D diagnostics.

Creates a two-panel top-view animation:

1. ETAN/free surface, using ``cmocean.cm.delta``;
2. M2 phase, using ``cmocean.cm.phase``.

The phase panel is based on a post-spinup harmonic fit to ETAN.  By default it
shows the evolving phase ``angle(eta_hat * exp(i omega t))`` at each animation
time.  Use ``--phase-mode static`` to show the fitted spatial phase lag only.

By default the 60-cell west/south/north sponge buffer is excluded, matching the
physical-domain movies.  Use ``--sponge-cells 0 --mark-sponge-boundary`` to see
the full telescopic domain with dashed physical-domain markers.
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
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter

try:
    import cmocean

    CMAP_DELTA = cmocean.cm.delta
    CMAP_PHASE = cmocean.cm.phase
except Exception:  # pragma: no cover
    CMAP_DELTA = "RdBu_r"
    CMAP_PHASE = "twilight"


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
class Frame:
    meta_path: Path
    data_path: Path
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


def read_grid_field(run_dir: Path, name: str) -> np.ndarray:
    meta = parse_meta(run_dir / f"{name}.meta")
    raw = np.fromfile(run_dir / f"{name}.data", dtype=meta.dtype)
    expected = meta.nz * meta.ny * meta.nx
    if raw.size != expected:
        raise ValueError(f"{name}: expected {expected} values, got {raw.size}")
    result = raw.reshape(meta.nz, meta.ny, meta.nx)
    return result if meta.nz > 1 else result[0]


def diagnostic_frames(output_dir: Path, prefix: str) -> list[Frame]:
    frames: list[Frame] = []
    for meta_path in sorted(output_dir.glob(f"{prefix}.[0-9]*.meta")):
        if meta_path.name.startswith("._"):
            continue
        data_path = meta_path.with_suffix(".data")
        if not data_path.exists():
            continue
        meta = parse_meta(meta_path)
        if meta.time is None:
            match = re.search(r"\.(\d+)\.meta$", meta_path.name)
            time = float(match.group(1)) if match else np.nan
        else:
            time = meta.time
        frames.append(Frame(meta_path, data_path, time))
    if not frames:
        raise FileNotFoundError(f"No {prefix}.[0-9]*.meta files found in {output_dir}")
    return frames


def memmap_frame(frame: Frame, meta: Meta) -> np.memmap:
    expected = meta.nrecords * meta.nz * meta.ny * meta.nx
    actual = frame.data_path.stat().st_size // meta.dtype.itemsize
    if actual != expected:
        raise ValueError(f"{frame.data_path}: expected {expected} values, got {actual}")
    return np.memmap(
        frame.data_path,
        dtype=meta.dtype,
        mode="r",
        shape=(meta.nrecords, meta.nz, meta.ny, meta.nx),
    )


def field_index(meta: Meta, field: str) -> int:
    if field not in meta.fields:
        raise ValueError(f"{field} absent from {meta.fields}")
    return meta.fields.index(field)


def crop_slices(nx: int, ny: int, sponge: int) -> tuple[slice, slice]:
    if sponge < 0:
        raise ValueError("sponge must be non-negative")
    if 2 * sponge >= ny or sponge >= nx:
        raise ValueError(f"sponge={sponge} incompatible with nx={nx}, ny={ny}")
    return slice(sponge, ny - sponge), slice(sponge, nx)


def load_etan_stack(
    frames: list[Frame],
    meta: Meta,
    etan_index: int,
    ysl: slice,
    xsl: slice,
) -> tuple[np.ndarray, np.ndarray]:
    times = np.asarray([f.time_s for f in frames], dtype=float)
    sample_shape = (ysl.stop - ysl.start, xsl.stop - xsl.start)
    stack = np.empty((len(frames), *sample_shape), dtype=np.float64)
    for n, frame in enumerate(frames):
        arr = memmap_frame(frame, meta)
        stack[n] = np.asarray(arr[etan_index, 0, ysl, xsl], dtype=float)
    return times, stack


def fit_harmonic(values: np.ndarray, times: np.ndarray, period: float) -> np.ndarray:
    """Fit constant, trend, cos, sin; return C for Re(C exp(i omega t))."""
    centered = (times - times.mean()) / period
    omega_t = 2 * np.pi * times / period
    design = np.column_stack(
        (np.ones(times.size), centered, np.cos(omega_t), np.sin(omega_t))
    )
    coefficients = np.linalg.pinv(design) @ values.reshape(times.size, -1)
    harmonic = coefficients[2] - 1j * coefficients[3]
    return harmonic.reshape(values.shape[1:])


def writer_for(path: Path, fps: int):
    if path.suffix.lower() == ".gif":
        return PillowWriter(fps=fps)
    if not matplotlib.animation.writers.is_available("ffmpeg"):
        raise RuntimeError(
            "ffmpeg is not available to Matplotlib. Install ffmpeg or use --format gif."
        )
    return FFMpegWriter(fps=fps, codec="libx264", bitrate=2400)


def robust_etan_limit(stack: np.ndarray, wet_mask: np.ndarray) -> float:
    valid = wet_mask[None, :, :] & np.isfinite(stack)
    values = stack[valid]
    limit = float(np.nanpercentile(np.abs(values), 98.0))
    if not np.isfinite(limit) or limit == 0:
        limit = float(np.nanmax(np.abs(values))) or 1.0
    return limit


def make_animation(args: argparse.Namespace) -> Path:
    run_dir = args.run_dir.resolve()
    diag_dir = (args.diag_dir or (run_dir / "output")).resolve()
    save_dir = (args.save_dir or (diag_dir / "animations_eta_phase")).resolve()
    frames = diagnostic_frames(diag_dir, args.prefix)
    meta = parse_meta(frames[0].meta_path)
    etan_index = field_index(meta, "ETAN")
    ysl, xsl = crop_slices(meta.nx, meta.ny, args.sponge_cells)

    xc = read_grid_field(run_dir, "XC")[ysl, xsl]
    yc = read_grid_field(run_dir, "YC")[ysl, xsl]
    hfac0 = read_grid_field(run_dir, "hFacC")[0, ysl, xsl]
    wet_mask = hfac0 > 0

    times, etan = load_etan_stack(frames, meta, etan_index, ysl, xsl)
    keep = times >= args.spinup_periods * args.period - 0.1
    if np.count_nonzero(keep) < 4:
        raise ValueError("Fewer than four post-spinup frames are available for phase fit")
    eta_hat = fit_harmonic(etan[keep], times[keep], args.period)
    eta_amp = np.abs(eta_hat)
    phase_static = np.angle(eta_hat)
    phase_static = np.ma.masked_where((~wet_mask) | (eta_amp <= args.min_amplitude), phase_static)

    eta_lim = robust_etan_limit(etan, wet_mask)
    etan_mask = np.broadcast_to(~wet_mask[None, :, :], etan.shape)
    etan = np.ma.masked_where(etan_mask, etan)

    selected_count = len(frames) if args.max_frames is None else min(args.max_frames, len(frames))
    selected = np.arange(selected_count)

    fig, axes = plt.subplots(1, 2, figsize=(12.4, 6.2), constrained_layout=True)
    eta_mesh = axes[0].pcolormesh(
        xc, yc, etan[0], shading="auto", cmap=CMAP_DELTA, vmin=-eta_lim, vmax=eta_lim
    )
    phase0 = phase_static
    if args.phase_mode == "instantaneous":
        phase0 = np.ma.masked_where(
            phase_static.mask,
            np.angle(eta_hat * np.exp(1j * 2 * np.pi * times[0] / args.period)),
        )
    phase_mesh = axes[1].pcolormesh(
        xc, yc, phase0, shading="auto", cmap=CMAP_PHASE, vmin=-np.pi, vmax=np.pi
    )

    axes[0].set_title("Free surface ETAN")
    axes[1].set_title("M2 phase")
    for ax in axes:
        ax.set_xlabel("Longitude")
        ax.set_aspect("equal", adjustable="box")
        if args.mark_sponge_boundary and args.sponge_boundary_cells > 0:
            full_xc = read_grid_field(run_dir, "XC")
            full_yc = read_grid_field(run_dir, "YC")
            west_lon = float(full_xc[0, args.sponge_boundary_cells])
            south_lat = float(full_yc[args.sponge_boundary_cells, 0])
            north_lat = float(full_yc[-args.sponge_boundary_cells, 0])
            ax.axvline(west_lon, color="k", lw=1.2, ls="--", alpha=0.75)
            ax.axhline(south_lat, color="k", lw=1.2, ls="--", alpha=0.75)
            ax.axhline(north_lat, color="k", lw=1.2, ls="--", alpha=0.75)
    axes[0].set_ylabel("Latitude")
    fig.colorbar(eta_mesh, ax=axes[0], shrink=0.88, label="m")
    cbar = fig.colorbar(phase_mesh, ax=axes[1], shrink=0.88, label="radians")
    cbar.set_ticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi])
    cbar.set_ticklabels(["$-\\pi$", "$-\\pi/2$", "0", "$\\pi/2$", "$\\pi$"])

    timeline = fig.add_axes([0.20, 0.018, 0.58, 0.035])
    spinup_s = args.spinup_periods * args.period
    run_end = max(times[:selected_count])
    timeline.axvspan(0.0, spinup_s / args.period, color="0.85", label="spin-up")
    timeline.axvspan(spinup_s / args.period, run_end / args.period, color="#d8ecff", label="analysis")
    current_line = timeline.axvline(times[0] / args.period, color="k", lw=2)
    timeline.set_xlim(0.0, run_end / args.period)
    timeline.set_yticks([])
    timeline.set_xlabel("time since model start (M2 periods)")
    timeline.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False)
    title = fig.suptitle("")

    def update(i: int):
        idx = selected[i]
        t = times[idx]
        eta_mesh.set_array(etan[idx].ravel())
        if args.phase_mode == "instantaneous":
            phase = np.angle(eta_hat * np.exp(1j * 2 * np.pi * t / args.period))
            phase = np.ma.masked_where(phase_static.mask, phase)
        else:
            phase = phase_static
        phase_mesh.set_array(phase.ravel())
        current_line.set_xdata([t / args.period, t / args.period])
        title.set_text(
            f"ETAN and M2 phase | t = {t / 86400:.2f} d = {t / args.period:.2f} M2; "
            f"analysis +{max(0.0, t - spinup_s) / args.period:.2f} M2"
        )
        return eta_mesh, phase_mesh, current_line, title

    suffix = ".gif" if args.format == "gif" else ".mp4"
    out_path = save_dir / f"etan_m2_phase_{args.phase_mode}{suffix}"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    anim = FuncAnimation(fig, update, frames=len(selected), blit=False)
    anim.save(out_path, writer=writer_for(out_path, args.fps), dpi=args.dpi)
    plt.close(fig)
    return out_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--diag-dir", type=Path, default=None)
    parser.add_argument("--save-dir", type=Path, default=None)
    parser.add_argument("--prefix", default="diag_state_2d")
    parser.add_argument("--sponge-cells", type=int, default=60)
    parser.add_argument("--sponge-boundary-cells", type=int, default=60)
    parser.add_argument("--mark-sponge-boundary", action="store_true")
    parser.add_argument("--period", type=float, default=44712.0)
    parser.add_argument("--spinup-periods", type=float, default=3.0)
    parser.add_argument(
        "--phase-mode",
        choices=["instantaneous", "static"],
        default="instantaneous",
        help="instantaneous: angle(eta_hat exp(i omega t)); static: angle(eta_hat).",
    )
    parser.add_argument(
        "--min-amplitude",
        type=float,
        default=1e-5,
        help="Mask phase where fitted ETAN M2 amplitude is below this value in m.",
    )
    parser.add_argument("--fps", type=int, default=5)
    parser.add_argument("--format", choices=["mp4", "gif"], default="mp4")
    parser.add_argument("--max-frames", type=int, default=None)
    parser.add_argument("--dpi", type=int, default=145)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    out = make_animation(args)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
