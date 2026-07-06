#!/usr/bin/env python3
"""Animate HR production depth slices with velocity quivers.

Creates one 3-panel top-view animation per requested depth:

1. vertical velocity, using cmocean.cm.balance;
2. temperature, using cmocean.cm.thermal;
3. horizontal speed, using cmocean.cm.speed;

Horizontal velocity quivers are overlaid on all panels.  By default the 60-cell
west/south/north sponge buffer is excluded, leaving the original 600 x 600
physical domain.

The script is designed for large MITgcm MDS diagnostics.  It memory-maps each
``diag_state_3d.*.data`` file and reads only the requested vertical level.
"""

from __future__ import annotations

import argparse
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter

try:
    import cmocean

    CMAPS = {
        "w": cmocean.cm.balance,
        "temp": cmocean.cm.thermal,
        "speed": cmocean.cm.speed,
    }
except Exception:  # pragma: no cover - fallback for lightweight environments
    CMAPS = {"w": "RdBu_r", "temp": "inferno", "speed": "viridis"}


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

    precision_match = re.search(r"dataprec\s*=\s*\[\s*'([^']+)'", text)
    if not precision_match:
        raise ValueError(f"Cannot parse dataprec from {path}")
    precision = precision_match.group(1).strip()
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
    data_path = run_dir / f"{name}.data"
    raw = np.fromfile(data_path, dtype=meta.dtype)
    expected = meta.nz * meta.ny * meta.nx
    if raw.size != expected:
        raise ValueError(f"{data_path}: expected {expected} values, got {raw.size}")
    array = raw.reshape(meta.nz, meta.ny, meta.nx)
    return array if meta.nz > 1 else array[0]


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
            # Fall back to iteration number if a hand-edited meta lacks time.
            match = re.search(r"\.(\d+)\.meta$", meta_path.name)
            time = float(match.group(1)) if match else np.nan
        else:
            time = meta.time
        frames.append(Frame(meta_path=meta_path, data_path=data_path, time_s=time))
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


def nearest_level(rc: np.ndarray, depth_m: float) -> tuple[int, float]:
    z = np.asarray(rc).reshape(-1)
    depths = np.abs(z)
    index = int(np.argmin(np.abs(depths - depth_m)))
    return index, float(depths[index])


def crop_slices(nx: int, ny: int, sponge: int) -> tuple[slice, slice]:
    """Return y,x slices that remove W/S/N sponge cells, not east."""
    if sponge < 0:
        raise ValueError("sponge must be non-negative")
    if 2 * sponge >= ny or sponge >= nx:
        raise ValueError(f"sponge={sponge} is incompatible with nx={nx}, ny={ny}")
    y = slice(sponge, ny - sponge)
    x = slice(sponge, nx)
    return y, x


def field_index(meta: Meta, field: str) -> int:
    if field not in meta.fields:
        raise ValueError(f"{field} absent from {meta.fields}")
    return meta.fields.index(field)


def robust_limits(
    frames: list[Frame],
    meta: Meta,
    field_indices: dict[str, int],
    k: int,
    ysl: slice,
    xsl: slice,
    wet_mask: np.ndarray,
    max_samples: int,
) -> dict[str, tuple[float, float]]:
    """Estimate useful color limits from a subset of frames."""
    if max_samples <= 0:
        chosen = frames
    else:
        n = min(max_samples, len(frames))
        chosen = [frames[i] for i in np.linspace(0, len(frames) - 1, n, dtype=int)]

    w_values, t_values, s_values = [], [], []
    for frame in chosen:
        arr = memmap_frame(frame, meta)
        u = np.asarray(arr[field_indices["UVEL"], k, ysl, xsl], dtype=float)
        v = np.asarray(arr[field_indices["VVEL"], k, ysl, xsl], dtype=float)
        w = np.asarray(arr[field_indices["WVEL"], k, ysl, xsl], dtype=float)
        t = np.asarray(arr[field_indices["THETA"], k, ysl, xsl], dtype=float)
        speed = np.hypot(u, v)
        valid = wet_mask & np.isfinite(w) & np.isfinite(t) & np.isfinite(speed)
        w_values.append(w[valid])
        t_values.append(t[valid])
        s_values.append(speed[valid])

    w_all = np.concatenate(w_values)
    t_all = np.concatenate(t_values)
    s_all = np.concatenate(s_values)

    w_abs = float(np.nanpercentile(np.abs(w_all), 98.0))
    if not np.isfinite(w_abs) or w_abs == 0:
        w_abs = float(np.nanmax(np.abs(w_all))) or 1.0

    t_lo, t_hi = np.nanpercentile(t_all, [2.0, 98.0])
    s_lo, s_hi = 0.0, np.nanpercentile(s_all, 98.0)
    if not np.isfinite(t_lo) or not np.isfinite(t_hi) or t_lo == t_hi:
        t_lo, t_hi = float(np.nanmin(t_all)), float(np.nanmax(t_all))
    if not np.isfinite(s_hi) or s_hi == 0:
        s_hi = float(np.nanmax(s_all)) or 1.0

    return {
        "w": (-w_abs, w_abs),
        "temp": (float(t_lo), float(t_hi)),
        "speed": (float(s_lo), float(s_hi)),
    }


def writer_for(path: Path, fps: int):
    suffix = path.suffix.lower()
    if suffix in {".gif"}:
        return PillowWriter(fps=fps)
    if not matplotlib.animation.writers.is_available("ffmpeg"):
        raise RuntimeError(
            "ffmpeg is not available to Matplotlib. Install ffmpeg or use --format gif."
        )
    return FFMpegWriter(fps=fps, codec="libx264", bitrate=2400)


def make_animation(
    run_dir: Path,
    output_dir: Path,
    frames: list[Frame],
    meta: Meta,
    depth_request: float,
    k: int,
    depth_actual: float,
    ysl: slice,
    xsl: slice,
    wet_mask: np.ndarray,
    q_stride: int,
    fps: int,
    fmt: str,
    period_s: float,
    spinup_s: float,
    limits: dict[str, tuple[float, float]],
    max_frames: int | None,
    dpi: int,
    mark_sponge_boundary: bool,
    sponge_boundary_cells: int,
) -> Path:
    fields = {name: field_index(meta, name) for name in ("UVEL", "VVEL", "WVEL", "THETA")}

    xc = read_grid_field(run_dir, "XC")[ysl, xsl]
    yc = read_grid_field(run_dir, "YC")[ysl, xsl]
    # Coordinates are lon/lat on the spherical-polar grid.
    x = np.asarray(xc)
    y = np.asarray(yc)

    selected = frames
    if max_frames is not None and max_frames > 0:
        selected = frames[:max_frames]

    qys = slice(None, None, q_stride)
    qxs = slice(None, None, q_stride)
    qx = x[qys, qxs]
    qy = y[qys, qxs]

    fig, axes = plt.subplots(1, 3, figsize=(17.5, 6.1), constrained_layout=True)
    titles = [
        "Vertical velocity WVEL",
        "Temperature THETA",
        "Horizontal speed",
    ]
    cmaps = [CMAPS["w"], CMAPS["temp"], CMAPS["speed"]]
    keys = ["w", "temp", "speed"]
    labels = ["m s$^{-1}$", "$^\\circ$C", "m s$^{-1}$"]

    first_arr = memmap_frame(selected[0], meta)
    u0 = np.asarray(first_arr[fields["UVEL"], k, ysl, xsl], dtype=float)
    v0 = np.asarray(first_arr[fields["VVEL"], k, ysl, xsl], dtype=float)
    panels = [
        np.ma.masked_where(~wet_mask, np.asarray(first_arr[fields["WVEL"], k, ysl, xsl], dtype=float)),
        np.ma.masked_where(~wet_mask, np.asarray(first_arr[fields["THETA"], k, ysl, xsl], dtype=float)),
        np.ma.masked_where(~wet_mask, np.hypot(u0, v0)),
    ]
    u0 = np.where(wet_mask, u0, np.nan)
    v0 = np.where(wet_mask, v0, np.nan)

    meshes = []
    quivers = []
    for ax, title, panel, cmap, key, label in zip(axes, titles, panels, cmaps, keys, labels):
        vmin, vmax = limits[key]
        mesh = ax.pcolormesh(x, y, panel, shading="auto", cmap=cmap, vmin=vmin, vmax=vmax)
        q = ax.quiver(
            qx,
            qy,
            u0[qys, qxs],
            v0[qys, qxs],
            color="k",
            alpha=0.58,
            pivot="middle",
            scale=2.0,
            width=0.0021,
            headwidth=3.0,
        )
        ax.set_title(title)
        ax.set_xlabel("Longitude")
        ax.set_aspect("equal", adjustable="box")
        if mark_sponge_boundary and sponge_boundary_cells > 0:
            # Mark the edge of the original 600 x 600 physical domain.  The
            # telescopic buffer occupies W/S/N: i < 60, j < 60, j >= 660.
            west_lon = float(read_grid_field(run_dir, "XC")[0, sponge_boundary_cells])
            south_lat = float(read_grid_field(run_dir, "YC")[sponge_boundary_cells, 0])
            north_lat = float(read_grid_field(run_dir, "YC")[-sponge_boundary_cells, 0])
            ax.axvline(west_lon, color="k", lw=1.2, ls="--", alpha=0.75)
            ax.axhline(south_lat, color="k", lw=1.2, ls="--", alpha=0.75)
            ax.axhline(north_lat, color="k", lw=1.2, ls="--", alpha=0.75)
        fig.colorbar(mesh, ax=ax, shrink=0.86, label=label)
        meshes.append(mesh)
        quivers.append(q)
    axes[0].set_ylabel("Latitude")

    timeline = fig.add_axes([0.18, 0.018, 0.64, 0.035])
    run_start = 0.0
    run_end = max(frame.time_s for frame in selected)
    timeline.axvspan(run_start / period_s, spinup_s / period_s, color="0.85", label="spin-up")
    timeline.axvspan(spinup_s / period_s, run_end / period_s, color="#d8ecff", label="analysis")
    current_line = timeline.axvline(selected[0].time_s / period_s, color="k", lw=2)
    timeline.set_xlim(run_start / period_s, run_end / period_s)
    timeline.set_yticks([])
    timeline.set_xlabel("time since model start (M2 periods)")
    timeline.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False)

    title_text = fig.suptitle("")

    def update(i: int):
        frame = selected[i]
        arr = memmap_frame(frame, meta)
        u = np.asarray(arr[fields["UVEL"], k, ysl, xsl], dtype=float)
        v = np.asarray(arr[fields["VVEL"], k, ysl, xsl], dtype=float)
        w = np.asarray(arr[fields["WVEL"], k, ysl, xsl], dtype=float)
        temp = np.asarray(arr[fields["THETA"], k, ysl, xsl], dtype=float)
        speed = np.hypot(u, v)
        w = np.ma.masked_where(~wet_mask, w)
        temp = np.ma.masked_where(~wet_mask, temp)
        speed = np.ma.masked_where(~wet_mask, speed)
        u = np.where(wet_mask, u, np.nan)
        v = np.where(wet_mask, v, np.nan)

        for mesh, panel in zip(meshes, (w, temp, speed)):
            mesh.set_array(panel.ravel())
        for q in quivers:
            q.set_UVC(u[qys, qxs], v[qys, qxs])

        current_line.set_xdata([frame.time_s / period_s, frame.time_s / period_s])
        analysis_time = max(0.0, frame.time_s - spinup_s)
        title_text.set_text(
            f"Depth request {depth_request:g} m; model level {k} at {depth_actual:.1f} m | "
            f"t = {frame.time_s / 86400:.2f} d = {frame.time_s / period_s:.2f} M2; "
            f"analysis +{analysis_time / period_s:.2f} M2"
        )
        return [*meshes, *quivers, current_line, title_text]

    suffix = ".gif" if fmt == "gif" else ".mp4"
    out_path = output_dir / f"depth_{int(round(depth_request)):03d}m_w_temp_speed_quiver{suffix}"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    anim = FuncAnimation(fig, update, frames=len(selected), blit=False)
    anim.save(out_path, writer=writer_for(out_path, fps=fps), dpi=dpi)
    plt.close(fig)
    return out_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument(
        "--diag-dir",
        type=Path,
        default=None,
        help="Directory containing diag_state_3d files; default: RUN_DIR/output.",
    )
    parser.add_argument(
        "--save-dir",
        type=Path,
        default=None,
        help="Directory for animations; default: DIAG_DIR/animations.",
    )
    parser.add_argument("--prefix", default="diag_state_3d")
    parser.add_argument("--depths", nargs="+", type=float, default=[5, 50, 100, 150, 300])
    parser.add_argument("--sponge-cells", type=int, default=60)
    parser.add_argument(
        "--sponge-boundary-cells",
        type=int,
        default=60,
        help="Cell count used only for dashed physical-domain markers.",
    )
    parser.add_argument(
        "--mark-sponge-boundary",
        action="store_true",
        help="Draw dashed lines at the original physical-domain edge.",
    )
    parser.add_argument("--quiver-stride", type=int, default=24)
    parser.add_argument("--fps", type=int, default=5)
    parser.add_argument("--format", choices=["mp4", "gif"], default="mp4")
    parser.add_argument("--period", type=float, default=44712.0, help="M2 period in seconds.")
    parser.add_argument("--spinup-periods", type=float, default=3.0)
    parser.add_argument(
        "--limit-samples",
        type=int,
        default=12,
        help="Number of frames used to estimate robust color limits; 0 means all frames.",
    )
    parser.add_argument(
        "--max-frames",
        type=int,
        default=None,
        help="Optional quick-test limit on number of animation frames.",
    )
    parser.add_argument("--dpi", type=int, default=145)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_dir = args.run_dir.resolve()
    diag_dir = (args.diag_dir or (run_dir / "output")).resolve()
    save_dir = (args.save_dir or (diag_dir / "animations")).resolve()
    frames = diagnostic_frames(diag_dir, args.prefix)
    first_meta = parse_meta(frames[0].meta_path)
    required = {"UVEL", "VVEL", "WVEL", "THETA"}
    missing = sorted(required.difference(first_meta.fields))
    if missing:
        raise ValueError(f"{frames[0].meta_path} is missing required fields: {missing}")

    rc = read_grid_field(run_dir, "RC")
    hfac = read_grid_field(run_dir, "hFacC")
    ysl, xsl = crop_slices(first_meta.nx, first_meta.ny, args.sponge_cells)
    field_indices = {name: field_index(first_meta, name) for name in required}
    spinup_s = args.spinup_periods * args.period

    print(f"Run directory: {run_dir}")
    print(f"Diagnostics:   {diag_dir}")
    print(f"Animations:    {save_dir}")
    print(f"Frames found:  {len(frames)}")
    print(
        "Crop:         "
        f"j={ysl.start}:{ysl.stop}, i={xsl.start}:{xsl.stop} "
        f"({ysl.stop - ysl.start} x {xsl.stop - xsl.start})"
    )
    print(f"Spin-up:      {spinup_s:.0f} s = {args.spinup_periods:g} M2")

    for depth in args.depths:
        k, actual = nearest_level(rc, depth)
        wet_mask = np.asarray(hfac[k, ysl, xsl] > 0)
        print(f"\nDepth {depth:g} m -> model level k={k}, |RC|={actual:.2f} m")
        limits = robust_limits(
            frames,
            first_meta,
            field_indices,
            k,
            ysl,
            xsl,
            wet_mask=wet_mask,
            max_samples=args.limit_samples,
        )
        print(
            "Color limits: "
            f"W={limits['w']}, THETA={limits['temp']}, speed={limits['speed']}"
        )
        out = make_animation(
            run_dir=run_dir,
            output_dir=save_dir,
            frames=frames,
            meta=first_meta,
            depth_request=depth,
            k=k,
            depth_actual=actual,
            ysl=ysl,
            xsl=xsl,
            wet_mask=wet_mask,
            q_stride=args.quiver_stride,
            fps=args.fps,
            fmt=args.format,
            period_s=args.period,
            spinup_s=spinup_s,
            limits=limits,
            max_frames=args.max_frames,
            dpi=args.dpi,
            mark_sponge_boundary=args.mark_sponge_boundary,
            sponge_boundary_cells=args.sponge_boundary_cells,
        )
        print(f"Wrote {out}")


if __name__ == "__main__":
    main()
