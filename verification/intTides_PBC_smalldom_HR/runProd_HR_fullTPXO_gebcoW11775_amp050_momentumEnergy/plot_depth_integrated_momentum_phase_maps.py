#!/usr/bin/env python3
"""Plot depth-integrated momentum-budget terms over an average M2 period.

The production run writes time-averaged momentum diagnostics every M2/8.  This
script groups those averages by tidal phase, averages across complete M2 cycles,
depth-integrates each U/V tendency pair, interpolates to tracer-cell centers,
and makes 8 maps per term.

By default the first averaged diagnostic interval is skipped.  In this run it
ends exactly at the 3-M2 ramp cutoff and includes the previous M2/8 interval,
so skipping it leaves 40 intervals = 5 full M2 cycles x 8 phase bins.

Plots show depth-integrated vector magnitude as color and a subsampled quiver
for direction.  Units are m^2 s^-2 after vertical integration.  TOTUTEND and
TOTVTEND are converted from MITgcm's m/s/day to m/s^2 before integration.

Default plots mask the 60-cell west/south/north buffer zones and show the
original physical domain.  Use --include-buffers to plot the full computational
domain.
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np


SECONDS_PER_DAY = 86400.0


@dataclass
class Meta:
    nx: int
    ny: int
    nz: int
    nrecords: int
    dtype: np.dtype
    fields: list[str]
    time_start: float | None
    time_end: float | None


@dataclass
class Term:
    name: str
    pretty: str
    prefix: str
    u_field: str
    v_field: str
    scale: float = 1.0


TERMS = [
    Term("total_tendency", "Total tendency", "diag_mom_budget", "TOTUTEND", "TOTVTEND", 1.0 / SECONDS_PER_DAY),
    Term("advection_total", "Advection, total", "diag_mom_budget", "Um_Advec", "Vm_Advec"),
    Term("pressure_gradient", "Pressure gradient", "diag_mom_budget", "Um_dPhiX", "Vm_dPhiY"),
    Term("external_forcing", "External forcing", "diag_mom_budget", "Um_Ext", "Vm_Ext"),
    Term("dissipation_total", "Dissipation, total", "diag_mom_budget", "Um_Diss", "Vm_Diss"),
    Term("adams_bashforth", "Adams-Bashforth tendency", "diag_mom_budget", "AB_gU", "AB_gV"),
    Term("coriolis", "Coriolis", "diag_mom_vecinv", "Um_Cori", "Vm_Cori"),
    Term("vertical_advection", "Vertical advection", "diag_mom_vecinv", "Um_AdvZ3", "Vm_AdvZ3"),
    Term("relative_vorticity_advection", "Relative-vorticity advection", "diag_mom_vecinv", "Um_AdvRe", "Vm_AdvRe"),
    Term("horizontal_dissipation", "Horizontal dissipation", "diag_mom_vecinv", "Um_hDis2", "Vm_hDis2"),
    Term("side_drag", "Side drag", "diag_mom_vecinv", "USidDrag", "VSidDrag"),
    Term("bottom_drag", "Bottom drag", "diag_mom_vecinv", "UBotDrag", "VBotDrag"),
]


def parse_meta(path: Path) -> Meta:
    text = path.read_text()
    dim_match = re.search(r"dimList\s*=\s*\[(.*?)\];", text, re.S)
    if not dim_match:
        raise ValueError(f"Cannot parse dimList from {path}")
    numbers = [int(x) for x in re.findall(r"[-+]?\d+", dim_match.group(1))]
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
    nrec_match = re.search(r"nrecords\s*=\s*\[\s*(\d+)", text)
    nrecords = int(nrec_match.group(1)) if nrec_match else 1
    field_match = re.search(r"fldList\s*=\s*\{(.*?)\};", text, re.S)
    fields = re.findall(r"'([^']+)'", field_match.group(1)) if field_match else []
    fields = [field.strip() for field in fields]
    time_match = re.search(r"timeInterval\s*=\s*\[\s*([^\]]+)", text)
    times = [float(x) for x in re.findall(r"[-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?", time_match.group(1))] if time_match else []
    if len(times) == 1:
        time_start = None
        time_end = times[0]
    elif len(times) >= 2:
        time_start, time_end = times[0], times[1]
    else:
        time_start = None
        time_end = None
    return Meta(nx, ny, nz, nrecords, dtype, fields, time_start, time_end)


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


def depth_integrate_to_c(
    u_term: np.ndarray,
    v_term: np.ndarray,
    drf: np.ndarray,
    hfac_w: np.ndarray,
    hfac_s: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Depth-integrate U/V tendencies and interpolate to tracer-cell centers."""
    u_face = np.sum(u_term * hfac_w * drf[:, None, None], axis=0)
    v_face = np.sum(v_term * hfac_s * drf[:, None, None], axis=0)

    u_c = np.empty_like(u_face, dtype=np.float64)
    v_c = np.empty_like(v_face, dtype=np.float64)
    u_c[:, :-1] = 0.5 * (u_face[:, :-1] + u_face[:, 1:])
    u_c[:, -1] = u_face[:, -1]
    v_c[:-1, :] = 0.5 * (v_face[:-1, :] + v_face[1:, :])
    v_c[-1, :] = v_face[-1, :]
    return u_c, v_c


def physical_mask(wet: np.ndarray, buffer_cells: int) -> np.ndarray:
    jj, ii = np.indices(wet.shape)
    return wet & (ii >= buffer_cells) & (jj >= buffer_cells) & (jj < wet.shape[0] - buffer_cells)


def collect_phase_means(
    run_dir: Path,
    diag_dir: Path,
    terms: list[Term],
    nphase: int,
    skip_first: bool,
) -> tuple[dict[str, dict[str, np.ndarray]], list[float]]:
    drf = read_grid(run_dir, "DRF").reshape(-1)
    hfac_w = read_grid(run_dir, "hFacW")
    hfac_s = read_grid(run_dir, "hFacS")
    ny, nx = hfac_w.shape[1:]

    sums: dict[str, dict[str, np.ndarray]] = {
        term.name: {
            "u": np.zeros((nphase, ny, nx), dtype=np.float64),
            "v": np.zeros((nphase, ny, nx), dtype=np.float64),
            "count": np.zeros(nphase, dtype=np.int64),
        }
        for term in terms
    }
    phase_mid_times = np.zeros(nphase, dtype=np.float64)

    # Group terms by source diagnostic stream to avoid opening files repeatedly.
    by_prefix: dict[str, list[Term]] = {}
    for term in terms:
        by_prefix.setdefault(term.prefix, []).append(term)

    for prefix, prefix_terms in by_prefix.items():
        files = sorted(diag_dir.glob(f"{prefix}.[0-9]*.meta"))
        if skip_first and files:
            files = files[1:]
        usable = (len(files) // nphase) * nphase
        if usable != len(files):
            print(f"{prefix}: using first {usable} of {len(files)} files to keep complete phase cycles")
            files = files[:usable]
        if not files:
            raise FileNotFoundError(f"No complete {prefix} diagnostics found in {diag_dir}")

        for sample_index, meta_path in enumerate(files):
            phase = sample_index % nphase
            array, meta = memmap_mds(meta_path)
            missing = [
                field
                for term in prefix_terms
                for field in (term.u_field, term.v_field)
                if field not in meta.fields
            ]
            if missing:
                raise ValueError(f"{meta_path}: missing fields {sorted(set(missing))}")
            mid_time = (
                0.5 * (meta.time_start + meta.time_end)
                if meta.time_start is not None and meta.time_end is not None
                else float(sample_index)
            )
            if sums[prefix_terms[0].name]["count"][phase] == 0:
                phase_mid_times[phase] = mid_time

            for term in prefix_terms:
                u = np.asarray(array[meta.fields.index(term.u_field)], dtype=np.float64) * term.scale
                v = np.asarray(array[meta.fields.index(term.v_field)], dtype=np.float64) * term.scale
                u_c, v_c = depth_integrate_to_c(u, v, drf, hfac_w, hfac_s)
                sums[term.name]["u"][phase] += u_c
                sums[term.name]["v"][phase] += v_c
                sums[term.name]["count"][phase] += 1
            print(f"{prefix}: processed {sample_index + 1:03d}/{len(files):03d} {meta_path.name}", flush=True)

    # Each term count was incremented once per term per phase sample.  Because
    # all terms in a prefix share files, counts are valid per term.
    for term in terms:
        count = sums[term.name]["count"]
        if np.any(count == 0):
            raise ValueError(f"{term.name}: missing phase counts {count}")
        sums[term.name]["u"] /= count[:, None, None]
        sums[term.name]["v"] /= count[:, None, None]
    return sums, list(phase_mid_times)


def write_term_plots(
    out_dir: Path,
    xc: np.ndarray,
    yc: np.ndarray,
    mask: np.ndarray,
    term: Term,
    u_phase: np.ndarray,
    v_phase: np.ndarray,
    quiver_stride: int,
) -> None:
    import os

    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    term_dir = out_dir / term.name
    term_dir.mkdir(parents=True, exist_ok=True)

    mag = np.sqrt(u_phase**2 + v_phase**2)
    mag[:, ~mask] = np.nan
    scale_max = np.nanpercentile(mag, 99.0)
    if not np.isfinite(scale_max) or scale_max <= 0:
        scale_max = 1.0

    xs = xc[::quiver_stride, ::quiver_stride]
    ys = yc[::quiver_stride, ::quiver_stride]
    mm = mask[::quiver_stride, ::quiver_stride]

    for phase in range(u_phase.shape[0]):
        u = u_phase[phase].copy()
        v = v_phase[phase].copy()
        m = mag[phase].copy()
        u[~mask] = np.nan
        v[~mask] = np.nan
        q_u = u[::quiver_stride, ::quiver_stride]
        q_v = v[::quiver_stride, ::quiver_stride]
        q_u[~mm] = np.nan
        q_v[~mm] = np.nan

        fig, ax = plt.subplots(figsize=(8, 7), constrained_layout=True)
        pcm = ax.pcolormesh(xc, yc, m, shading="auto", cmap="magma", vmin=0, vmax=scale_max)
        ax.quiver(xs, ys, q_u, q_v, color="k", pivot="mid", scale=None, width=0.0015, alpha=0.65)
        cb = fig.colorbar(pcm, ax=ax, shrink=0.85)
        cb.set_label(r"Depth-integrated tendency magnitude (m$^2$ s$^{-2}$)")
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_title(f"{term.pretty}: phase {phase + 1}/8")
        ax.set_aspect("equal", adjustable="box")
        fig.savefig(term_dir / f"{term.name}_phase{phase + 1:02d}.png", dpi=180)
        plt.close(fig)


def write_summary_csv(path: Path, terms: list[Term], phase_means: dict[str, dict[str, np.ndarray]], mask: np.ndarray) -> None:
    import csv

    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["term", "phase", "mean_magnitude_m2_s2", "max_magnitude_m2_s2"])
        for term in terms:
            u = phase_means[term.name]["u"]
            v = phase_means[term.name]["v"]
            mag = np.sqrt(u**2 + v**2)
            for phase in range(mag.shape[0]):
                values = mag[phase][mask]
                writer.writerow([
                    term.name,
                    phase + 1,
                    f"{np.nanmean(values):.16e}",
                    f"{np.nanmax(values):.16e}",
                ])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--output-dir", type=Path, default=None,
                        help="Default: RUN_DIR/output/postprocessing/momentum_phase_maps")
    parser.add_argument("--nphase", type=int, default=8)
    parser.add_argument("--buffer-cells", type=int, default=60)
    parser.add_argument("--include-buffers", action="store_true",
                        help="Plot the full computational domain instead of masking buffer zones.")
    parser.add_argument("--no-skip-first", action="store_true",
                        help="Do not skip the first averaged interval ending at the ramp cutoff.")
    parser.add_argument("--quiver-stride", type=int, default=22)
    parser.add_argument("--terms", default="all",
                        help="Comma-separated term names, or 'all'.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_dir = args.run_dir.resolve()
    diag_dir = run_dir / "output"
    out_dir = args.output_dir or (diag_dir / "postprocessing" / "momentum_phase_maps")
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.terms == "all":
        terms = TERMS
    else:
        wanted = {item.strip() for item in args.terms.split(",") if item.strip()}
        lookup = {term.name: term for term in TERMS}
        missing = sorted(wanted - set(lookup))
        if missing:
            raise ValueError(f"Unknown term(s): {missing}. Valid terms: {sorted(lookup)}")
        terms = [lookup[name] for name in wanted]

    xc = read_grid(run_dir, "XC")
    yc = read_grid(run_dir, "YC")
    depth = read_grid(run_dir, "Depth")
    hfac_c = read_grid(run_dir, "hFacC")
    wet = (hfac_c[0] > 0) & (depth > 0)
    plot_mask = wet if args.include_buffers else physical_mask(wet, args.buffer_cells)

    phase_means, phase_mid_times = collect_phase_means(
        run_dir=run_dir,
        diag_dir=diag_dir,
        terms=terms,
        nphase=args.nphase,
        skip_first=not args.no_skip_first,
    )

    for term in terms:
        print(f"plotting {term.name}", flush=True)
        write_term_plots(
            out_dir=out_dir,
            xc=xc,
            yc=yc,
            mask=plot_mask,
            term=term,
            u_phase=phase_means[term.name]["u"],
            v_phase=phase_means[term.name]["v"],
            quiver_stride=args.quiver_stride,
        )

    summary_path = out_dir / "momentum_phase_map_summary.csv"
    write_summary_csv(summary_path, terms, phase_means, plot_mask)

    info_path = out_dir / "README.txt"
    with info_path.open("w") as handle:
        handle.write("Depth-integrated momentum-budget phase maps\n")
        handle.write("===========================================\n\n")
        handle.write(f"Run directory: {run_dir}\n")
        handle.write(f"Phase bins: {args.nphase}\n")
        handle.write(f"Skipped first averaged interval: {not args.no_skip_first}\n")
        handle.write(f"Included buffer zones: {args.include_buffers}\n")
        handle.write(f"Buffer cells: {args.buffer_cells}\n")
        handle.write("Units: m^2 s^-2 after depth integration. TOTUTEND/TOTVTEND converted from m/s/day.\n\n")
        handle.write("Terms:\n")
        for term in terms:
            handle.write(f"  {term.name}: {term.pretty} ({term.u_field}, {term.v_field})\n")
        handle.write("\nRepresentative midpoint times for phases, seconds:\n")
        for phase, time in enumerate(phase_mid_times, start=1):
            handle.write(f"  phase {phase}: {time:.1f}\n")

    print("\nWrote momentum phase maps to:")
    print(f"  {out_dir}")
    print(f"Summary:")
    print(f"  {summary_path}")
    print(f"Notes:")
    print(f"  {info_path}")


if __name__ == "__main__":
    main()

