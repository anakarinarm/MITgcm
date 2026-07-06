#!/usr/bin/env python3
"""Plot M2 ETAN amplitude/phase summary panels for TPXO boundary tests.

This is a visual companion to ``compare_bc_eta_amphidromes_to_tpxo.py``.
It fits M2 SSH harmonics from each boundary-test ``bc_wave_eta`` stream,
interpolates TPXO9 M2 elevation onto the same grid, and makes summary figures
with one column per case and two rows:

1. M2 ETAN amplitude
2. M2 ETAN phase

Two amplitude versions are written:

* shared amplitude scale: best for absolute comparison against TPXO;
* local amplitude scale: best for seeing spatial structure in weak-response
  runs where a TPXO-scale colorbar would make everything look flat.
"""

from __future__ import annotations

import argparse
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
except Exception:  # pragma: no cover
    CMAP_AMP = "viridis"
    CMAP_PHASE = "twilight"

from compare_bc_eta_amphidromes_to_tpxo import (
    collect_eta_files,
    crop_slices,
    default_runs,
    fit_eta_harmonic,
    interpolate_tpxo_eta,
    masked,
    parse_meta,
    read_grid,
)


@dataclass
class Case:
    name: str
    eta: np.ndarray
    xc: np.ndarray
    yc: np.ndarray
    depth: np.ndarray
    wet: np.ndarray
    amp_p95: float
    amp_max: float


def short_name(name: str) -> str:
    replacements = {
        "runBC_LR_fullTPXO_": "",
        "obeta_balanced": "SSH+vel balanced",
        "obeta_weakvel025": "SSH + 25% vel",
        "obeta_weakvel010": "SSH + 10% vel",
        "obeta_zerovel": "SSH only, vel=0",
        "obeta": "SSH+vel raw",
        "balanced": "vel balanced, no SSH",
    }
    out = name
    for old, new in replacements.items():
        out = out.replace(old, new)
    return out


def load_case(run_dir: Path, args: argparse.Namespace) -> tuple[Case, np.ndarray]:
    eta_files = collect_eta_files(run_dir)
    meta0 = parse_meta(eta_files[0])
    ysl, xsl = crop_slices(meta0.nx, meta0.ny, 0 if args.include_sponge else args.sponge_cells)
    xc = read_grid(run_dir, "XC")[ysl, xsl]
    yc = read_grid(run_dir, "YC")[ysl, xsl]
    depth = read_grid(run_dir, "Depth")[ysl, xsl]
    wet = read_grid(run_dir, "hFacC")[0, ysl, xsl] > 0
    eta, _times = fit_eta_harmonic(eta_files, ysl, xsl, args.period, args.spinup_periods)
    tpxo_eta = interpolate_tpxo_eta(args.tpxo_elevation, args.constituent, xc, yc)
    wet = wet & np.isfinite(eta) & np.isfinite(tpxo_eta)
    amp = np.abs(eta)
    return (
        Case(
            name=short_name(run_dir.name),
            eta=eta,
            xc=xc,
            yc=yc,
            depth=depth,
            wet=wet,
            amp_p95=float(np.nanpercentile(amp[wet], 95.0)),
            amp_max=float(np.nanmax(amp[wet])),
        ),
        tpxo_eta,
    )


def make_tpxo_case(template: Case, tpxo_eta: np.ndarray) -> Case:
    amp = np.abs(tpxo_eta)
    return Case(
        name="TPXO9 M2",
        eta=tpxo_eta,
        xc=template.xc,
        yc=template.yc,
        depth=template.depth,
        wet=template.wet & np.isfinite(tpxo_eta),
        amp_p95=float(np.nanpercentile(amp[template.wet & np.isfinite(tpxo_eta)], 95.0)),
        amp_max=float(np.nanmax(amp[template.wet & np.isfinite(tpxo_eta)])),
    )


def plot_cases(cases: list[Case], save_path: Path, args: argparse.Namespace, shared_amp_scale: bool) -> None:
    ncols = len(cases)
    fig, axes = plt.subplots(
        2,
        ncols,
        figsize=(max(3.25 * ncols, 15.0), 7.4),
        constrained_layout=True,
        squeeze=False,
    )

    shared_vmax = np.nanpercentile(np.concatenate([np.abs(case.eta[case.wet]).ravel() for case in cases]), 99.0)
    amp_mappable = None
    phase_mappable = None
    for col, case in enumerate(cases):
        amp = np.abs(case.eta)
        phase = np.angle(case.eta)
        amp_vmax = shared_vmax if shared_amp_scale else max(np.nanpercentile(amp[case.wet], 99.0), 1e-12)

        ax = axes[0, col]
        ax.set_facecolor(args.land_color)
        ax.set_aspect("equal", adjustable="box")
        amp_mappable = ax.pcolormesh(
            case.xc,
            case.yc,
            masked(amp, case.wet),
            shading="auto",
            cmap=CMAP_AMP,
            vmin=0,
            vmax=amp_vmax,
        )
        ax.contour(case.xc, case.yc, masked(case.depth, case.wet), levels=args.bathy_contours, colors="k", linewidths=0.22, alpha=0.35)
        scale_note = f"p95={case.amp_p95:.3g} m"
        if not shared_amp_scale:
            scale_note += f", vmax={amp_vmax:.3g}"
        ax.set_title(f"{case.name}\n|η| ({scale_note})", fontsize=8)
        if col == 0:
            ax.set_ylabel("Amplitude\nLatitude")

        ax = axes[1, col]
        ax.set_facecolor(args.land_color)
        ax.set_aspect("equal", adjustable="box")
        phase_mappable = ax.pcolormesh(
            case.xc,
            case.yc,
            masked(phase, case.wet),
            shading="auto",
            cmap=CMAP_PHASE,
            vmin=-np.pi,
            vmax=np.pi,
        )
        ax.contour(case.xc, case.yc, masked(case.depth, case.wet), levels=args.bathy_contours, colors="k", linewidths=0.22, alpha=0.35)
        ax.set_title("phase", fontsize=8)
        if col == 0:
            ax.set_ylabel("Phase\nLatitude")

        for ax in axes[:, col]:
            ax.set_xlabel("Longitude")
            ax.tick_params(labelsize=7)

    label = "M2 ETAN amplitude (m), shared scale" if shared_amp_scale else "M2 ETAN amplitude (m), local per-row scale"
    fig.colorbar(amp_mappable, ax=axes[0, :], shrink=0.84, pad=0.01, label=label)
    cb = fig.colorbar(phase_mappable, ax=axes[1, :], shrink=0.84, pad=0.01, label="M2 ETAN phase (rad)")
    cb.set_ticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi])
    cb.set_ticklabels(["$-\\pi$", "$-\\pi/2$", "0", "$\\pi/2$", "$\\pi$"])

    title_scale = "shared amplitude scale" if shared_amp_scale else "local amplitude scale"
    fig.suptitle(f"Boundary-test M2 ETAN amplitude and phase vs TPXO ({title_scale})", fontsize=14)
    save_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path, dpi=args.dpi)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", nargs="*", type=Path, help="Run directories. Defaults to LR fullTPXO boundary-test runs.")
    parser.add_argument("--save-dir", type=Path, default=Path("intTides_PBC_smalldom_HR/bc_eta_tpxo_amphidrome_comparison"))
    parser.add_argument("--tpxo-elevation", type=Path, default=Path("/Users/karina/Research/Canyons/data/tides/h_tpxo9.v5a.nc"))
    parser.add_argument("--constituent", default="M2")
    parser.add_argument("--period", type=float, default=44712.0)
    parser.add_argument("--spinup-periods", type=float, default=1.0)
    parser.add_argument("--sponge-cells", type=int, default=20)
    parser.add_argument("--include-sponge", action="store_true")
    parser.add_argument("--land-color", default="tan")
    parser.add_argument("--dpi", type=int, default=170)
    parser.add_argument("--bathy-contours", nargs="+", type=float, default=[50, 100, 200, 500, 1000, 1500, 2000])
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    runs = args.runs or default_runs()
    cases: list[Case] = []
    tpxo_case: Case | None = None
    for run in runs:
        print(f"Loading {run}")
        case, tpxo_eta = load_case(run.resolve(), args)
        cases.append(case)
        if tpxo_case is None:
            tpxo_case = make_tpxo_case(case, tpxo_eta)

    if tpxo_case is None:
        raise RuntimeError("No cases loaded")
    cases.append(tpxo_case)

    shared = args.save_dir / "bc_eta_all_cases_plus_tpxo_amp_phase_shared_scale_landscape.png"
    local = args.save_dir / "bc_eta_all_cases_plus_tpxo_amp_phase_local_scale_landscape.png"
    plot_cases(cases, shared, args, shared_amp_scale=True)
    plot_cases(cases, local, args, shared_amp_scale=False)
    print(f"Wrote {shared}")
    print(f"Wrote {local}")


if __name__ == "__main__":
    main()
