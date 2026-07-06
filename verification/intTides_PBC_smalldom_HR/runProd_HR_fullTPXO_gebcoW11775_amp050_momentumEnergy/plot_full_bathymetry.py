#!/usr/bin/env python3
"""Plot full GEBCO W117.75 HR bathymetry, including the W/S/N sponge."""

from __future__ import annotations

from pathlib import Path

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import cmocean

    CMAP = cmocean.cm.deep
except Exception:  # pragma: no cover
    CMAP = "viridis_r"


RUN_DIR = Path(__file__).resolve().parent
EXPERIMENT = RUN_DIR.parent
INPUT = EXPERIMENT / "input"

TAG = "HR750x720_gebcoW11775_physW11775_tel60"
NX = 750
NY = 720
SPONGE = 60
XG_ORIGIN = 241.982924784919
YG_ORIGIN = 31.022894044779


def read_be64(path: Path, shape: tuple[int, ...]) -> np.ndarray:
    values = np.fromfile(path, dtype=">f8")
    expected = int(np.prod(shape))
    if values.size != expected:
        raise ValueError(f"{path}: expected {expected} values, got {values.size}")
    return values.reshape(shape)


def edges(origin: float, widths: np.ndarray) -> np.ndarray:
    return origin + np.r_[0.0, np.cumsum(widths)]


def main() -> None:
    dx = read_be64(INPUT / f"smalldom_{TAG}_dx.bin", (NX,))
    dy = read_be64(INPUT / f"smalldom_{TAG}_dy.bin", (NY,))
    bathy = read_be64(INPUT / f"smalldom_{TAG}_bat_min10m_obeta.bin", (NY, NX))

    x_edges = edges(XG_ORIGIN, dx)
    y_edges = edges(YG_ORIGIN, dy)
    x_plot = np.where(x_edges > 180.0, x_edges - 360.0, x_edges)
    depth = np.ma.masked_where(bathy >= 0.0, -bathy)

    fig, ax = plt.subplots(figsize=(11.5, 8.5), constrained_layout=True)
    ax.set_facecolor("tan")
    mesh = ax.pcolormesh(
        x_plot,
        y_edges,
        depth,
        shading="auto",
        cmap=CMAP,
        vmin=0.0,
        vmax=float(depth.max()),
    )
    cbar = fig.colorbar(mesh, ax=ax, shrink=0.88)
    cbar.set_label("Depth (m)")

    west_inner = x_plot[SPONGE]
    south_inner = y_edges[SPONGE]
    north_inner = y_edges[NY - SPONGE]
    ax.axvline(west_inner, color="k", lw=1.4, ls="--", alpha=0.8, label="inner sponge edge")
    ax.axhline(south_inner, color="k", lw=1.4, ls="--", alpha=0.8)
    ax.axhline(north_inner, color="k", lw=1.4, ls="--", alpha=0.8)

    ax.plot(-117.75, 31.25, marker="*", ms=12, color="magenta", mec="k", mew=0.7)
    ax.text(
        -117.75,
        31.25,
        " physical SW corner",
        color="magenta",
        ha="left",
        va="bottom",
        fontsize=9,
        weight="bold",
    )

    ax.set_title("HR GEBCO W117.75 bathymetry with W/S/N sponge")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(loc="upper left", frameon=True)

    out = RUN_DIR / "full_bathymetry_with_sponge.png"
    fig.savefig(out, dpi=220)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
