#!/usr/bin/env python3
"""Generate analytical south-entering M2 tidal boundary forcing.

Option A:

- south boundary receives a prescribed normal M2 velocity with amplitude A0
  and a cosine taper near the corners;
- north boundary receives a spatially uniform compensating normal velocity so
  instantaneous complex barotropic transport closes;
- west boundary receives zero prescribed normal velocity;
- no SSH is prescribed.

The files are MITgcm OBCS tidal amplitude/phase files.  Phase is stored as a
time lag in seconds over the M2 period; this experiment uses zero phase
everywhere for the imposed south and north velocities.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


NX = 220
NY = 240
NR = 60
PERIOD = 44712.0
A0 = 0.02
TAPER_CELLS = 10
TAG = "M2_analyticSouth_A002_LR220x240_telescopic20"

HERE = Path(__file__).resolve().parent
EXPERIMENT = HERE.parent
INPUT = EXPERIMENT / "input"
GRID_DIR = EXPERIMENT / "runBC_LR_fullTPXO_balanced"


def read_be32(path: Path, shape: tuple[int, ...]) -> np.ndarray:
    values = np.fromfile(path, dtype=">f4")
    expected = int(np.prod(shape))
    if values.size != expected:
        raise ValueError(f"{path}: expected {expected} values, got {values.size}")
    return values.reshape(shape).astype(np.float64)


def write_be64(path: Path, array: np.ndarray) -> None:
    array.astype(">f8").tofile(path)
    print(f"wrote {path.name}: n={array.size}, bytes={path.stat().st_size}")


def cosine_edge_taper(n: int, n_taper: int) -> np.ndarray:
    taper = np.ones(n, dtype=np.float64)
    if n_taper <= 0:
        return taper
    ramp = 0.5 * (1.0 - np.cos(np.pi * np.arange(n_taper) / n_taper))
    taper[:n_taper] = ramp
    taper[-n_taper:] = ramp[::-1]
    return taper


def main() -> None:
    drf = read_be32(GRID_DIR / "DRF.data", (NR,))
    dxg = read_be32(GRID_DIR / "DXG.data", (NY, NX))
    dyg = read_be32(GRID_DIR / "DYG.data", (NY, NX))
    hfac_w = read_be32(GRID_DIR / "hFacW.data", (NR, NY, NX))
    hfac_s = read_be32(GRID_DIR / "hFacS.data", (NR, NY, NX))

    # Depth-integrated normal area per horizontal open-boundary point.
    area_w = (hfac_w[:, :, 1] * drf[:, None] * dyg[None, :, 1]).sum(axis=0)
    area_s = (hfac_s[:, 1, :] * drf[:, None] * dxg[None, 1, :]).sum(axis=0)
    area_n = (hfac_s[:, -1, :] * drf[:, None] * dxg[None, -1, :]).sum(axis=0)

    wet_s = area_s > 0.0
    wet_n = area_n > 0.0

    south_amp = A0 * cosine_edge_taper(NX, TAPER_CELLS)
    south_amp = np.where(wet_s, south_amp, 0.0)
    q_s = float(np.sum(south_amp * area_s))

    north_amp_uniform = q_s / float(area_n[wet_n].sum())
    north_amp = np.where(wet_n, north_amp_uniform, 0.0)

    west_amp = np.zeros(NY, dtype=np.float64)
    zero_w = np.zeros(NY, dtype=np.float64)
    zero_x = np.zeros(NX, dtype=np.float64)

    q_w = float(np.sum(west_amp * area_w))
    q_n = float(np.sum(north_amp * area_n))
    residual = q_w + q_s - q_n

    write_be64(INPUT / f"OBWuam_{TAG}.bin", west_amp)
    write_be64(INPUT / f"OBWuph_{TAG}.bin", zero_w)
    write_be64(INPUT / f"OBSvam_{TAG}.bin", south_amp)
    write_be64(INPUT / f"OBSvph_{TAG}.bin", zero_x)
    write_be64(INPUT / f"OBNvam_{TAG}.bin", north_amp)
    write_be64(INPUT / f"OBNvph_{TAG}.bin", zero_x)

    print("A0 [m/s]:", A0)
    print("taper cells:", TAPER_CELLS)
    print("area W/S/N [m2]:", float(area_w.sum()), float(area_s.sum()), float(area_n.sum()))
    print("Q_W/Q_S/Q_N [m3/s]:", q_w, q_s, q_n)
    print("Q_W + Q_S - Q_N residual [m3/s]:", residual)
    print("south amp min/max/rms [m/s]:", float(south_amp.min()), float(south_amp.max()), float(np.sqrt(np.mean(south_amp**2))))
    print("north uniform return amp [m/s]:", float(north_amp_uniform))
    print("west amp max [m/s]:", float(west_amp.max()))


if __name__ == "__main__":
    main()
