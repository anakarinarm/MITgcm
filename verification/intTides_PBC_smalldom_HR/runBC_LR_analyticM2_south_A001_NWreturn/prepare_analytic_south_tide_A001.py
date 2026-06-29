#!/usr/bin/env python3
"""Generate A0=0.01 analytical south M2 forcing with north+west return."""

from __future__ import annotations

from pathlib import Path

import numpy as np


NX = 220
NY = 240
NR = 60
PERIOD = 44712.0
A0 = 0.01
TAPER_CELLS = 10
TAG = "M2_analyticSouth_A001_NWreturn_LR220x240_telescopic20"

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

    area_w = (hfac_w[:, :, 1] * drf[:, None] * dyg[None, :, 1]).sum(axis=0)
    area_s = (hfac_s[:, 1, :] * drf[:, None] * dxg[None, 1, :]).sum(axis=0)
    area_n = (hfac_s[:, -1, :] * drf[:, None] * dxg[None, -1, :]).sum(axis=0)
    wet_w = area_w > 0.0
    wet_s = area_s > 0.0
    wet_n = area_n > 0.0

    south_amp = np.where(wet_s, A0 * cosine_edge_taper(NX, TAPER_CELLS), 0.0)
    q_s = float(np.sum(south_amp * area_s))

    # OBCS sign convention in the transport diagnostic:
    # Q_in = Q_W + Q_S - Q_N.  A west return/outflow is negative west velocity,
    # represented as positive amplitude with a half-period phase lag.
    return_speed = q_s / float(area_w[wet_w].sum() + area_n[wet_n].sum())
    west_amp = np.where(wet_w, return_speed, 0.0)
    west_phase = np.where(wet_w, PERIOD / 2.0, 0.0)
    north_amp = np.where(wet_n, return_speed, 0.0)

    zero_x = np.zeros(NX)
    q_w = -float(np.sum(west_amp * area_w))
    q_n = float(np.sum(north_amp * area_n))
    residual = q_w + q_s - q_n

    write_be64(INPUT / f"OBWuam_{TAG}.bin", west_amp)
    write_be64(INPUT / f"OBWuph_{TAG}.bin", west_phase)
    write_be64(INPUT / f"OBSvam_{TAG}.bin", south_amp)
    write_be64(INPUT / f"OBSvph_{TAG}.bin", zero_x)
    write_be64(INPUT / f"OBNvam_{TAG}.bin", north_amp)
    write_be64(INPUT / f"OBNvph_{TAG}.bin", zero_x)

    print("return architecture: north + west")
    print("A0 [m/s]:", A0)
    print("area W/S/N [m2]:", float(area_w.sum()), float(area_s.sum()), float(area_n.sum()))
    print("Q_W/Q_S/Q_N [m3/s]:", q_w, q_s, q_n)
    print("residual Q_W+Q_S-Q_N [m3/s]:", residual)
    print("south max/rms [m/s]:", float(south_amp.max()), float(np.sqrt(np.mean(south_amp**2))))
    print("north return max/rms [m/s]:", float(north_amp.max()), float(np.sqrt(np.mean(north_amp**2))))
    print("west return max/rms [m/s]:", float(west_amp.max()), float(np.sqrt(np.mean(west_amp**2))))
    print("west return phase [s]:", PERIOD / 2.0)


if __name__ == "__main__":
    main()
