#!/usr/bin/env python3
"""Build LR no-west balanced M2 tidal boundary files.

Starting point is the LR full-TPXO balanced harmonic forcing.  This script sets
the west normal velocity tide to zero, then applies a depth-uniform complex
barotropic correction to the south and north normal velocity tides so

    Q_W + Q_S - Q_N = 0

with Q_W = 0.  The output remains in MITgcm OBCS amplitude/phase files.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


NX = 220
NY = 240
NR = 60
PERIOD = 44712.0
OMEGA = 2.0 * np.pi / PERIOD

TAG_IN = "M2_fullTPXO_balanced_LR220x240_telescopic20"
TAG_OUT = "M2_noW_balanced_LR220x240_telescopic20"

HERE = Path(__file__).resolve().parent
EXPERIMENT = HERE.parent
INPUT = EXPERIMENT / "input"
GRID_DIR = EXPERIMENT / "runBC_LR_fullTPXO_balanced"


def read_be64(path: Path, shape: tuple[int, ...]) -> np.ndarray:
    values = np.fromfile(path, dtype=">f8")
    expected = int(np.prod(shape))
    if values.size != expected:
        raise ValueError(f"{path}: expected {expected} values, got {values.size}")
    return values.reshape(shape).astype(np.float64)


def read_be32(path: Path, shape: tuple[int, ...]) -> np.ndarray:
    values = np.fromfile(path, dtype=">f4")
    expected = int(np.prod(shape))
    if values.size != expected:
        raise ValueError(f"{path}: expected {expected} values, got {values.size}")
    return values.reshape(shape).astype(np.float64)


def amp_phase_to_complex(amp: np.ndarray, phase_s: np.ndarray) -> np.ndarray:
    # OBCS tidal phase files are phase lags in seconds over the tidal period.
    return amp * np.exp(-1j * OMEGA * phase_s)


def complex_to_amp_phase(z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    amp = np.abs(z)
    phase = (-np.angle(z) / OMEGA) % PERIOD
    phase = np.where(amp > 0.0, phase, 0.0)
    return amp, phase


def write_be64(path: Path, array: np.ndarray) -> None:
    array.astype(">f8").tofile(path)
    print(f"wrote {path.name}: n={array.size}, bytes={path.stat().st_size}")


def main() -> None:
    drf = read_be32(GRID_DIR / "DRF.data", (NR,))
    dxg = read_be32(GRID_DIR / "DXG.data", (NY, NX))
    dyg = read_be32(GRID_DIR / "DYG.data", (NY, NX))
    hfac_w = read_be32(GRID_DIR / "hFacW.data", (NR, NY, NX))
    hfac_s = read_be32(GRID_DIR / "hFacS.data", (NR, NY, NX))

    area_w = (hfac_w[:, :, 1] * drf[:, None] * dyg[None, :, 1]).sum(axis=0)
    area_s = (hfac_s[:, 1, :] * drf[:, None] * dxg[None, 1, :]).sum(axis=0)
    area_n = (hfac_s[:, -1, :] * drf[:, None] * dxg[None, -1, :]).sum(axis=0)

    aw = read_be64(INPUT / f"OBWuam_{TAG_IN}.bin", (NY,))
    pw = read_be64(INPUT / f"OBWuph_{TAG_IN}.bin", (NY,))
    a_s = read_be64(INPUT / f"OBSvam_{TAG_IN}.bin", (NX,))
    p_s = read_be64(INPUT / f"OBSvph_{TAG_IN}.bin", (NX,))
    a_n = read_be64(INPUT / f"OBNvam_{TAG_IN}.bin", (NX,))
    p_n = read_be64(INPUT / f"OBNvph_{TAG_IN}.bin", (NX,))

    cw = amp_phase_to_complex(aw, pw)
    cs = amp_phase_to_complex(a_s, p_s)
    cn = amp_phase_to_complex(a_n, p_n)

    q_before = np.sum(cw * area_w) + np.sum(cs * area_s) - np.sum(cn * area_n)

    cw_new = np.zeros_like(cw)
    q_no_w = np.sum(cs * area_s) - np.sum(cn * area_n)
    correction = q_no_w / (area_s.sum() + area_n.sum())
    cs_new = cs - correction
    cn_new = cn + correction

    q_after = np.sum(cw_new * area_w) + np.sum(cs_new * area_s) - np.sum(cn_new * area_n)

    aw_new, pw_new = complex_to_amp_phase(cw_new)
    as_new, ps_new = complex_to_amp_phase(cs_new)
    an_new, pn_new = complex_to_amp_phase(cn_new)

    write_be64(INPUT / f"OBWuam_{TAG_OUT}.bin", aw_new)
    write_be64(INPUT / f"OBWuph_{TAG_OUT}.bin", pw_new)
    write_be64(INPUT / f"OBSvam_{TAG_OUT}.bin", as_new)
    write_be64(INPUT / f"OBSvph_{TAG_OUT}.bin", ps_new)
    write_be64(INPUT / f"OBNvam_{TAG_OUT}.bin", an_new)
    write_be64(INPUT / f"OBNvph_{TAG_OUT}.bin", pn_new)

    print("area W/S/N [m2]:", float(area_w.sum()), float(area_s.sum()), float(area_n.sum()))
    print("full balanced residual |Q| [m3/s]:", float(abs(q_before)))
    print("no-west raw |Q| [m3/s]:", float(abs(q_no_w)))
    print("S/N correction amplitude [m/s]:", float(abs(correction)))
    print("no-west balanced residual |Q| [m3/s]:", float(abs(q_after)))
    print("west max amplitude [m/s]:", float(aw_new.max()))
    print("south max amplitude [m/s]:", float(as_new.max()))
    print("north max amplitude [m/s]:", float(an_new.max()))


if __name__ == "__main__":
    main()
