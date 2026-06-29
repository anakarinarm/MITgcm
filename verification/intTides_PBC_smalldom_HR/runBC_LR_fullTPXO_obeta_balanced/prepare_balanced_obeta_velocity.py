#!/usr/bin/env python3
"""Balance the normal velocity records for the SSH+velocity OBCS test.

This script starts from the raw TPXO M2 OBCS time-record files generated for
`runBC_LR_fullTPXO_obeta`.  It leaves SSH and tangential velocities unchanged,
but applies a depth-uniform barotropic correction to the normal velocities so
the instantaneous open-boundary volume transport satisfies

    Q_W + Q_S - Q_N = 0

using MITgcm/OBCS sign convention:

    west/south positive velocities are inflow, north positive velocity is outflow.

The correction mirrors the new-method logic in pkg/obcs/obcs_balance_flow.F,
but bakes the correction into the files so `useOBCSbalance` can remain false.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


NX = 220
NY = 240
NR = 60
N_RECORDS = 16

RAW_TAG = "M2_fullTPXO_LR220x240_telescopic20_obeta16"
BALANCED_TAG = "M2_fullTPXO_balanced_LR220x240_telescopic20_obeta16"

HERE = Path(__file__).resolve().parent
EXPERIMENT = HERE.parent
GRID_DIR = EXPERIMENT / "runBC_LR_fullTPXO_obeta"
INPUT = EXPERIMENT / "input"


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


def write_be64(path: Path, array: np.ndarray) -> None:
    array.astype(">f8").tofile(path)
    print(f"wrote {path.name}  shape={array.shape}  bytes={path.stat().st_size}")


def main() -> None:
    drf = read_be32(GRID_DIR / "DRF.data", (NR,))
    dxg = read_be32(GRID_DIR / "DXG.data", (NY, NX))
    dyg = read_be32(GRID_DIR / "DYG.data", (NY, NX))
    hfac_w = read_be32(GRID_DIR / "hFacW.data", (NR, NY, NX))
    hfac_s = read_be32(GRID_DIR / "hFacS.data", (NR, NY, NX))

    area_w = hfac_w[:, :, 1] * drf[:, None] * dyg[None, :, 1]
    area_s = hfac_s[:, 1, :] * drf[:, None] * dxg[None, 1, :]
    area_n = hfac_s[:, -1, :] * drf[:, None] * dxg[None, -1, :]

    area_w_sum = float(area_w.sum())
    area_s_sum = float(area_s.sum())
    area_n_sum = float(area_n.sum())
    area_total = area_w_sum + area_s_sum + area_n_sum
    if area_total <= 0.0:
        raise ValueError("total open-boundary normal area is zero")

    obwu = read_be64(INPUT / f"OBWu_{RAW_TAG}.bin", (N_RECORDS, NR, NY))
    obsv = read_be64(INPUT / f"OBSv_{RAW_TAG}.bin", (N_RECORDS, NR, NX))
    obnv = read_be64(INPUT / f"OBNv_{RAW_TAG}.bin", (N_RECORDS, NR, NX))

    before = []
    after = []
    corrections = []

    for n in range(N_RECORDS):
        q_w = float(np.sum(obwu[n] * area_w))
        q_s = float(np.sum(obsv[n] * area_s))
        q_n = float(np.sum(obnv[n] * area_n))
        q_in = q_w + q_s - q_n
        correction = q_in / area_total

        # Same signs as OBCS_BALANCE_FLOW new method:
        # flowW=-inFlow, flowS=-inFlow, flowN=+inFlow.
        obwu[n] = obwu[n] - correction
        obsv[n] = obsv[n] - correction
        obnv[n] = obnv[n] + correction

        q_w2 = float(np.sum(obwu[n] * area_w))
        q_s2 = float(np.sum(obsv[n] * area_s))
        q_n2 = float(np.sum(obnv[n] * area_n))
        before.append(q_in)
        after.append(q_w2 + q_s2 - q_n2)
        corrections.append(correction)

    write_be64(INPUT / f"OBWu_{BALANCED_TAG}.bin", obwu)
    write_be64(INPUT / f"OBSv_{BALANCED_TAG}.bin", obsv)
    write_be64(INPUT / f"OBNv_{BALANCED_TAG}.bin", obnv)

    before = np.asarray(before)
    after = np.asarray(after)
    corrections = np.asarray(corrections)
    print("area W/S/N/total [m2]:", area_w_sum, area_s_sum, area_n_sum, area_total)
    print("raw Q_W + Q_S - Q_N max abs [m3/s]:", float(np.max(np.abs(before))))
    print("balanced residual max abs [m3/s]:", float(np.max(np.abs(after))))
    print("correction velocity max abs [m/s]:", float(np.max(np.abs(corrections))))
    print("correction velocity rms [m/s]:", float(np.sqrt(np.mean(corrections**2))))


if __name__ == "__main__":
    main()
