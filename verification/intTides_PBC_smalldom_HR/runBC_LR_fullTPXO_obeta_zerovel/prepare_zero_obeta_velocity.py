#!/usr/bin/env python3
"""Generate exact-zero prescribed velocity files for the SSH-only OBCS test."""

from __future__ import annotations

from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
EXPERIMENT = HERE.parent
INPUT = EXPERIMENT / "input"

BASE_FILES = {
    "OBWu": "OBWu_M2_fullTPXO_balanced_LR220x240_telescopic20_obeta16.bin",
    "OBWv": "OBWv_M2_fullTPXO_LR220x240_telescopic20_obeta16.bin",
    "OBSu": "OBSu_M2_fullTPXO_LR220x240_telescopic20_obeta16.bin",
    "OBSv": "OBSv_M2_fullTPXO_balanced_LR220x240_telescopic20_obeta16.bin",
    "OBNu": "OBNu_M2_fullTPXO_LR220x240_telescopic20_obeta16.bin",
    "OBNv": "OBNv_M2_fullTPXO_balanced_LR220x240_telescopic20_obeta16.bin",
}


def main() -> None:
    for name, base in BASE_FILES.items():
        src = INPUT / base
        values = np.fromfile(src, dtype=">f8")
        if values.size == 0:
            raise ValueError(f"{src} is empty")
        out_name = base.replace(
            "_LR220x240_telescopic20_obeta16.bin",
            "_zerovel_LR220x240_telescopic20_obeta16.bin",
        )
        out = INPUT / out_name
        np.zeros(values.size, dtype=">f8").tofile(out)
        print(f"{name}: wrote {out.name} n={values.size} bytes={out.stat().st_size}")


if __name__ == "__main__":
    main()
