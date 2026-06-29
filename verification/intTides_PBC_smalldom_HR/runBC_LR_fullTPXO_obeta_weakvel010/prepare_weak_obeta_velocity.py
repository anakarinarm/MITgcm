#!/usr/bin/env python3
"""Generate weak prescribed-velocity records for the SSH-clamped OBCS tests.

This is the same helper included in the alpha=0.25 run directory.  It writes
both weak-velocity families, alpha=0.25 and alpha=0.10, into `../input/`.
"""

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


def write_scaled(alpha: float, tag: str) -> None:
    print(f"alpha={alpha:g} tag=weakvel{tag}")
    for name, base in BASE_FILES.items():
        src = INPUT / base
        values = np.fromfile(src, dtype=">f8")
        if values.size == 0:
            raise ValueError(f"{src} is empty")
        if not np.isfinite(values).all():
            raise ValueError(f"{src} contains non-finite values")

        out_name = base.replace(
            "_LR220x240_telescopic20_obeta16.bin",
            f"_weakvel{tag}_LR220x240_telescopic20_obeta16.bin",
        )
        out = INPUT / out_name
        (values * alpha).astype(">f8").tofile(out)
        print(
            f"  {name}: {out.name} n={values.size} "
            f"maxabs={np.max(np.abs(values)) * alpha:.6e}"
        )


def main() -> None:
    write_scaled(0.25, "025")
    write_scaled(0.10, "010")


if __name__ == "__main__":
    main()
