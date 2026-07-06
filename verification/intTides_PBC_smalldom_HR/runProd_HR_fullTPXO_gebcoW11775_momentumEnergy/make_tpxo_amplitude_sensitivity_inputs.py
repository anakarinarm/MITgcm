#!/usr/bin/env python
"""Create amplitude-scaled TPXO M2 OBCS files for sensitivity runs.

This script scales only the time-dependent tidal boundary records:

- OB[W,S,N]eta
- OB[W,S,N]u
- OB[W,S,N]v

The linear T/S boundary files are intentionally not changed.  Scaling SSH and
barotropic velocity together preserves the TPXO phase structure and the
barotropic consistency of the accepted production configuration while changing
the incident tidal amplitude.
"""

from __future__ import annotations

import array
import sys
from pathlib import Path


RUN_DIR = Path(__file__).resolve().parent
INPUT_DIR = RUN_DIR.parent / "input"

BASE_TAG = "M2_fullTPXO_HR750x720_gebcoW11775_physW11775_tel60_obeta16"

SCALES = {
    "amp050": 0.5,
    "amp200": 2.0,
}

TIDAL_PREFIXES = (
    "OBWeta",
    "OBSeta",
    "OBNeta",
    "OBWu",
    "OBWv",
    "OBSu",
    "OBSv",
    "OBNu",
    "OBNv",
)


def scaled_name(prefix: str, scale_tag: str) -> str:
    return f"{prefix}_{BASE_TAG}_{scale_tag}.bin"


def base_name(prefix: str) -> str:
    return f"{prefix}_{BASE_TAG}.bin"


def scale_file(src: Path, dst: Path, factor: float) -> None:
    data = array.array("d")
    with src.open("rb") as handle:
        data.frombytes(handle.read())
    if data.itemsize != 8:
        raise RuntimeError("Expected 8-byte native doubles")
    if len(data) == 0:
        raise ValueError(f"{src} is empty")

    # MITgcm binary input here is big-endian float64.  Convert to native
    # doubles, scale in place, then convert back to big-endian before writing.
    if sys.byteorder == "little":
        data.byteswap()
    for i, value in enumerate(data):
        data[i] = factor * value
    if sys.byteorder == "little":
        data.byteswap()
    with dst.open("wb") as handle:
        data.tofile(handle)
    print(f"{dst.name}: {factor:g} x {src.name} ({len(data)} values)")


def main() -> None:
    for scale_tag, factor in SCALES.items():
        print(f"\nCreating {scale_tag} tidal OBCS files (factor={factor:g})")
        for prefix in TIDAL_PREFIXES:
            src = INPUT_DIR / base_name(prefix)
            dst = INPUT_DIR / scaled_name(prefix, scale_tag)
            if not src.exists():
                raise FileNotFoundError(src)
            scale_file(src, dst, factor)


if __name__ == "__main__":
    main()
