#!/usr/bin/env python3
"""Create a bathymetry variant suitable for prescribed-SSH OBCS tests.

The SSH+velocity OBCS experiment uses MITgcm's nonlinear free-surface/r*
formulation so that OB*etaFile can be prescribed.  The original LR telescopic
bathymetry contains sub-meter wet slivers.  Those cells are fine-ish for a
linear free-surface velocity-forced test, but not for a prescribed SSH tide of
O(0.4 m): the free surface can remove too much of the water column and trigger
CALC_R_STAR.

For this boundary-condition architecture test, we are not interested in
tidal-flat/shallow-sliver dynamics, so keep all existing wet cells wet but
deepen any wet cell shallower than MIN_WET_DEPTH to that minimum depth.
"""

from __future__ import annotations

import struct
from pathlib import Path


NX = 220
NY = 240
MIN_WET_DEPTH = 10.0  # m

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "input"

SOURCE = INPUT / "smalldom_LR220x240_telescopic20_bat.bin"
TARGET = INPUT / "smalldom_LR220x240_telescopic20_bat_min10m_obeta.bin"


def main() -> None:
    raw = SOURCE.read_bytes()
    expected_bytes = NX * NY * 8
    if len(raw) != expected_bytes:
        raise ValueError(f"{SOURCE}: expected {expected_bytes} bytes, got {len(raw)}")

    bathy = list(struct.unpack(f">{NX * NY}d", raw))
    out = bathy.copy()

    wet_depth_before = []
    wet_depth_after = []
    shallow_count = 0
    for index, value in enumerate(bathy):
        if value < 0.0:
            depth = -value
            wet_depth_before.append(depth)
            if depth < MIN_WET_DEPTH:
                out[index] = -MIN_WET_DEPTH
                shallow_count += 1
            wet_depth_after.append(-out[index])

    TARGET.write_bytes(struct.pack(f">{NX * NY}d", *out))

    print(f"source: {SOURCE}")
    print(f"target: {TARGET}")
    print(f"wet cells: {len(wet_depth_before)}")
    print(f"modified wet cells shallower than {MIN_WET_DEPTH:g} m: {shallow_count}")
    print(
        "minimum wet depth before/after: "
        f"{min(wet_depth_before):.6g} / {min(wet_depth_after):.6g} m"
    )


if __name__ == "__main__":
    main()
