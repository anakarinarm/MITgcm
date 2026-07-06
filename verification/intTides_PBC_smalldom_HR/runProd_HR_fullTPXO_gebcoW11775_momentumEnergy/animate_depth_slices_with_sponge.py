#!/usr/bin/env python3
"""Animate depth slices including the W/S/N sponge buffer.

This is a convenience wrapper around ``animate_depth_slices.py``.  It keeps the
same 3-panel layout and default depths, but sets:

```text
--sponge-cells 0
--mark-sponge-boundary
--save-dir output/animations_with_sponge
```

The dashed lines mark the edge of the original physical 600 x 600 domain:

- west physical edge: i = 60;
- south physical edge: j = 60;
- north physical edge: j = 660.
"""

from __future__ import annotations

import sys
from pathlib import Path

from animate_depth_slices import main


def has_option(*names: str) -> bool:
    return any(arg in names or any(arg.startswith(name + "=") for name in names) for arg in sys.argv[1:])


if __name__ == "__main__":
    run_dir = Path(__file__).resolve().parent

    if not has_option("--sponge-cells"):
        sys.argv.extend(["--sponge-cells", "0"])

    if not has_option("--mark-sponge-boundary"):
        sys.argv.append("--mark-sponge-boundary")

    if not has_option("--save-dir"):
        sys.argv.extend(["--save-dir", str(run_dir / "output" / "animations_with_sponge")])

    main()
