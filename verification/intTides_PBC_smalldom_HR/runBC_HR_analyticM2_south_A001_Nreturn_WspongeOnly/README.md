# HR analytical M2 forcing: A0=0.01 m/s, north return, west-only strong sponge

This run is a controlled follow-up to:

- `../runBC_HR_analyticM2_south_A001_Nreturn`
- `../runBC_HR_analyticM2_south_A001_Nreturn_spongeGentle`

The analytical forcing, grid, diagnostics, and strong sponge timescales are
unchanged from the baseline. The only intended change is that the momentum
sponge is active on the west boundary only.

## Purpose

The baseline HR analytical run had good south/north behavior but an
inward-dominated west diagnostic. The gentle-sponge test made the west signal
worse, suggesting the west side needs more effective absorption, not weaker
absorption.

This test asks whether the strong sponge should be focused only on the west
boundary, leaving the south source boundary and north return boundary to behave
more radiatively through Stevens without local sponge damping.

## Configuration

Same as the baseline HR analytical run:

- no prescribed SSH/elevation;
- south normal M2 velocity with `A0 = 0.01 m/s`;
- 30-cell cosine taper near south corners;
- north uniform compensating return flow;
- west prescribed velocity exactly zero;
- Stevens/radiative W/S/N;
- boundary/reflection diagnostics only.

Changed OBCS sponge side switches:

| parameter | baseline strong sponge | this test |
|---|---|---|
| `OBCSsponge_W` | `.TRUE.` | `.TRUE.` |
| `OBCSsponge_S` | `.TRUE.` | `.FALSE.` |
| `OBCSsponge_N` | `.TRUE.` | `.FALSE.` |
| `OBCSsponge_E` | `.FALSE.` | `.FALSE.` |

Sponge timescales remain the strong baseline values:

| parameter | value |
|---|---:|
| `spongeThickness` | 60 cells |
| `Urelaxobcsbound` | `1397.25 s` |
| `Urelaxobcsinner` | `44712 s` |
| `Vrelaxobcsbound` | `1397.25 s` |
| `Vrelaxobcsinner` | `44712 s` |

## Run and analyze

```bash
mpiexec -n 16 ./mitgcmuv
python3 check_boundary_reflections.py
python3 modal_characteristic_decomposition.py
```

## What to compare

Compare against:

- `../runBC_HR_analyticM2_south_A001_Nreturn`
- `../runBC_HR_analyticM2_south_A001_Nreturn_spongeGentle`

Useful diagnostics:

- final ETAN min/max/mean/RMS;
- west `outward_flux_W` at offsets 65, 70, 80, 100;
- west mode-1 outward/inward characteristic energy;
- west mode-1 incoming energy fraction;
- south mode-1 outward/inward energy and incoming fraction;
- north mode-1 energy, to check that removing north sponge does not pollute
  the return boundary.

Interpretation:

- If west improves and south/north remain clean, this is a strong candidate
  boundary setup.
- If west improves but north gets noisy, the north return may need sponge or a
  different return-flow design.
- If west does not improve, the issue is likely not controlled by OBCS sponge
  side switches alone.

