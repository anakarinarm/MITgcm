# LR analytical M2 forcing: south A0=0.02 m/s, north return

This is the first controlled analytical-forcing experiment.

It uses the no-SSH/radiative/strong-sponge architecture, but replaces the TPXO
boundary velocity files with idealized M2 harmonic forcing:

- south boundary: prescribed normal velocity, `A0 = 0.02 m/s`;
- south amplitude is cosine-tapered over 10 cells near both corners;
- north boundary: spatially uniform compensating return velocity;
- west boundary: zero prescribed normal velocity;
- no prescribed SSH/elevation;
- Stevens/radiative behavior remains enabled on W/S/N;
- 20-cell strong momentum sponge is retained.

The goal is not to reproduce TPXO exactly.  The goal is to create a controlled
barotropic tide that can generate internal tides without importing a full,
coarse, externally phase-locked tidal solution into this small domain.

## Generated files

Created by:

```bash
python3 prepare_analytic_south_tide.py
```

Output files:

- `../input/OBWuam_M2_analyticSouth_A002_LR220x240_telescopic20.bin`
- `../input/OBWuph_M2_analyticSouth_A002_LR220x240_telescopic20.bin`
- `../input/OBSvam_M2_analyticSouth_A002_LR220x240_telescopic20.bin`
- `../input/OBSvph_M2_analyticSouth_A002_LR220x240_telescopic20.bin`
- `../input/OBNvam_M2_analyticSouth_A002_LR220x240_telescopic20.bin`
- `../input/OBNvph_M2_analyticSouth_A002_LR220x240_telescopic20.bin`

Generation diagnostics:

- south max amplitude: `0.020 m/s`
- south RMS amplitude: `0.01937 m/s`
- north uniform return amplitude: `0.04290 m/s`
- west max amplitude: `0.0 m/s`
- transport residual `Q_W + Q_S - Q_N`: `0.0 m3/s`

Important caveat: the north open-section area is much smaller than the south
open-section area, so the north return velocity is larger than the imposed
south amplitude.  This is still useful as a controlled first diagnostic, but it
is best interpreted as a south-forced/north-return barotropic tide, not a pure
one-way Kelvin wave.

## Run and analyze

```bash
mpiexec -n 16 ./mitgcmuv
python3 check_boundary_reflections.py
python3 modal_characteristic_decomposition.py
```

Compare against:

- `../runBC_LR_noSSH_noW_strongSponge`
- `../runBC_LR_noSSH_sponge_strong`
- `../runBC_LR_recommended_noSSH_balanced_radiative`

Useful first questions:

- Does ETAN remain small?
- Does the south boundary produce an outgoing internal tide with cleaner modal
  characteristics than TPXO-derived forcing?
- Does west remain quieter when west prescribed velocity is zero?
- Does the north return flow contaminate the north boundary too strongly?

## Result

The run completed normally; STDERR only contained the harmless retired
`EXACT_CONSERV` warning.

ETAN was the cleanest of the tested boundary-forcing family:

| run | final ETAN min/max | mean | RMS |
|---|---:|---:|---:|
| analytic south A002 | `-0.0237 / +0.0141 m` | `+0.00018 m` | `0.00589 m` |
| noW + strong sponge | `-0.0315 / +0.0039 m` | `-0.00784 m` | `0.01153 m` |
| TPXO full + strong sponge | `-0.0208 / +0.0195 m` | `+0.00035 m` | `0.00902 m` |

At offset 22 cells, mode-1 characteristic energy was:

| run | W mode 1 | S mode 1 | N mode 1 |
|---|---:|---:|---:|
| analytic south A002 | `17.46 / 17.22 MW`, frac `0.497` | `44.49 / 13.21 MW`, frac `0.229` | `1.68 / 1.03 MW`, frac `0.381` |
| noW + strong sponge | `4.00 / 3.98 MW`, frac `0.498` | `10.39 / 2.91 MW`, frac `0.219` | `0.46 / 0.22 MW`, frac `0.326` |
| TPXO full + strong sponge | `7.26 / 7.63 MW`, frac `0.512` | `20.58 / 4.58 MW`, frac `0.182` | `0.38 / 0.22 MW`, frac `0.362` |

West mode 1 remained nearly balanced at the near-boundary section but became
inward dominated farther in:

| offset | W mode 1 outward/inward | incoming fraction |
|---:|---:|---:|
| 22 | `17.46 / 17.22 MW` | `0.497` |
| 23 | `17.28 / 17.57 MW` | `0.504` |
| 27 | `16.57 / 19.00 MW` | `0.534` |
| 33 | `15.71 / 22.11 MW` | `0.585` |

Interpretation:

- The analytical forcing is barotropically clean: ETAN RMS is lower than the
  TPXO-derived cases and there is no mean ETAN drift.
- The south boundary produces a strong, clean outgoing internal tide.
- The experiment is probably too energetic for a first production candidate:
  the south mode-1 outgoing energy is about twice the TPXO full+strong case.
- The large energy level is consistent with the imposed transport: south
  `A0=0.02 m/s` requires a `0.0429 m/s` uniform north return because the north
  open area is much smaller than the south open area.

Recommended next analytical tests:

1. `A0=0.01 m/s` with the same north-return architecture.  This should reduce
   energy roughly by a factor of four and is likely closer to the TPXO-energy
   range.
2. A return-flow variant that spreads the compensation over both north and
   west while keeping west much weaker than TPXO, to reduce the large north
   return velocity.

## A001 follow-up comparison

Two `A0=0.01 m/s` follow-up tests were run:

- `runBC_LR_analyticM2_south_A001_Nreturn`
- `runBC_LR_analyticM2_south_A001_NWreturn`

Summary at offset 22:

| run | W mode 1 | S mode 1 | N mode 1 | ETAN RMS |
|---|---:|---:|---:|---:|
| A002 N-return | `17.46 / 17.22 MW` | `44.49 / 13.21 MW` | `1.68 / 1.03 MW` | `0.00589 m` |
| A001 N-return | `4.38 / 4.55 MW` | `11.41 / 3.72 MW` | `0.41 / 0.26 MW` | `0.00289 m` |
| A001 NW-return | `1.76 / 1.95 MW` | `5.20 / 1.39 MW` | `0.089 / 0.086 MW` | `0.00185 m` |

The `A001 N-return` case is the best current analytical candidate if we want a
south internal-tide energy scale comparable to the noW TPXO test.  The
`A001 NW-return` case is the cleanest barotropically but may be too weak unless
the amplitude is increased.
