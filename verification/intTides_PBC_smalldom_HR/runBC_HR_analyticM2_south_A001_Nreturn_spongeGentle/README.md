# HR analytical M2 forcing: A0=0.01 m/s, north return, gentle sponge

This run is a controlled follow-up to:

- `../runBC_HR_analyticM2_south_A001_Nreturn`

The forcing and grid are unchanged. The only intended change is the momentum
sponge strength.

## Purpose

The baseline HR analytical test looked good at the south and north boundaries,
but the west diagnostic remained inward dominated:

- west bulk flux was negative at all tested offsets;
- west mode-1 incoming energy fraction increased from about `0.575` to `0.631`
  moving farther into the domain.

Because the prescribed west tide is exactly zero, this test asks whether the
west signal is partly caused by a too-abrupt/too-strong sponge rather than by
the analytical forcing itself.

## Configuration

Same as the baseline HR analytical run:

- no prescribed SSH/elevation;
- south normal M2 velocity with `A0 = 0.01 m/s`;
- 30-cell cosine taper near south corners;
- north uniform compensating return flow;
- west prescribed velocity exactly zero;
- Stevens/radiative W/S/N;
- boundary/reflection diagnostics only.

Changed sponge:

| parameter | baseline strong sponge | this gentle sponge |
|---|---:|---:|
| `spongeThickness` | 60 cells | 60 cells |
| `Urelaxobcsbound` | `1397.25 s` | `11178 s` |
| `Urelaxobcsinner` | `44712 s` | `178848 s` |
| `Vrelaxobcsbound` | `1397.25 s` | `11178 s` |
| `Vrelaxobcsinner` | `44712 s` | `178848 s` |

`useLinearSponge` remains `.FALSE.`. In the MITgcm OBCS sponge code, this keeps
a blended target that fades from the boundary value to the interior value,
which is smoother than directly relaxing every sponge cell to the boundary
value.

## Run and analyze

```bash
mpiexec -n 16 ./mitgcmuv
python3 check_boundary_reflections.py
python3 modal_characteristic_decomposition.py
```

## What to compare

Compare against:

- `../runBC_HR_analyticM2_south_A001_Nreturn`

Useful diagnostics:

- final ETAN min/max/mean/RMS;
- west `outward_flux_W` at offsets 65, 70, 80, 100;
- west mode-1 outward/inward characteristic energy;
- west mode-1 incoming energy fraction;
- south mode-1 outward/inward energy to ensure the source tide was not degraded.

Interpretation:

- If the west inward signal weakens, the strong sponge was probably too abrupt
  or too reflective.
- If the west signal is unchanged or worse, the next suspect is the west-side
  bathymetry / oblique internal-wave field / diagnostic geometry rather than
  sponge strength alone.

## Result

The run completed normally. STDERR only contained the harmless retired
`EXACT_CONSERV` warning.

Final ETAN snapshot:

| run | ETAN min m | ETAN max m | ETAN mean m | ETAN RMS m |
|---|---:|---:|---:|---:|
| baseline strong sponge | -0.03718 | 0.00819 | -0.00674 | 0.01218 |
| gentle sponge | -0.03607 | 0.00776 | -0.00675 | 0.01182 |

The free surface is essentially unchanged and very slightly cleaner by RMS, so
the test does not indicate a barotropic mass-balance problem.

Boundary/reflection comparison:

| run | boundary | offset cells | bulk outward flux MW | progressive index |
|---|---|---:|---:|---:|
| baseline strong sponge | W | 65 | -2.03 | -0.093 |
| baseline strong sponge | W | 70 | -2.54 | -0.113 |
| baseline strong sponge | W | 80 | -3.55 | -0.149 |
| baseline strong sponge | W | 100 | -4.74 | -0.184 |
| gentle sponge | W | 65 | -3.17 | -0.145 |
| gentle sponge | W | 70 | -3.61 | -0.161 |
| gentle sponge | W | 80 | -4.43 | -0.190 |
| gentle sponge | W | 100 | -5.44 | -0.219 |
| baseline strong sponge | S | 65 | 10.37 | 0.489 |
| gentle sponge | S | 65 | 10.82 | 0.529 |
| baseline strong sponge | N | 65 | 1.70 | 0.338 |
| gentle sponge | N | 65 | 1.62 | 0.355 |

Mode-1 characteristic comparison:

| run | boundary | offset cells | mode-1 outward MW | mode-1 inward MW | incoming fraction |
|---|---|---:|---:|---:|---:|
| baseline strong sponge | W | 65 | 4.16 | 5.63 | 0.575 |
| baseline strong sponge | W | 100 | 4.39 | 7.51 | 0.631 |
| gentle sponge | W | 65 | 3.68 | 5.94 | 0.617 |
| gentle sponge | W | 100 | 3.76 | 7.54 | 0.667 |
| baseline strong sponge | S | 65 | 9.53 | 1.60 | 0.144 |
| gentle sponge | S | 65 | 9.35 | 1.06 | 0.101 |
| baseline strong sponge | N | 65 | 0.85 | 0.36 | 0.296 |
| gentle sponge | N | 65 | 0.72 | 0.29 | 0.288 |

Interpretation:

- Gentle sponge did **not** improve the west boundary. West bulk flux became
  more negative and west mode-1 incoming fraction increased.
- South became slightly cleaner by the modal diagnostic, with lower incoming
  mode-1 fraction.
- North stayed modest.
- This points away from the idea that the baseline sponge was too abrupt or too
  reflective at the west boundary. Instead, the west boundary likely needs more
  effective absorption, more targeted west-side treatment, or a diagnostic that
  better handles oblique/topographic wave structure.

Decision:

- Do not use the gentle-sponge variant as the production baseline.
- Prefer the original strong-sponge HR analytical run unless a better
  west-focused test improves the west diagnostics.

