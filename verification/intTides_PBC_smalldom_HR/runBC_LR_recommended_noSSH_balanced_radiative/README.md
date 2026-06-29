# Recommended LR boundary setup: no SSH, balanced/radiative velocity

This is the LR configuration branch to carry forward after the boundary
architecture tests.

The important choice is deliberate:

- do **not** prescribe `OB*etaFile`;
- prescribe balanced W/S/N M2 normal velocity amplitudes/phases;
- keep Stevens/radiative behavior on W/S/N;
- keep a 20-cell momentum-only sponge in the telescopic buffer;
- keep T/S boundary profiles available for inflow/radiative tracer handling.

## Why this branch

The prescribed-SSH experiments showed that full-amplitude TPXO SSH clamping
creates a large near-standing modal characteristic signal at the open
boundaries, even when velocity is exactly zero.

Key comparison at offset 22 cells:

| test | W mode-1 outward/inward | S mode-1 outward/inward | N mode-1 outward/inward |
|---|---:|---:|---:|
| SSH + zero velocity | 182.47 / 182.47 MW | 84.98 / 84.98 MW | 53.12 / 53.12 MW |
| SSH + 10% velocity | 182.68 / 182.62 MW | 85.02 / 85.18 MW | 53.21 / 53.26 MW |
| no-SSH balanced/radiative | 7.11 / 7.75 MW | 19.84 / 4.07 MW | 0.38 / 0.22 MW |

The zero-velocity SSH test proved that the large standing signal comes from
SSH clamping itself, not from velocity imbalance.

ETAN was also much calmer without SSH clamping:

- SSH + zero velocity final ETAN: min `-0.425 m`, RMS `0.0407 m`
- no-SSH balanced/radiative final ETAN: min `-0.0206 m`, RMS `0.0090 m`

So this branch is the best current production candidate.

## Configuration

- Grid: `220 x 240 x 60`
- Physical domain: `200 x 200`, coarsened by 3x from HR
- Buffer: 20 telescopic cells west/south/north
- Eastern side: closed
- Boundary forcing: W/S/N M2 normal velocity only, volume-balanced in the
  precomputed amplitude/phase files
- No prescribed SSH/elevation
- Radiation: Stevens enabled on W/S/N
- Sponge: momentum-only, 20 cells
- Horizontal viscosity: grid-aware `viscAhGrid=0.05`
- Vertical viscosity: weak background `viscAz=1e-4`
- Bottom drag: quadratic, `0.002`
- Diagnostics: boundary-wave fields for reflection/modal analysis

## Run and analyze

```bash
mpiexec -n 16 ./mitgcmuv
python3 check_boundary_reflections.py
python3 modal_characteristic_decomposition.py
```

## Current known behavior

This clean recommended directory has now been run and reproduces the source
run, `runBC_LR_fullTPXO_balanced`, exactly.

The run completed normally; STDERR only contained the harmless retired
`EXACT_CONSERV` warning.

At offset 22 cells:

- bulk fluxes: W `+0.84 MW`, S `+19.68 MW`, N `+0.62 MW`
- W mode 1: `7.11 / 7.75 MW`, incoming fraction `0.522`
- S mode 1: `19.84 / 4.07 MW`, incoming fraction `0.170`
- N mode 1: `0.38 / 0.22 MW`, incoming fraction `0.372`

Final saved ETAN:

- min/max: `-0.0206 / +0.0193 m`
- mean: `+0.00036 m`
- RMS: `0.00895 m`

For comparison, the SSH-only zero-velocity test had final ETAN min
`-0.425 m` and RMS `0.0407 m`, plus much larger standing modal energies.

This is not perfect — the west side still deserves sponge/radiation tuning —
but it is physically cleaner than full-boundary SSH clamping and much less
contaminated by an imposed standing boundary signal.

## Recommended next tuning tests

Start from this directory and vary only one thing at a time:

1. Sponge strength/shape:
   - current: boundary relaxation one M2/16 record interval, inner relaxation
     two M2 periods;
   - next mild test: keep `spongeThickness=20`, make inner relaxation weaker
     or boundary relaxation slightly stronger.
2. Sponge width:
   - test `30` cells if/when the LR grid is extended accordingly;
   - avoid adding new bathymetric features in the sponge.
3. Forcing architecture:
   - if west remains problematic, test south-dominant/no-west velocity forcing
     within the no-SSH/radiative family, not the SSH-clamped family.

Once an LR no-SSH/radiative member is satisfactory, upscale that architecture
to HR and then add the full momentum-budget diagnostics.

## Sponge tuning update

Two one-change sponge tests were run:

- `runBC_LR_noSSH_sponge_gentle`
- `runBC_LR_noSSH_sponge_strong`

The stronger sponge slightly improved west mode 1 near the boundary:

| run | W mode-1 offset 22 | incoming fraction |
|---|---:|---:|
| gentle | 7.05 / 7.64 MW | 0.520 |
| baseline | 7.11 / 7.75 MW | 0.522 |
| strong | 7.26 / 7.63 MW | 0.512 |

At offset 33, all cases still have inward-dominated west mode 1:

| run | W mode-1 offset 33 | incoming fraction |
|---|---:|---:|
| gentle | 5.37 / 10.10 MW | 0.653 |
| baseline | 5.51 / 10.51 MW | 0.656 |
| strong | 5.70 / 10.55 MW | 0.649 |

ETAN RMS stayed near `0.009 m` in all cases.  Conclusion: sponge strength alone
is not the dominant control.  The strong sponge is modestly better near the
west boundary, but the next tests should examine forcing/radiation architecture
within the no-SSH family.

## No-west forcing update

`runBC_LR_noSSH_noW_strongSponge` removed the prescribed west normal velocity
while keeping west open/radiative/sponge and rebalanced transport over S/N.

This substantially improved the west mode-1 diagnostic:

| run | W mode-1 offset 22 | W mode-1 offset 33 |
|---|---:|---:|
| noW + strong sponge | 4.00 / 3.98 MW | 3.48 / 5.45 MW |
| full W/S/N + strong sponge | 7.26 / 7.63 MW | 5.70 / 10.55 MW |
| baseline | 7.11 / 7.75 MW | 5.51 / 10.51 MW |

But it also reduced south mode-1 outward energy at offset 22 from about
`20.6 MW` to `10.4 MW` and increased final ETAN RMS from `0.0090 m` to
`0.0115 m` with a mean bias around `-0.0078 m`.

Conclusion: west prescribed velocity is a major contributor to the west
inward-energy issue, but the simple noW rebalance is not yet a clean final
production choice.  The next logical tests are partial-west amplitude or a
south-dominant/no-west family with a more careful transport correction.
