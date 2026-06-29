# Low-resolution boundary test: full W/S/N TPXO-like forcing

This is one member of a low-resolution pair for testing boundary-condition
architecture before spending HR production time.

The HR 600 x 600 physical domain is coarsened by 3x to 200 x 200 and embedded
in a 220 x 240 telescopic grid with 20 stretched cells west, south, and north.
The eastern side remains closed. The vertical grid and stratification are kept
the same as HR. Bathymetry is block-averaged from the finite HR bathymetry and
then normally extended into the LR buffers.

This member keeps coarsened W/S/N M2 velocity patterns and applies a small
global barotropic correction to close complex transport. It is the LR analogue
of the corrected full-boundary telescopic-sponge test.

## Cost

The LR grid has 1/9 as many horizontal columns as HR and uses `deltaT=27 s`
instead of 9 s. For the same eight-M2 duration, this is about 27x cheaper in
time-step-column count. The executable is `../build_lr/mitgcmuv`.

## Configuration

- Grid: 220 x 240 x 60
- Physical domain: 200 x 200, coarsened by 3x from HR
- Buffer: 20 cells west/south/north
- Sponge: momentum-only, 20 cells, same nominal relaxation times as HR
- Diagnostics: same boundary-wave fields, four samples per M2 period
- Analysis offsets: 22, 23, 27, 33 cells from the numerical boundary, roughly
  matching HR old-domain distances 5, 10, 20, and 40 cells inward

The boundary transport balance was computed from meter-scaled spherical metrics
using the LR input grid and bathymetry. If this pair becomes the basis for a
final LR protocol, a zero-step initialization can be used to recompute the same
forcing with exact MITgcm `DXG/DYG/hFac` metrics, but the current files are
already suitable for boundary-architecture screening.

## Run and analyze

```bash
mpiexec -n 16 ./mitgcmuv
python check_boundary_reflections.py
python modal_characteristic_decomposition.py
```

Compare against `../runBC_LR_south_dominant`. The key questions are whether the
LR full-boundary case reproduces the HR symptoms: small domain-mean ETAN but
west inward modal characteristic excess.

## Completed result

The run completed normally. The post-ramp domain-mean ETAN M2 amplitude is
0.000570 m, so this LR full-boundary forcing is barotropically well behaved.

The LR case reproduces the qualitative HR west-boundary issue. Near the inner
edge of the old physical domain, west modes 1--5 are close to balanced:
18.47 MW outward and 17.22 MW inward at offset 22. Farther inward, the west
side becomes inward-dominated: 15.51 MW outward and 20.89 MW inward at offset
33. The bulk pressure-work diagnostic likewise changes from weakly positive
near the buffer to negative farther into the old domain.

South remains strongly outward in baroclinic energy, with modes 1--5 giving
about 33--37 MW outward and 9 MW inward across the sampled sections. North is
weaker and moderately reflective.
