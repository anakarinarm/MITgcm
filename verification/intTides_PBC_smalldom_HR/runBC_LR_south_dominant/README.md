# Low-resolution boundary test: south-dominant Kelvin-wave forcing

This is the south-dominant member of the low-resolution boundary-architecture
pair. It is identical to `../runBC_LR_fullTPXO_balanced` except for the M2 tidal
velocity files.

The southern M2 normal-velocity pattern is retained as the dominant incoming
Kelvin-wave forcing. The detailed west and north TPXO-like patterns are
omitted; west and north contain only spatially uniform barotropic return flow
to close complex transport.

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

The south-only transport amplitude before balancing is about 2.43e6 m3/s using
the LR meter-scaled spherical metrics. Closing that transport through uniform
west+north return flow requires amplitude about 0.00924 m/s. This is a
controlled diagnostic forcing, not a claim that the return flow is physically
resolved TPXO structure.

## Run and analyze

```bash
mpiexec -n 16 ./mitgcmuv
python check_boundary_reflections.py
python modal_characteristic_decomposition.py
```

Compare against `../runBC_LR_fullTPXO_balanced`. If this member keeps
domain-mean ETAN small while reducing west inward modal characteristic energy,
then south-dominant forcing is a strong candidate to port back to HR.

## Completed result

The run completed normally, but the post-ramp domain-mean ETAN M2 amplitude is
0.022803 m, about 40 times larger than the LR full-boundary member
(0.000570 m). This is not the catastrophic 0.76 m response from the unbalanced
one-boundary isolation tests, but it shows that the uniform west+north return
flow is dynamically less clean than the coarsened full-boundary tide.

South-dominant forcing reduces the absolute baroclinic energy levels but does
not remove the west-boundary pattern. At offset 22, west modes 1--5 are
15.61 MW outward and 14.44 MW inward. At offset 33, they are 12.97 MW outward
and 17.45 MW inward, still inward-dominated. The bulk pressure-work diagnostic
also remains weakly positive near the buffer and negative farther into the old
physical domain.

Conclusion from the LR pair: south-dominant forcing alone is not an obvious
improvement. It slightly lowers the west energy levels, but keeps the same
west inward-excess structure and increases domain-mean ETAN substantially.
