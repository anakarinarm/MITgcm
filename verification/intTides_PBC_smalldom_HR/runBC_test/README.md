# Boundary-condition validation run

This is an eight-M2-period test of the Stevens open boundaries and the
three-period smooth tidal ramp. It intentionally omits the momentum-budget
diagnostics and permanent checkpoints.

The field output is `bc_wave_3d` (U, V and PHIHYD at all 60 levels) plus
`bc_wave_eta` in single precision at four phases per M2 period. Run
`check_boundary_reflections.py` afterward. It excludes the first three ramp
periods and evaluates baroclinic pressure-work flux on sections parallel to
each open boundary.

Run the analysis with:

```bash
python3 check_boundary_reflections.py
```

It writes `boundary_reflection_summary.csv` and
`boundary_reflection_summary.png`. Positive flux is outward. A negative flux
is inward-dominated. The reflection-amplitude ratio is an equal-impedance,
single-mode proxy; judge it by its consistency across several offsets rather
than by one section alone.

For the completed no-sponge test, the domain-mean M2 ETAN amplitude is
0.00194 m. Outward pressure-work flux is positive at all sampled sections,
but the reflection-amplitude proxy is large: approximately 0.69--0.95 west,
0.50--0.71 south, and 0.55--0.77 north. These values are not exact multimode
reflection coefficients, but their consistency indicates that the no-sponge
configuration is not sufficiently absorbing.

This no-sponge run is retained as the baseline. The follow-up configuration
is in `../runBC_sponge_test`, using a documented weak 30-cell sponge; it does
not overwrite any files or conclusions from this test.

Run with the rebuilt 16-rank executable from `../build/mitgcmuv`.

## Tidal transport preparation

`prepare_tidal_bc.py` interpolates TPXO velocity as a complex harmonic onto
the actual MITgcm U/V boundary faces. It then applies the minimum area-weighted
uniform correction shared by all wet open-boundary faces so that the complex
net M2 transport is zero. The generated `*_M2_zeroNet.bin` files are used by
this run.

The correction amplitude is 0.00183 m/s. The residual complex transport is
approximately 1.3e-10 m3/s. This is an intentionally idealized current-forced
experiment: the original TPXO velocities implied about 0.5 m of M2 elevation,
consistent with TPXO's local elevation solution, but that domain-mean tide has
been removed here by design.
