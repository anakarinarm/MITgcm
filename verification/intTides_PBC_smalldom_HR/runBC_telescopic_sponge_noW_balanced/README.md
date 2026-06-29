# Balanced boundary-isolation test: south+north-pattern forcing

This is the cleaner follow-up to `../runBC_telescopic_sponge_noW`. The
unbalanced no-west run produced a post-ramp domain-mean ETAN M2 amplitude of
0.762594 m because the isolated south+north transport amplitude is about
1.746e6 m3/s.

Here the south+north M2 normal-velocity patterns are retained, but a uniform
complex barotropic correction over the wet W/S/N open-boundary faces restores
zero net transport. The correction amplitude is 0.003667 m/s and the residual
complex transport is about 0.046 m3/s.

The grid, bathymetry, hydrography, Stevens settings, 60-cell sponge, run length,
and diagnostics are identical to the corrected telescopic-sponge control.

Run and analyze:

```bash
mpiexec -n 16 ./mitgcmuv
python check_boundary_reflections.py
python modal_characteristic_decomposition.py
```

Interpretation: this is the cleaner test for west-boundary reflection when the
original west tidal pattern is absent. If west remains inward-dominated here,
then south/north forcing, corner interactions, or the west radiation/sponge
configuration are feeding the west-side inward characteristic signal.
