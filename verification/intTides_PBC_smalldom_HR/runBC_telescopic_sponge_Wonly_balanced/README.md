# Balanced boundary-isolation test: west-pattern forcing

This is the cleaner follow-up to `../runBC_telescopic_sponge_Wonly`. The
unbalanced west-only run produced a post-ramp domain-mean ETAN M2 amplitude of
0.763448 m because the isolated west transport amplitude is about 1.746e6 m3/s.

Here the west M2 normal-velocity pattern is retained, but a uniform complex
barotropic correction over the wet W/S/N open-boundary faces restores zero net
transport. The correction amplitude is 0.003667 m/s and the residual complex
transport is about 0.046 m3/s.

The grid, bathymetry, hydrography, Stevens settings, 60-cell sponge, run length,
and diagnostics are identical to the corrected telescopic-sponge control.

Run and analyze:

```bash
mpiexec -n 16 ./mitgcmuv
python check_boundary_reflections.py
python modal_characteristic_decomposition.py
```

Interpretation: compare this run with `../runBC_telescopic_sponge_noW_balanced`
and `../runBC_telescopic_sponge_test`. If the west modal characteristic energy
is outward-dominated here, the full-run west inward excess is not simply caused
by the west forcing pattern alone.
