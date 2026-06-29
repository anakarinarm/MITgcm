# Boundary-isolation test: south+north forcing, no west tide

This run is cloned from `../runBC_telescopic_sponge_test` after the corrected
telescopic forcing and 60-cell sponge were established.

South and north M2 normal velocities are forced, but the west M2 normal
velocity is zero:

- west `u` tide: zero-amplitude file
- south `v` tide: `../input/OBSvam_M2_zeroNet_telescopic60.bin`
- north `v` tide: `../input/OBNvam_M2_zeroNet_telescopic60.bin`

The grid, bathymetry, hydrography, Stevens boundary settings, sponge settings,
run length, and diagnostics are otherwise identical to the corrected
telescopic-sponge control.

Run and analyze:

```bash
mpiexec -n 16 ./mitgcmuv
python check_boundary_reflections.py
python modal_characteristic_decomposition.py
```

Interpretation: this is the cleanest test for reflection/radiation at the west
boundary when the west tide is not directly prescribed. If modal characteristic
energy near the west side is outward-dominated here, then the full-run inward
west signal is mainly the prescribed west tide. If it is still inward-dominated,
the west boundary/sponge/corner configuration is likely creating or trapping an
incoming internal-wave signal even without west tidal forcing.

## Completed result and caveat

This intentionally isolated run is strongly volume-unbalanced by itself. The
post-ramp domain-mean ETAN M2 amplitude is 0.762594 m, nearly opposite in phase
to the west-only case. The isolated south+north transport amplitude is about
1.746e6 m3/s; in the full corrected run this is nearly cancelled by the west
transport. Therefore this run is useful as evidence that the boundary groups
are mutually balancing the barotropic tide, but it is not a clean internal-wave
reflection test.

The modal diagnostic shows near-balance at the west 5-cell old-domain section:
summed modes 1--5 give 23.92 MW outward and 23.80 MW inward. At the 40-cell
section the west side is inward-dominated, with 21.89 MW outward and 26.08 MW
inward. Because the large barotropic basin response contaminates the
interpretation, use `../runBC_telescopic_sponge_noW_balanced` for the cleaner
follow-up.
