# Boundary-isolation test: west-only forcing

This run is cloned from `../runBC_telescopic_sponge_test` after the corrected
telescopic forcing and 60-cell sponge were established.

Only the west M2 normal velocity is forced:

- west `u` tide: `../input/OBWuam_M2_zeroNet_telescopic60.bin`
- south `v` tide: zero-amplitude file
- north `v` tide: zero-amplitude file

The grid, bathymetry, hydrography, Stevens boundary settings, sponge settings,
run length, and diagnostics are otherwise identical to the corrected
telescopic-sponge control.

Run and analyze:

```bash
mpiexec -n 16 ./mitgcmuv
python check_boundary_reflections.py
python modal_characteristic_decomposition.py
```

Interpretation: if this west-only run still shows inward-dominated modal
characteristic energy at the west side, then the negative west pressure-work
signal in the full run is mostly prescribed tide injection from the west. If
the west side becomes outward-dominated here but not in the full run, then
south/north forcing or corner interactions are feeding the west-boundary issue.

## Completed result and caveat

This intentionally isolated run is strongly volume-unbalanced by itself. The
post-ramp domain-mean ETAN M2 amplitude is 0.763448 m. The isolated west
transport amplitude is about 1.746e6 m3/s; in the full corrected run this is
nearly cancelled by the south+north transport. Therefore this run is useful as
evidence that the boundary groups are mutually balancing the barotropic tide,
but it is not a clean internal-wave reflection test.

The modal diagnostic nevertheless shows that the west section is
outward-dominated in this unbalanced west-only run: summed modes 1--5 at the
5-cell old-domain section give 3.94 MW outward and 2.69 MW inward. Because the
large barotropic basin response contaminates the interpretation, use
`../runBC_telescopic_sponge_Wonly_balanced` for the cleaner follow-up.
