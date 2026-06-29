# Boundary-condition control: telescopic grid without sponge

This is the no-sponge member of a controlled pair. It is identical to
`../runBC_telescopic_sponge_test` in grid, bathymetry, hydrography, open-boundary
formulation, corrected tidal forcing, run length, and diagnostics. Its only
experimental difference is `useOBCSsponge=.FALSE.` and `spongeThickness=0`.

The complete old 600 x 600 physical domain remains unchanged inside a
660 x 720 grid. Sixty geometrically stretched cells form numerical buffers at
west, south, and north, with a maximum adjacent-size ratio of 1.025. Bathymetry
is extruded normally into these buffers, so no new topographic feature exists.

The M2 harmonics retain the old forcing over the original coordinate spans.
Over added corner segments they taper smoothly to zero with a raised cosine,
preventing duplicated full-strength forcing around the southwest and northwest
corners. The subsequent balancing correction leaves only about 0.000376 m/s
at the outer corner. The uniform wet-face correction is 0.000375877 m/s. It was
computed using this initialized grid's exact `hFacW`, `hFacS`, `DXG`, `DYG`,
and `DRF`; residual complex transport is 2.4e-10 m3/s.

Run and analyze:

```bash
mpiexec -n 16 ./mitgcmuv
python check_boundary_reflections.py
```

The checker samples old-domain locations 5, 10, 20, and 40 cells inward from
the buffer and uses cumulative physical distance on the stretched grid. Compare
its CSV directly with the sponge-pair result. This run isolates effects of the
grid transition and corrected boundary forcing from sponge damping.

## Completed test result

The eight-M2 run completed normally. The corrected boundary forcing removed
the large artificial basin-filling signal: the post-ramp domain-mean ETAN M2
amplitude is 0.000873 m, compared with 0.008976 m in the invalid first
telescopic setup.

The pressure-work diagnostic still shows inward-dominated baroclinic flux at
the west side of the original physical domain, with outward flux from -5.6 to
-9.8 MW over the 5--40-cell sections and negative progressive-wave index. The
south boundary is outward and more progressive than the old-grid sponge case
(reflection-amplitude proxy about 0.49--0.53), while the north boundary remains
moderately reflective (about 0.73--0.74). Because the west metric is undefined
when flux is inward, the next boundary-design work should focus on west-side
forcing/radiation rather than only sponge strength.

## Modal characteristic diagnostic

`modal_characteristic_decomposition.py` was added here for consistency with the
sponge member. It uses the initial T/S stratification and a WKB/cosine
hydrostatic-mode approximation. For the no-sponge telescopic run, summing modes
1--5 gives west outward/inward/net characteristic energy of 17.1/23.4/-6.3 MW
at the 5-cell old-domain section and 14.5/26.5/-12.0 MW at the 40-cell section.
The 60-cell sponge therefore reduces the west inward excess slightly, but does
not change the qualitative diagnosis.
