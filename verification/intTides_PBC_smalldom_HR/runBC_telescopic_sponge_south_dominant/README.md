# South-dominant Kelvin-wave boundary test

This run tests the hypothesis that the dominant incoming barotropic tide enters
through the southern boundary, consistent with a Northern Hemisphere coastal
Kelvin wave propagating northward with the closed/eastern coast on its right.

It is cloned from `../runBC_telescopic_sponge_test`: same telescopic grid,
bathymetry, hydrography, Stevens settings, 60-cell momentum sponge, run length,
and diagnostics.

## Forcing design

The full southern M2 normal-velocity pattern from the corrected telescopic
forcing is retained as the primary imposed tide:

- `OBSvam_M2_southDominant_balanced_telescopic60.bin`
- `OBSvph_M2_southDominant_balanced_telescopic60.bin`

The detailed TPXO-derived west and north velocity patterns are omitted. Instead,
west and north carry only a spatially uniform complex barotropic return
correction that closes the southern transport:

- west `u` amplitude: 0.008148 m/s, uniform along wet boundary points
- north `v` amplitude: 0.008148 m/s on wet points, zero on dry points

The southern outward transport amplitude before balancing is 2.575e6 m3/s. The
uniform west+north return correction leaves residual complex transport of about
0.23 m3/s, effectively zero compared with the imposed tide.

This is therefore "south-dominant", not strictly "south-only". A truly
south-only velocity forcing is strongly volume-unbalanced and would excite a
large artificial domain-mean ETAN response, as seen in the earlier boundary
isolation tests.

## Sea-surface elevation

No `OB*etaFile` is prescribed in this first test. MITgcm can read OBCS eta
files, but in this configuration that is not a clean one-line addition:

- the current run uses `nonlinFreeSurf=0`; OBCS eta files are rejected unless
  nonlinear free surface is enabled;
- OBCS eta files use the generic `externForcingPeriod/externForcingCycle`
  machinery, which is currently set to 21600/43200 s for the two-record
  constant T/S boundary files, not to the M2 period.

So adding eta forcing would require a separate, explicitly documented
free-surface/forcing-frequency experiment. This run first tests whether the
south-dominant velocity forcing, with closed barotropic transport and radiative
west/north behavior, improves the internal-wave boundary diagnostics.

## Run and analyze

```bash
mpiexec -n 16 ./mitgcmuv
python check_boundary_reflections.py
python modal_characteristic_decomposition.py
```

Key comparisons:

- domain-mean ETAN M2 amplitude should remain small, unlike the unbalanced
  one-boundary isolation tests;
- west modal characteristic energy should tell us whether removing the detailed
  west TPXO pattern reduces the full-run west inward excess;
- south should remain the primary barotropic source, but south baroclinic flux
  may still be outward because locally generated internal tides can radiate
  back through the incoming barotropic boundary.
