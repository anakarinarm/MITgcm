# Boundary-condition validation: weak 30-cell sponge

This is the second boundary-condition experiment. It preserves the completed
no-sponge baseline in `../runBC_test` and uses the same grid, zero-net M2
forcing, three-period tidal ramp, Stevens boundaries, eight-period duration,
and reflection diagnostics. The only experimental change is the sponge.

## Baseline retained: no sponge

The first test completed successfully with domain-mean M2 ETAN amplitude
0.00194 m and positive outward pressure-work flux at all sampled sections.
However, its approximate reflection-amplitude proxy was high:

- West: 0.69--0.95
- South: 0.50--0.71
- North: 0.55--0.77

Those files and results remain unchanged in `../runBC_test`.

## Sponge tested here

The sponge covers 30 cells, approximately 5--6 km. Relaxation is smoothly
tapered from 2 M2 periods (89424 s) at the inner edge to M2/4 (11178 s) at
the open boundary. Both velocity components are relaxed toward the evolving
boundary values. `useLinearSponge` is false, so the relaxation target also
blends smoothly into the interior. Temperature and salinity sponge switches
are deliberately off; Stevens continues to handle tracers at the open faces.

This checkout contained an asymmetric south-salinity sponge expression. The
experiment-local `../code/obcs_sponge.F` corrects its missing
`spongeThickness-jsl` factor so all boundaries use the documented taper.

## Failed first sponge run: tracer sponge produced NaNs

The first attempt enabled sponge relaxation for U, V, temperature, and
salinity. It is retained here as part of the experiment history because it
failed decisively. The pressure-solver RHS began growing near model step 100,
reached approximately 1e15 at step 118 and 1e78 at step 119, and became NaN
at step 120 (1080 s). By 990 s, maximum U had already reached about 15 m/s.

Two isolated 150-step tests identified the trigger:

- Momentum-only sponge completed normally. At 1350 s, ETAN remained within
  approximately +/-2.1e-5 m and velocities were below 1.5e-5 m/s during the
  still-weak tidal ramp.
- Tracer-only sponge reproduced the rapid growth, overflow, and NaNs.

The boundary tracer files were also checked directly. They are finite
big-endian float64, their two time records are identical, and every boundary
value exactly matches the horizontally uniform initial T/S profile. Thus the
failure is not caused by malformed files or an initial boundary mismatch; it
is associated with the OBCS tracer-sponge code path in this configuration.
For this internal-tide test, the robust and physically appropriate choice is
therefore a momentum-only sponge. Do not interpret or compare the output from
the failed tracer-sponge run; remove or archive it before rerunning.

## Run and analyze

Run from this directory with the rebuilt executable:

```bash
mpiexec -n 16 ../build/mitgcmuv
python check_boundary_reflections.py
```

Changing the runtime sponge switches does not require recompilation. The
existing executable already includes `ALLOW_OBCS_SPONGE`.

The run writes about 8.3 GB of single-precision U, V, PHIHYD, and ETAN output.
Compare its `boundary_reflection_summary.csv` directly with the baseline file
in `../runBC_test` before changing the production configuration.

## Completed result

The weak 30-cell sponge did not materially reduce reflection. At the
40-cell section, the first sampled section outside the sponge, the approximate
reflection-amplitude proxy changed as follows:

- West: 0.953 without sponge -> 0.955 with sponge
- South: 0.709 -> 0.703
- North: 0.770 -> 0.755

The changes are negligible relative to the limitations of the proxy. Values at
offsets 2--20 lie inside the damped region and should not be interpreted as
conservative incoming/outgoing-wave decompositions.

The likely reason is scale separation. A WKB estimate from the initial
stratification gives M2 wavelengths of approximately 174, 87, 58, 43, and
35 km for modes 1--5. The 30-cell sponge is only 5--6 km wide, and its original
time scale is too slow for low modes to attenuate appreciably while crossing
it.

The next useful test should therefore not be a small adjustment to this sponge.
Use a 60--80-cell buffer outside the scientific analysis region and a stronger,
smooth Rayleigh taper. If the existing domain cannot be extended, a practical
native-OBCS test is 60 cells with a zero-at-inner-edge taper and an outer
relaxation time around M2/16 (2794.5 s), followed by modal—not bulk—reflection
diagnostics outside the sponge. The preferable long-term solution is to extend
the domain so the absorbing layer does not consume the region of interest.
