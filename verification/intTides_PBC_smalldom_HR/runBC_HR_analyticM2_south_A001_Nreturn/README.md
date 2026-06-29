# HR analytical M2 forcing: south A0=0.01 m/s, north-only return

This is the first HR boundary-test version of the preferred LR analytical
forcing architecture.

It is copied from `runBC_telescopic_sponge_test`, but replaces the TPXO-derived
tidal velocity files with controlled analytical M2 forcing:

- no prescribed SSH/elevation;
- south normal M2 velocity with `A0 = 0.01 m/s`;
- 30-cell cosine taper near south corners, physically comparable to the 10-cell
  LR taper;
- north uniform compensating return flow;
- west prescribed velocity exactly zero;
- Stevens/radiative W/S/N;
- 60-cell strong momentum sponge, matching the preferred LR analytical family;
- boundary/reflection diagnostics only, no momentum-budget diagnostics yet.

The purpose is to verify the HR boundary behavior before promoting the
configuration to a production momentum-budget run.

## Generated forcing

Created by:

```bash
python3 prepare_analytic_south_tide_HR.py
```

Output files:

- `../input/OBWuam_M2_analyticSouth_A001_Nreturn_HR660x720_telescopic60.bin`
- `../input/OBWuph_M2_analyticSouth_A001_Nreturn_HR660x720_telescopic60.bin`
- `../input/OBSvam_M2_analyticSouth_A001_Nreturn_HR660x720_telescopic60.bin`
- `../input/OBSvph_M2_analyticSouth_A001_Nreturn_HR660x720_telescopic60.bin`
- `../input/OBNvam_M2_analyticSouth_A001_Nreturn_HR660x720_telescopic60.bin`
- `../input/OBNvph_M2_analyticSouth_A001_Nreturn_HR660x720_telescopic60.bin`

Generation diagnostics:

- south max/RMS amplitude: `0.0100 / 0.00970 m/s`
- north return max/RMS amplitude: `0.01848 / 0.01390 m/s`
- west max/RMS amplitude: `0.0 / 0.0 m/s`
- transport residual `Q_W + Q_S - Q_N`: `2.3e-10 m3/s`

## Run and analyze

```bash
mpiexec -n 16 ./mitgcmuv
python3 check_boundary_reflections.py
python3 modal_characteristic_decomposition.py
```

Compare primarily against the LR analytical baseline:

- `../runBC_LR_analyticM2_south_A001_Nreturn`

Decision criteria:

- ETAN should remain small and unbiased.
- South should produce a clear outgoing internal-tide signal.
- West should be much less contaminated than the old TPXO/SSH setups.
- North return should not dominate the modal diagnostics.
