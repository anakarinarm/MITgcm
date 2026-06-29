# LR full-TPXO prescribed SSH with volume-balanced velocity records

This is the follow-up to `../runBC_LR_fullTPXO_obeta`.

The raw SSH+velocity test showed that prescribing SSH helped the west
modal-characteristic split, but raw full-boundary TPXO velocities created a
large barotropic transport imbalance and very large absolute boundary energies.
This run keeps the same prescribed TPXO SSH records, keeps tangential velocity
records unchanged, and only replaces the normal velocity records with
volume-balanced versions.

## What changed relative to `../runBC_LR_fullTPXO_obeta`

- Same nonlinear-free-surface executable: `../build_lr_obeta/mitgcmuv`.
- Same minimum-wet-depth bathymetry:
  `../input/smalldom_LR220x240_telescopic20_bat_min10m_obeta.bin`.
- Same SSH files:
  - `OBWeta_M2_fullTPXO_LR220x240_telescopic20_obeta16.bin`
  - `OBSeta_M2_fullTPXO_LR220x240_telescopic20_obeta16.bin`
  - `OBNeta_M2_fullTPXO_LR220x240_telescopic20_obeta16.bin`
- Same tangential velocity files:
  - `OBWv_M2_fullTPXO_LR220x240_telescopic20_obeta16.bin`
  - `OBSu_M2_fullTPXO_LR220x240_telescopic20_obeta16.bin`
  - `OBNu_M2_fullTPXO_LR220x240_telescopic20_obeta16.bin`
- New balanced normal velocity files:
  - `OBWu_M2_fullTPXO_balanced_LR220x240_telescopic20_obeta16.bin`
  - `OBSv_M2_fullTPXO_balanced_LR220x240_telescopic20_obeta16.bin`
  - `OBNv_M2_fullTPXO_balanced_LR220x240_telescopic20_obeta16.bin`

`useOBCSbalance` remains `.FALSE.` because the correction is already baked into
the files.

## Balancing method

`prepare_balanced_obeta_velocity.py` reads the raw time-record velocities and
the initialized MITgcm grid metrics from `../runBC_LR_fullTPXO_obeta`.

For each M2 record, it computes

```text
Q_in = Q_W + Q_S - Q_N
```

where west/south positive velocities are inflow and north positive velocity is
outflow, matching MITgcm/OBCS convention.  It then applies one depth-uniform
velocity correction over all wet W/S/N normal boundary faces:

```text
OBWu <- OBWu - Q_in/A_total
OBSv <- OBSv - Q_in/A_total
OBNv <- OBNv + Q_in/A_total
```

This preserves the SSH forcing and most of the TPXO velocity structure while
removing the basin-filling barotropic transport defect.

Generated correction diagnostics:

```text
area W/S/N/total [m2]: 2.015649e8 / 1.396747e8 / 6.120911e7 / 4.024487e8
raw Q_W + Q_S - Q_N max abs: 8.3400106e5 m3/s
balanced residual max abs: 6.98e-10 m3/s
correction velocity max abs: 0.002072 m/s
correction velocity rms: 0.001494 m/s
```

## Run

```bash
python3 prepare_balanced_obeta_velocity.py
mpiexec -n 16 ./mitgcmuv
python3 check_boundary_reflections.py
python3 modal_characteristic_decomposition.py
```

Compare primarily against:

- `../runBC_LR_fullTPXO_obeta`
- `../runBC_LR_fullTPXO_balanced`
- `../runBC_LR_south_dominant`

The success criterion is not only a clean west characteristic split.  We also
need the absolute characteristic energies and north/south bulk flux metrics to
fall back near the earlier balanced LR controls.

## Result

The run completed cleanly.  There was no `CALC_R_STAR` failure and no model
`STOP`; the only `STDERR.*` entries were the old retired `EXACT_CONSERV`
warnings.

Balancing the normal velocity records helped substantially compared with the
raw SSH+velocity test, but it did not make this configuration production-ready.

At 22-cell offset:

```text
case                         W bulk flux    S bulk flux    N bulk flux
LR fullTPXO balanced          +0.8 MW       +19.7 MW        +0.6 MW
LR raw SSH+velocity         +162.6 MW      +473.7 MW       -85.6 MW
LR balanced SSH+velocity     +15.3 MW       -20.0 MW        -7.3 MW
```

For west mode 1:

```text
case                         inward/outward amp   incoming frac   out/in energy
LR fullTPXO balanced              1.044              0.522          7.1/7.8 MW
LR raw SSH+velocity               0.990              0.495       9278/9101 MW
LR balanced SSH+velocity          0.984              0.492        204/198 MW
```

Interpretation: volume-balancing the prescribed velocity records removes much
of the pathological raw-obeta amplification and preserves the improved west
characteristic split.  But the absolute modal energies are still one to two
orders of magnitude too large, and south/north bulk fluxes have the wrong
inward tendency.  So this test is useful diagnostically, but I would not port it
to HR as the production boundary condition.

The likely issue is over-constraint: prescribing SSH and the full barotropic
velocity pattern at every open boundary is still too stiff for this cropped
domain, even after net transport is closed.  The next test should use the SSH
constraint more gently, e.g. prescribed SSH with velocity radiation/relaxation
instead of fully prescribed velocity records.
