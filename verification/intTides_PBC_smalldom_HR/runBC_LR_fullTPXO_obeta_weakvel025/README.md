# LR full-TPXO SSH + weak balanced velocity test, alpha = 0.25

Purpose: test whether prescribed TPXO sea-surface height can provide the useful
barotropic phase constraint without over-constraining the model with full TPXO
boundary velocities.

This run keeps the same LR telescopic grid, minimum-10 m obeta bathymetry,
prescribed TPXO M2 SSH, T/S boundary files, diagnostics, and 20-cell momentum
sponge used by `runBC_LR_fullTPXO_obeta_balanced`.

Difference from `runBC_LR_fullTPXO_obeta_balanced`:

- all prescribed boundary velocity records are scaled by `alpha = 0.25`;
- normal velocity records start from the volume-balanced obeta fields, so the
  net prescribed barotropic transport remains zero at every forcing record;
- tangential velocity records are also scaled by 0.25 so the boundary is not
  forced with a full-amplitude shear while the normal flow is weakened;
- `useOBCSbalance=.FALSE.` because balancing has already been applied in the
  binary files.

Generated velocity files:

- `../input/OBWu_M2_fullTPXO_balanced_weakvel025_LR220x240_telescopic20_obeta16.bin`
- `../input/OBWv_M2_fullTPXO_weakvel025_LR220x240_telescopic20_obeta16.bin`
- `../input/OBSu_M2_fullTPXO_weakvel025_LR220x240_telescopic20_obeta16.bin`
- `../input/OBSv_M2_fullTPXO_balanced_weakvel025_LR220x240_telescopic20_obeta16.bin`
- `../input/OBNu_M2_fullTPXO_weakvel025_LR220x240_telescopic20_obeta16.bin`
- `../input/OBNv_M2_fullTPXO_balanced_weakvel025_LR220x240_telescopic20_obeta16.bin`

They can be regenerated with:

```sh
python3 prepare_weak_obeta_velocity.py
```

Run normally from this directory with the existing `mitgcmuv` symlink.

After completion, run:

```sh
python3 check_boundary_reflections.py
python3 modal_characteristic_decomposition.py
```

Interpretation target: compared with the full-amplitude balanced obeta run, we
want much smaller absolute modal characteristic energies while retaining a
reasonable outgoing/incoming split and avoiding large ETAN growth.

## Result after run

The run completed cleanly; STDERR only contained the harmless retired
`EXACT_CONSERV` warning.

At offset 22 cells, bulk outward fluxes from `check_boundary_reflections.py`
were:

- W: `+0.958 MW`
- S: `-1.251 MW`
- N: `-0.443 MW`

This is a large reduction relative to the full-velocity SSH run
(`+15.3, -20.0, -7.26 MW`), close to the expected `alpha^2` scaling of flux.

The modal characteristic decomposition still showed almost equal outward and
inward mode-1 energies:

- W mode 1: `183.8 / 183.4 MW`, incoming fraction `0.499`
- S mode 1: `85.2 / 86.2 MW`, incoming fraction `0.503`
- N mode 1: `53.7 / 54.0 MW`, incoming fraction `0.501`

Final saved ETAN range was approximately `-0.425 to +0.0075 m`.  The large
negative value matches the prescribed TPXO SSH phase/amplitude, not a runaway
transport imbalance.

Interpretation: reducing velocity amplitude fixed the excessive net boundary
flux, but the prescribed full-amplitude SSH still creates a near-standing
baroclinic characteristic signal at the boundaries.
