# LR full-TPXO SSH + weak balanced velocity test, alpha = 0.10

Purpose: a gentler companion to the `alpha = 0.25` experiment.  This tests
whether TPXO SSH alone supplies most of the useful barotropic phase information
while only a very weak velocity constraint is needed to keep the boundary from
drifting.

This run keeps the same LR telescopic grid, minimum-10 m obeta bathymetry,
prescribed TPXO M2 SSH, T/S boundary files, diagnostics, and 20-cell momentum
sponge used by `runBC_LR_fullTPXO_obeta_balanced`.

Difference from `runBC_LR_fullTPXO_obeta_balanced`:

- all prescribed boundary velocity records are scaled by `alpha = 0.10`;
- normal velocity records start from the volume-balanced obeta fields, so the
  net prescribed barotropic transport remains zero at every forcing record;
- tangential velocity records are also scaled by 0.10;
- `useOBCSbalance=.FALSE.` because balancing has already been applied in the
  binary files.

Generated velocity files:

- `../input/OBWu_M2_fullTPXO_balanced_weakvel010_LR220x240_telescopic20_obeta16.bin`
- `../input/OBWv_M2_fullTPXO_weakvel010_LR220x240_telescopic20_obeta16.bin`
- `../input/OBSu_M2_fullTPXO_weakvel010_LR220x240_telescopic20_obeta16.bin`
- `../input/OBSv_M2_fullTPXO_balanced_weakvel010_LR220x240_telescopic20_obeta16.bin`
- `../input/OBNu_M2_fullTPXO_weakvel010_LR220x240_telescopic20_obeta16.bin`
- `../input/OBNv_M2_fullTPXO_balanced_weakvel010_LR220x240_telescopic20_obeta16.bin`

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

Interpretation target: if this run is much calmer than `alpha = 0.25` while
still keeping ETAN controlled, it is probably the better SSH-constrained
boundary architecture to upscale to HR.

## Result after run

The run completed cleanly; STDERR only contained the harmless retired
`EXACT_CONSERV` warning.

At offset 22 cells, bulk outward fluxes from `check_boundary_reflections.py`
were:

- W: `+0.153 MW`
- S: `-0.200 MW`
- N: `-0.070 MW`

This is much quieter than both the full-velocity SSH run and the `alpha = 0.25`
run, again close to the expected `alpha^2` scaling of flux.

The modal characteristic decomposition still showed almost exactly equal
outward and inward mode-1 energies:

- W mode 1: `182.7 / 182.6 MW`, incoming fraction `0.500`
- S mode 1: `85.0 / 85.2 MW`, incoming fraction `0.500`
- N mode 1: `53.2 / 53.3 MW`, incoming fraction `0.500`

Final saved ETAN range was approximately `-0.425 to +0.0030 m`.  The large
negative value matches the prescribed TPXO SSH phase/amplitude, not a runaway
transport imbalance.

Interpretation: this is the calmest SSH+velocity test in terms of net flux, but
the full-amplitude prescribed SSH still dominates the boundary-characteristic
diagnostic and looks like a nearly standing imposed signal.  It is not obvious
that this boundary architecture is safer than the no-SSH balanced/radiative
setup for production.
