# LR prescribed SSH with Orlanski/radiative velocity boundaries

This is the gentler follow-up to:

- `../runBC_LR_fullTPXO_obeta`
- `../runBC_LR_fullTPXO_obeta_balanced`

Those tests showed that prescribing TPXO SSH improves the west
modal-characteristic split, but prescribing full barotropic velocity records on
all open boundaries over-constrains the cropped domain and inflates the absolute
boundary energies.

This test keeps the SSH constraint but does **not** prescribe U/V files.
Instead, velocities are computed by Orlanski radiation at the W/S/N open
boundaries.  T/S are still prescribed from the same 16-record constant boundary
files, so the stratification at the boundary remains controlled.

## Dedicated build

This run uses:

```text
../code_lr_obeta_orlanski
../build_lr_obeta_orlanski/mitgcmuv
```

Relative to `../code_lr_obeta`, the only intended code-option change is:

```text
#define ALLOW_ORLANSKI
```

The nonlinear free-surface options from the SSH-prescribed tests are retained.
Stevens remains disabled at runtime because Stevens + nonlinear free surface is
not implemented in this MITgcm/OBCS path.

## Boundary setup

In `data.obcs`:

- `OBWetaFile`, `OBSetaFile`, and `OBNetaFile` are prescribed.
- T/S boundary files are prescribed.
- U/V files are deliberately omitted.
- `useOrlanskiWest/South/North=.TRUE.`
- `useOBCStides=.FALSE.`
- the momentum sponge is retained.

This means `OBCS_CALC` computes radiative U/V/T/S first, then
`OBCS_PRESCRIBE_READ` overwrites only the supplied SSH and T/S fields.

## Run

```bash
mpiexec -n 16 ./mitgcmuv
python3 check_boundary_reflections.py
python3 modal_characteristic_decomposition.py
```

Compare against:

- `../runBC_LR_fullTPXO_balanced`
- `../runBC_LR_fullTPXO_obeta`
- `../runBC_LR_fullTPXO_obeta_balanced`

Success would mean:

- west modal-characteristic split remains near balanced;
- absolute modal energies drop near the earlier LR fullTPXO balanced control;
- south/north bulk fluxes no longer show the inward tendency seen in
  `runBC_LR_fullTPXO_obeta_balanced`.

## Result: blocked by MITgcm/OBCS compatibility

This run stops during `OBCS_CHECK` before time stepping:

```text
OBCS_CHECK: nonlinFreeSurf not yet implemented in Orlanski OBC
OBCS_CHECK: detected  1 fatal error(s)
```

So this boundary combination is not available in the stock MITgcm/OBCS path we
are using.  Prescribed `OB*etaFile` requires nonlinear free surface, but
Orlanski is explicitly disabled with nonlinear free surface.  Stevens is also
disabled with nonlinear free surface in this code path.  Therefore the clean
"prescribed SSH + radiative velocity" experiment cannot be run without modifying
MITgcm source code more deeply.

This is useful information: the practical choices left in the current code are
either:

1. prescribed SSH plus prescribed velocity records, tuned to be less
   over-constraining; or
2. no prescribed SSH, with balanced/radiative velocity forcing as in the earlier
   LR controls.
