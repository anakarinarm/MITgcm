# HR full-TPXO GEBCO W117.75 production baseline

This is the production version of the validated HR pilot:

```text
../runBC_HR_fullTPXO_obeta_gebcoW11775_pilot
```

The pilot was accepted because the ETAN/phase and W/T/speed animations looked
clean, the western generation site was physically consistent and far enough
from the open boundary, and the boundary/modal diagnostics did not show a large
configuration-induced artifact.

## Geometry and forcing

- Grid: `750 x 720 x 60`.
- Physical domain: `690 x 600`, plus a `60`-cell W/S/N sponge.
- Physical west edge: `242.25 E` = `117.75 W`.
- Physical south edge: `31.25 N`.
- East boundary: closed land wall.
- Bathymetry: GEBCO, with W/S/N sponge bathymetry flattened normal to the
  boundaries.
- Boundary forcing: full TPXO M2 SSH + barotropic velocity records on W/S/N.
- OBCS: prescribed SSH/velocity, nonlinear free surface/r*, momentum sponge,
  no Stevens for this prescribed-SSH configuration.

This run uses the same input binaries and dedicated build as the HR pilot:

```text
../code_hr_obeta_gebcoW11775
../build_hr_obeta_gebcoW11775
```

## Diagnostics

Production diagnostics are written to `output/`. If possible, replace this
directory with a symlink to external storage before running, e.g.

```bash
mkdir -p /Volumes/LaCie/MITgcm_outputs/intTides_PBC_smalldom_HR/runProd_HR_fullTPXO_gebcoW11775_momentumEnergy
rm -rf output
ln -s /Volumes/LaCie/MITgcm_outputs/intTides_PBC_smalldom_HR/runProd_HR_fullTPXO_gebcoW11775_momentumEnergy output
```

Diagnostic streams:

- `output/diag_state_3d`: phase-resolving 3-D state, pressure, and velocity-flux
  terms for energy-flux/conversion work.
- `output/diag_state_2d`: ETAN and PHIBOT.
- `output/diag_mom_budget`: averaged direct momentum-budget terms.
- `output/diag_mom_vecinv`: averaged vector-invariant decomposition terms.

The output begins after the 3-M2 forcing ramp:

```text
timePhase = 134136 s = 3 M2
frequency = M2/8 = 5589 s
endTime   = 357696 s = 8 M2
```

So this first production baseline saves five M2 periods of analysis output.

## Build/run

If the pilot binary was rebuilt after adding
`code_hr_obeta_gebcoW11775/DIAGNOSTICS_SIZE.h`, no new rebuild is needed.
Otherwise rebuild:

```bash
cd intTides_PBC_smalldom_HR/build_hr_obeta_gebcoW11775
make clean
make depend
make -j
```

Run:

```bash
cd ../runProd_HR_fullTPXO_gebcoW11775_momentumEnergy
mpiexec -n 20 ./mitgcmuv
```

## Post-processing

Animations:

```bash
python3 animate_etan_phase.py
python3 animate_depth_slices.py
python3 animate_depth_slices_with_sponge.py
```

Budget/energy scripts copied from the earlier HR production workflow:

```bash
python3 compute_barotropic_to_baroclinic_conversion.py
python3 plot_m2_barotropic_tide_diagnostics.py
python3 plot_depth_integrated_momentum_phase_maps.py
```

As usual: if one of these scripts needs a small path/prefix update after the
first run, update the script rather than changing the model configuration.
