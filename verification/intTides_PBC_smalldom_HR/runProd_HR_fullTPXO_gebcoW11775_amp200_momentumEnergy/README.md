# HR full-TPXO GEBCO W117.75 amplitude sensitivity: 2.0x tide

This run is the double-amplitude member of the preliminary tidal-amplitude
sensitivity pair.

It is identical to the accepted production baseline:

```text
../runProd_HR_fullTPXO_gebcoW11775_momentumEnergy
```

except that the prescribed TPXO M2 SSH and barotropic velocity boundary records
on W/S/N are multiplied by `2.0`.

Unchanged:

- grid, bathymetry, sponge, vertical coordinates, T/S initial and boundary
  fields;
- OBCS architecture;
- viscosity/diffusivity/numerics;
- run length and diagnostics;
- executable/build.

Changed:

- `OBWeta`, `OBSeta`, `OBNeta`;
- `OBWu`, `OBWv`, `OBSu`, `OBSv`, `OBNu`, `OBNv`.

The scaled files were generated from the baseline files by:

```bash
cd ../runProd_HR_fullTPXO_gebcoW11775_momentumEnergy
python make_tpxo_amplitude_sensitivity_inputs.py
```

## Output location

Recommended before running:

```bash
mkdir -p /Volumes/LaCie/MITgcm_outputs/intTides_PBC_smalldom_HR/runProd_HR_fullTPXO_gebcoW11775_amp200_momentumEnergy
rm -rf output
ln -s /Volumes/LaCie/MITgcm_outputs/intTides_PBC_smalldom_HR/runProd_HR_fullTPXO_gebcoW11775_amp200_momentumEnergy output
```

## Run

```bash
cd intTides_PBC_smalldom_HR/runProd_HR_fullTPXO_gebcoW11775_amp200_momentumEnergy
mpiexec -n 20 ./mitgcmuv
```

## Post-processing

```bash
python animate_etan_phase.py
python animate_depth_slices.py
python animate_depth_slices_with_sponge.py
python compute_barotropic_to_baroclinic_conversion.py
python plot_m2_barotropic_tide_diagnostics.py
python plot_depth_integrated_momentum_phase_maps.py
```
