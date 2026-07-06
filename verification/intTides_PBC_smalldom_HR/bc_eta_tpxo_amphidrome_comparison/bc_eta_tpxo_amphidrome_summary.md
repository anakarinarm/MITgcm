# Boundary-test M2 SSH amphidrome comparison to TPXO9

TPXO elevation file: `/Users/karina/Research/Canyons/data/tides/h_tpxo9.v5a.nc`
Constituent: M2
Spin-up excluded before harmonic fit: 1 M2 periods
Sponge crop W/S/N: 20 cells
Minimum depth for min-|η| search: 50 m

## Key interpretation columns

- `model_min_eta_amp_m`: how close the model gets to an SSH node.
- `model_min_winding_r8_rad`: values near ±6.283 imply amphidromic phase winding.
- `tpxo_at_model_min_eta_amp_m` and `tpxo_at_model_min_winding_r8_rad`: whether TPXO has a comparable feature at the model node.
- `model_to_tpxo_min_distance_km`: separation between model and TPXO low-|η| candidates.

## Compact result

| run | model |η| p95 (m) | model min lon,lat | model |η|min (m) | model r8 winding (rad) | TPXO |η| at model min (m) | TPXO r8 winding at model min | distance to TPXO min (km) | canyon model/TPXO |η| (m) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| runBC_LR_fullTPXO_balanced | 0.01936 | 242.9545, 31.6732 | 8.8e-05 | 6.283 | 0.4874 | 0 | 38.83 | 0.003935/0.4931 |
| runBC_LR_fullTPXO_obeta | 0.4877 | 242.5028, 32.2492 | 0.4618 | 0 | 0.4896 | 0 | 70.45 | 0.4748/0.4931 |
| runBC_LR_fullTPXO_obeta_balanced | 0.0313 | 243.2024, 31.5981 | 8.914e-05 | 6.283 | 0.4899 | 0 | 41.22 | 0.008206/0.4931 |
| runBC_LR_fullTPXO_obeta_weakvel025 | 0.007819 | 243.2024, 31.5981 | 2.122e-05 | -6.661e-16 | 0.4899 | 0 | 41.22 | 0.002054/0.4931 |
| runBC_LR_fullTPXO_obeta_weakvel010 | 0.003127 | 243.2024, 31.5981 | 8.432e-06 | -2.22e-16 | 0.4899 | 0 | 41.22 | 0.0008224/0.4931 |
| runBC_LR_fullTPXO_obeta_zerovel | 7.348e-09 | 242.9765, 31.6181 | 2.207e-11 | 6.283 | 0.4867 | 0 | 43.18 | 7.458e-09/0.4931 |

## Reading the result

If a TPXO-forced test has a near-zero model |η| with ~2π winding, but TPXO has large |η| and no winding at the same point, then the node is being produced by the model configuration rather than inherited from TPXO.
