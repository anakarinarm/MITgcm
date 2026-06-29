# LR no-SSH balanced/radiative test: stronger sponge

This is a one-change tuning experiment from
`runBC_LR_recommended_noSSH_balanced_radiative`.

Everything is unchanged except the momentum sponge relaxation time:

| setting | recommended baseline | this run |
|---|---:|---:|
| `spongeThickness` | 20 cells | 20 cells |
| `U/Vrelaxobcsbound` | 2794.5 s | 1397.25 s |
| `U/Vrelaxobcsinner` | 89424 s | 44712 s |

Interpretation:

- If west inward modal energy decreases, the baseline sponge was likely
  under-absorbing outgoing disturbances.
- If west inward modal energy increases, the sponge is likely too strong or too
  abrupt and is reflecting energy back into the physical domain.

Run and analyze:

```bash
mpiexec -n 16 ./mitgcmuv
python3 check_boundary_reflections.py
python3 modal_characteristic_decomposition.py
```

Compare against `../runBC_LR_recommended_noSSH_balanced_radiative`.

## Result

The run completed normally; STDERR only contained the harmless retired
`EXACT_CONSERV` warning.

Compared with the baseline, the stronger sponge modestly improved west mode 1
near the boundary, but did not remove the far-field inward dominance.

West mode-1 behavior:

| offset | outward/inward MW | incoming fraction |
|---:|---:|---:|
| 22 | 7.26 / 7.63 | 0.512 |
| 23 | 7.13 / 7.82 | 0.523 |
| 27 | 6.49 / 8.80 | 0.575 |
| 33 | 5.70 / 10.55 | 0.649 |

Final ETAN RMS was `0.00902 m`, essentially unchanged from the baseline
(`0.00895 m`).

Interpretation: stronger sponge is slightly preferable to the baseline for the
near-west sections, but sponge strength alone is not the dominant control on
the west-side inward modal-energy excess.
