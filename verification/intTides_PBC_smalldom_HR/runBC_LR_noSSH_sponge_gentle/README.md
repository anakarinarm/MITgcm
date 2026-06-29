# LR no-SSH balanced/radiative test: gentle sponge

This is a one-change tuning experiment from
`runBC_LR_recommended_noSSH_balanced_radiative`.

Everything is unchanged except the momentum sponge relaxation time:

| setting | recommended baseline | this run |
|---|---:|---:|
| `spongeThickness` | 20 cells | 20 cells |
| `U/Vrelaxobcsbound` | 2794.5 s | 11178 s |
| `U/Vrelaxobcsinner` | 89424 s | 178848 s |

Interpretation:

- If west inward modal energy decreases, the baseline sponge was likely too
  abrupt/reflective.
- If west inward modal energy increases, the sponge was likely too weak to
  absorb outgoing disturbances.

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

Compared with the baseline, the gentle sponge produced only very small changes.
West mode-1 behavior:

| offset | outward/inward MW | incoming fraction |
|---:|---:|---:|
| 22 | 7.05 / 7.64 | 0.520 |
| 23 | 6.92 / 7.81 | 0.530 |
| 27 | 6.25 / 8.64 | 0.580 |
| 33 | 5.37 / 10.10 | 0.653 |

Final ETAN RMS was `0.00897 m`, essentially unchanged from the baseline
(`0.00895 m`).

Interpretation: weakening/smoothing the 20-cell sponge does not fix the
west-boundary inward modal-energy excess.  It helps very slightly at the far
offset but not enough to change the decision.
