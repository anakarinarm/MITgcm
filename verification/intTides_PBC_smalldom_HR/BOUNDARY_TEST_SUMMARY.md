# Boundary-condition test summary for `intTides_PBC_smalldom_HR`

Last updated: 2026-06-27

This document is a single “memory file” for the boundary-condition work done before turning on the full momentum-budget diagnostics. The practical goal was to find a boundary setup that lets us study internally generated tides/waves without the open boundaries dominating the solution through unbalanced barotropic transport, excessive SSH clamping, or strong artificial reflections.

## Short version

The cleanest direction so far is **analytical south-entering M2 velocity, no prescribed SSH, north-only balanced return, Stevens/radiative open boundaries, and a momentum sponge**.

The current preferred low-resolution test is:

- `runBC_LR_analyticM2_south_A001_Nreturn`

The corresponding high-resolution test prepared for the first HR run is:

- `runBC_HR_analyticM2_south_A001_Nreturn`

The setup to avoid for production is:

- full-boundary TPXO velocity plus prescribed SSH/elevation (`OB*etaFile`) in this very small domain.

The main lesson was not simply “TPXO is bad.” It was more subtle: TPXO is a dynamically balanced large-scale tidal solution, but once it is interpolated onto a small regional MITgcm domain and imposed on every open boundary, the boundary data can over-constrain the model. Prescribed SSH was especially problematic: it imposed a large standing barotropic/free-surface signal that contaminated the internal-wave diagnostics.

## Why we paused the momentum budget work

The original plan was to configure `runTSini_mmtm/` for vector-invariant momentum-budget diagnostics. Before doing that, we tested the boundary conditions because a momentum budget is only useful if the boundary forcing is physically and numerically sane.

Early tests showed unintended ETAN growth and strong boundary signals. That meant the boundary problem needed to be solved first; otherwise the momentum diagnostics would mostly diagnose boundary imbalance rather than internal-tide dynamics.

## Numerical/physical configuration changes made earlier

The production-style configuration was moved away from large constant viscosities/diffusivities and toward grid-aware values:

- Horizontal viscosity:
  - `viscAh = 0`
  - `viscAhGrid = 0.05`
  - `viscAhGridMax = 0.05`
- Vertical viscosity:
  - `viscAz = 1e-4 m2/s`
- Horizontal tracer diffusivity:
  - `diffKhT = 0`
  - `diffKhS = 0`
- Vertical tracer diffusivity:
  - `diffKzT = 1e-5 m2/s`
  - `diffKzS = 1e-5 m2/s`
- Bottom drag:
  - `bottomDragQuadratic = 2e-3`

The philosophy was: use enough damping for stability and unresolved shear, but avoid damping the internal tide with large arbitrary constants before sensitivity tests.

## Boundary-condition families tested

### 1. Original TPXO-derived boundary forcing

Initial boundary files were made by interpolating TPXO9 tidal velocities onto the MITgcm mesh. The first tests used W/S/N tidal velocity forcing with T/S held fixed at the boundaries.

What we found:

- The solution had unintended ETAN behavior.
- Boundary fluxes were not cleanly outward/radiative.
- Reflection metrics were hard to interpret because the imposed boundary tide itself dominated the near-boundary signal.

This led to the first boundary-only diagnostic runs.

### 2. No-sponge boundary test

Run family:

- `runBC_test`

Purpose:

- Minimal-output boundary test before enabling the full momentum diagnostics.

What we found:

- The no-sponge case showed problematic boundary behavior.
- It was useful as a baseline, but not acceptable as a production boundary setup.

### 3. Sponge tests

Run families:

- `runBC_sponge_test`
- `runBC_telescopic_sponge_test`
- `runBC_LR_noSSH_sponge_gentle`
- `runBC_LR_noSSH_sponge_strong`

Sponge variants included:

- gentle sponge:
  - boundary relaxation timescale: `11178 s`
  - interior relaxation timescale: `178848 s`
- strong sponge:
  - boundary relaxation timescale: `1397.25 s`
  - interior relaxation timescale: `44712 s`

What we found:

- Sponge layers helped somewhat, especially near the boundaries.
- Sponge strength alone did not solve the problem.
- The far-west section still had inward-dominated mode-1 energy in several TPXO-forced cases.

Conclusion:

- Sponge is useful, but it cannot rescue inconsistent or over-constrained tidal forcing.

### 4. Telescopic sponge / larger buffer-zone idea

Run families:

- `runBC_telescopic_test`
- `runBC_telescopic_sponge_test`

Purpose:

- Add a buffer zone using telescopic cells.
- Keep the physical domain intact.
- Extend bathymetry into the sponge without adding new bathymetric features.
- Treat the sponge as a numerical buffer, not as part of the science domain.

What we found:

- This made physical sense as a numerical strategy.
- However, the first telescopic+sponge setup produced unexpectedly large boundary fluxes and did not reduce amplification below the desired level.
- The result pushed us toward separating two questions:
  1. Is the sponge/buffer working?
  2. Is the prescribed tidal boundary forcing dynamically consistent?

Conclusion:

- Telescopic sponge is still a good strategy for the HR configuration, but forcing consistency matters more than the sponge geometry.

### 5. Full-boundary TPXO with prescribed SSH/elevation

Run families:

- `runBC_LR_fullTPXO_obeta`
- `runBC_LR_fullTPXO_obeta_balanced`
- `runBC_LR_fullTPXO_obeta_weakvel010`
- `runBC_LR_fullTPXO_obeta_weakvel025`
- `runBC_LR_fullTPXO_obeta_zerovel`
- `runBC_LR_fullTPXO_obeta_orlanski`

Purpose:

- Test whether prescribing TPXO-consistent SSH together with velocity would make the barotropic tide more dynamically balanced.

Important MITgcm limitation found:

- Stock MITgcm/OBCS does **not** support Orlanski/Stevens-style radiation together with nonlinear free surface and prescribed OBCS eta in this configuration.
- The failed run reported:

```text
OBCS_CHECK: nonlinFreeSurf not yet implemented in Orlanski OBC
```

What we found:

- Prescribed SSH created very large standing modal-characteristic signals.
- The most revealing case was SSH-only with velocity zero:
  - ETAN RMS: about `0.0407 m`
  - W mode-1 outward/inward at offset 22: `182.47 / 182.47 MW`
  - S mode-1 outward/inward at offset 22: `84.98 / 84.98 MW`
  - N mode-1 outward/inward at offset 22: `53.12 / 53.12 MW`
- Equal outward and inward characteristic energy is a signature of a standing imposed signal, not clean outgoing internal-wave radiation.

Conclusion:

- For this small domain, full-boundary SSH clamping is inappropriate for the production internal-tide/momentum-budget run.

### 6. No-SSH, balanced/radiative TPXO velocity

Run family:

- `runBC_LR_recommended_noSSH_balanced_radiative`

Architecture:

- no prescribed SSH;
- balanced W/S/N M2 normal velocities;
- Stevens/radiative W/S/N boundaries;
- momentum sponge;
- boundary diagnostics only.

What we found:

- ETAN RMS was much cleaner than the SSH cases, about `0.009 m`.
- But west mode-1 energy remained problematic:
  - W offset 33 mode-1 outward/inward: `5.51 / 10.51 MW`
  - incoming fraction: `0.656`

Conclusion:

- Removing SSH was a major improvement.
- Full W/S/N TPXO velocity forcing still appeared too constraining for the small domain.

### 7. No-west TPXO velocity

Run family:

- `runBC_LR_noSSH_noW_strongSponge`

Architecture:

- no prescribed SSH;
- west tidal velocity set to zero;
- S/N rebalanced;
- Stevens/radiative W/S/N;
- strong sponge.

What we found:

- West mode-1 contamination improved compared with full W/S/N TPXO forcing.
- W offset 33 mode-1 outward/inward became `3.48 / 5.45 MW`, instead of `5.70 / 10.55 MW` for the strong-sponge full-TPXO case.
- South outgoing mode-1 energy was reduced:
  - S offset 22 mode-1 outward/inward: `10.39 / 2.91 MW`
- ETAN RMS increased slightly to about `0.0115 m`.

Conclusion:

- Removing west forcing helped, but this still depended on TPXO-derived forcing and did not feel like the cleanest controlled experiment.

### 8. Analytical south-entering M2 forcing

Motivation:

- TPXO is a large-scale tidal solution. In this small regional model, imposing TPXO on all open boundaries can over-constrain the local dynamics.
- A controlled analytical tide entering from the south is less realistic as a tidal hindcast, but more useful as an idealized internal-tide experiment.
- The physical idea is plausible: the barotropic tide likely enters the regional domain dominantly from the south, consistent with coastal/Kelvin-wave intuition.

#### `A0 = 0.02 m/s`, north-only return

Run family:

- `runBC_LR_analyticM2_south_A002`

Architecture:

- no SSH;
- south normal M2 velocity with `A0 = 0.02 m/s`;
- 10-cell cosine taper near south corners;
- uniform north compensating return flow;
- west velocity zero;
- Stevens/radiative W/S/N;
- strong 20-cell sponge.

What we found:

- ETAN RMS: `0.0059 m`
- South mode-1 energy was strong:
  - S offset 22 outward/inward: `44.49 / 13.21 MW`
- This was likely too energetic for the first controlled baseline.

#### `A0 = 0.01 m/s`, north-only return

Run family:

- `runBC_LR_analyticM2_south_A001_Nreturn`

Architecture:

- no SSH;
- south normal M2 velocity with `A0 = 0.01 m/s`;
- 10-cell cosine taper near south corners;
- uniform north compensating return flow;
- west velocity exactly zero;
- Stevens/radiative W/S/N;
- strong 20-cell sponge.

Generated forcing diagnostics:

- south max/RMS amplitude: `0.0100 / 0.00968 m/s`
- north return max/RMS amplitude: `0.02145 / 0.01617 m/s`
- west max/RMS amplitude: `0.0 / 0.0 m/s`
- transport residual `Q_W + Q_S - Q_N`: `0.0 m3/s`

What we found:

- ETAN RMS: `0.00289 m`
- South internal tide remained energetic but not excessive:
  - S offset 22 mode-1 outward/inward: `11.41 / 3.72 MW`
- West was not perfect, but cleaner than many TPXO cases:
  - W offset 22 mode-1 outward/inward: `4.38 / 4.55 MW`
  - W offset 33 mode-1 outward/inward: `3.91 / 5.80 MW`
- North return was visible but not catastrophic:
  - N offset 22 mode-1 outward/inward: `0.41 / 0.26 MW`

Conclusion:

- This is the current best analytical baseline.
- It preserves a useful south-generated internal tide while keeping ETAN much cleaner than TPXO/SSH experiments.

#### `A0 = 0.01 m/s`, north+west return

Run family:

- `runBC_LR_analyticM2_south_A001_NWreturn`

Architecture:

- no SSH;
- south normal M2 velocity with `A0 = 0.01 m/s`;
- 10-cell cosine taper near south corners;
- compensating return distributed over north and west;
- west return represented by a half-period phase lag;
- Stevens/radiative W/S/N;
- strong 20-cell sponge.

Generated forcing diagnostics:

- south max/RMS amplitude: `0.0100 / 0.00968 m/s`
- north return max/RMS amplitude: `0.004996 / 0.00377 m/s`
- west return max/RMS amplitude: `0.004996 / 0.004996 m/s`
- west return phase lag: `22356 s`, half an M2 period
- transport residual `Q_W + Q_S - Q_N`: `0.0 m3/s`

What we found:

- ETAN RMS was the cleanest of the analytical family: `0.00185 m`
- North-boundary contamination was strongly reduced:
  - N offset 22 mode-1 outward/inward: `0.089 / 0.086 MW`
- But the south internal tide was also weaker:
  - S offset 22 mode-1 outward/inward: `5.20 / 1.39 MW`
- West far-section mode 1 remained inward dominated:
  - W offset 33 mode-1 outward/inward: `1.38 / 2.63 MW`

Conclusion:

- North+west return is excellent for barotropic cleanliness, but may under-force the internal tide at `A0 = 0.01 m/s`.
- For the first HR test, the north-only return is a better compromise.

## Key numeric comparison

This table summarizes the most useful scalar diagnostics from the low-resolution experiments. Fluxes are shown in MW. Modal values are mode-1 outward/inward characteristic energy rates.

| run | ETAN RMS m | W flux22 MW | W flux33 MW | W m1 22 out/in MW frac | W m1 33 out/in MW frac | S m1 22 out/in MW frac | N m1 22 out/in MW frac |
|---|---:|---:|---:|---:|---:|---:|---:|
| LR no-SSH full TPXO balanced | 0.0090 | 0.84 | -4.47 | 7.11/7.75 (0.522) | 5.51/10.51 (0.656) | 19.84/4.07 (0.170) | 0.38/0.22 (0.372) |
| LR strong sponge | 0.0090 | 1.18 | -4.18 | 7.26/7.63 (0.512) | 5.70/10.55 (0.649) | 20.58/4.58 (0.182) | 0.38/0.22 (0.362) |
| LR gentle sponge | 0.0090 | 0.89 | -4.24 | 7.05/7.64 (0.520) | 5.37/10.10 (0.653) | 19.28/3.59 (0.157) | 0.37/0.23 (0.384) |
| LR no-west TPXO balanced | 0.0115 | 1.01 | -1.58 | 4.00/3.97 (0.498) | 3.48/5.45 (0.610) | 10.39/2.91 (0.219) | 0.46/0.22 (0.326) |
| LR SSH zero velocity | 0.0407 | -0.00 | 0.00 | 182.47/182.47 (0.500) | 169.50/169.50 (0.500) | 84.98/84.98 (0.500) | 53.12/53.12 (0.500) |
| LR SSH weak 10% velocity | 0.0408 | 0.15 | 0.19 | 182.68/182.62 (0.500) | 169.74/169.64 (0.500) | 85.02/85.18 (0.500) | 53.21/53.26 (0.500) |
| LR SSH weak 25% velocity | 0.0409 | 0.96 | 1.16 | 183.81/183.41 (0.499) | 170.95/170.36 (0.499) | 85.18/86.23 (0.503) | 53.68/53.98 (0.501) |
| LR SSH full balanced velocity | 0.0425 | 15.28 | 18.71 | 203.92/197.51 (0.492) | 192.64/183.10 (0.487) | 88.12/104.93 (0.544) | 62.09/66.91 (0.519) |
| LR analytic A002 N return | 0.0059 | 4.39 | -4.63 | 17.46/17.22 (0.497) | 15.71/22.11 (0.585) | 44.49/13.21 (0.229) | 1.67/1.03 (0.381) |
| LR analytic A001 N return | 0.0029 | 0.84 | -1.45 | 4.38/4.55 (0.509) | 3.91/5.80 (0.597) | 11.41/3.72 (0.246) | 0.41/0.26 (0.390) |
| LR analytic A001 NW return | 0.0018 | 0.22 | -1.09 | 1.76/1.95 (0.526) | 1.38/2.63 (0.656) | 5.20/1.39 (0.211) | 0.09/0.09 (0.491) |

Interpretation of the table:

- The SSH cases are immediately suspicious because they create huge equal in/out modal signals, even with zero velocity forcing.
- The no-SSH TPXO cases are much cleaner in ETAN, but still show west-boundary inward contamination.
- The analytical `A0 = 0.01 m/s`, north-only return case gives a good balance: low ETAN, useful south-generated internal tide, and no prescribed west tide.
- The analytical north+west return case is barotropically very clean, but likely too weak for the first production-style internal-tide baseline.

## What the Python diagnostics compute

Two diagnostic scripts were used repeatedly:

- `check_boundary_reflections.py`
- `modal_characteristic_decomposition.py`

They are complementary. The first gives a bulk section pressure-work/reflection proxy. The second tries to split the boundary signal into vertical-mode inward and outward characteristic components.

### `check_boundary_reflections.py`

Purpose:

- Diagnose internal-wave reflection at open boundaries from model output.
- Provide fast scalar summaries for comparing boundary configurations.

Input diagnostics:

- `bc_wave_3d.*`
  - `UVEL`
  - `VVEL`
  - `PHIHYD`
- `bc_wave_eta.*`
  - `ETAN`
- grid files:
  - `DRF`
  - `hFacC`
  - `DXC`
  - `DYC`
  - `RAC`

Important implementation detail:

- Earlier, ETAN and 3-D fields were accidentally placed in the same diagnostics stream. MITgcm then wrote only their common surface level, causing the first script error:

```text
expected 1080000 values, got 64800000
```

That was fixed by separating 3-D and 2-D diagnostics into `bc_wave_3d` and `bc_wave_eta`.

Algorithm:

1. Read all post-ramp boundary diagnostics.
2. Use only times after the ramp:
   - default M2 period: `44712 s`
   - default ramp cutoff: `3` M2 periods
3. Fit the M2 harmonic at each boundary section using a least-squares model with:
   - constant term;
   - linear trend;
   - cosine M2 term;
   - sine M2 term.
4. Reconstruct complex M2 amplitudes for pressure and normal velocity.
5. Form pressure as:

```text
pressure potential = PHIHYD + g * ETAN
```

6. Remove the depth-mean pressure and normal velocity to focus on the baroclinic component.
7. Compute outward baroclinic pressure-work flux on sections parallel to the open boundaries.

Boundary sign convention:

- `W`: outward velocity is westward, so the script uses `-U`.
- `S`: outward velocity is southward, so the script uses `-V`.
- `N`: outward velocity is northward, so the script uses `+V`.

Important output columns:

- `boundary`
  - `W`, `S`, or `N`.
- `offset_cells`
  - number of cells inward from the boundary section.
- `distance_km`
  - physical distance from the open boundary to the section. This was corrected for telescopic grids by summing actual grid spacing instead of assuming uniform `dx`/`dy`.
- `outward_flux_W`
  - outward baroclinic pressure-work flux.
  - Positive means net outward energy flux through that section.
  - Negative means net inward-directed flux.
- `progressive_index`
  - normalized pressure/velocity covariance:

```text
chi = Re(p u*) / sqrt(|p|^2 |u|^2)
```

  - `chi` near `+1`: strongly outward/progressive.
  - `chi` near `0`: standing or nearly quadrature pressure/velocity.
  - `chi < 0`: inward-dominated by this bulk metric.
- `reflection_energy_proxy`
  - under a local, single-mode, equal-impedance approximation:

```text
R_E = (1 - chi) / (1 + chi)
```

  - only computed when `chi > 0`.
- `reflection_amplitude_proxy`
  - square root of the energy proxy:

```text
R_A = sqrt(R_E)
```

Limitations:

- This is **not** a true multimode reflection coefficient.
- It is a bulk baroclinic section diagnostic.
- It can be badly contaminated by imposed SSH standing signals.
- It is most useful for comparing cases and checking consistency across nearby offsets.

### `modal_characteristic_decomposition.py`

Purpose:

- Separate boundary-normal M2 energy into approximate outward and inward modal characteristics.
- Help distinguish an outgoing internal tide from an imposed incoming/standing boundary signal.

Input diagnostics:

- Same model diagnostics as `check_boundary_reflections.py`:
  - `UVEL`
  - `VVEL`
  - `PHIHYD`
  - `ETAN`
- grid files:
  - `DRF`
  - `RC`
  - `hFacC`
  - `DXC`
  - `DYC`
  - `RAC`
- initial hydrography files from `data`:
  - `hydrogThetaFile`
  - `hydrogSaltFile`
- linear EOS parameters:
  - `tAlpha`
  - `sBeta`

Algorithm:

1. Read the same post-ramp M2 boundary sections.
2. Fit the M2 harmonic for pressure and normal velocity.
3. Remove the depth mean to isolate the baroclinic component.
4. Estimate stratification from the initial horizontally uniform T/S profiles using the model linear EOS.
5. Build simple WKB/cosine hydrostatic vertical modes.
6. Estimate long-wave modal phase speeds:

```text
c_n ≈ integral(N dz) / (n pi)
```

7. Project pressure and normal velocity onto modes.
8. For each mode, split into outward and inward characteristic velocity-like amplitudes:

```text
u_out = 0.5 * (u + p/c)
u_in  = 0.5 * (p/c - u)
```

9. Convert those amplitudes into outward and inward characteristic energy-rate proxies.

Important output columns:

- `boundary`
  - `W`, `S`, or `N`.
- `offset_cells`
  - section offset from the boundary.
- `distance_km`
  - physical section distance from the boundary.
- `mode`
  - vertical mode number. We generally inspected modes 1-5, with mode 1 used for the headline comparisons.
- `phase_speed_m_s`
  - estimated WKB long-wave speed for that mode.
- `outward_characteristic_W`
  - outward characteristic energy-rate proxy.
- `inward_characteristic_W`
  - inward characteristic energy-rate proxy.
- `net_outward_W`
  - outward minus inward.
- `inward_to_outward_amplitude`
  - approximate amplitude ratio:

```text
sqrt(inward / outward)
```

- `incoming_energy_fraction`
  - fraction of modal characteristic energy that is inward:

```text
inward / (outward + inward)
```

Interpretation:

- `incoming_energy_fraction < 0.5`
  - outward dominated.
- `incoming_energy_fraction ≈ 0.5`
  - standing/equal in-out signal.
- `incoming_energy_fraction > 0.5`
  - inward dominated.

Limitations:

- This is an approximate diagnostic decomposition, not a formal modal radiation boundary condition.
- The vertical modes use a simplified WKB/cosine approximation.
- Sections close to forced boundaries can still include the imposed barotropic/tidal signal.
- The absolute numbers should be treated with care; the relative comparison among runs is the most robust use.

## Why prescribed SSH was rejected

The SSH tests were extremely informative. If velocity is zero but SSH is prescribed, the internal-wave diagnostic should ideally be quiet. Instead, the SSH-only case produced huge equal outward/inward modal signals at all boundaries.

That means the prescribed elevation field was not simply “balancing” the velocity. It was acting like a strong imposed standing free-surface condition. In a small domain, this over-constrained the barotropic tide and polluted the baroclinic/modal diagnostics.

The result is why the recommended path is:

- no `OB*etaFile` for now;
- balanced normal velocity forcing;
- radiative/Stevens behavior where possible;
- sponge/buffer to absorb outgoing waves;
- analytical forcing for controlled process experiments.

## Why analytical south forcing is scientifically reasonable

The analytical south-entering tide is less realistic than TPXO as a tidal hindcast, but it is better matched to the scientific question: internal-tide generation and momentum budgets in a controlled regional setup.

Reasons it is worth using:

- It avoids imposing a coarse global/regional tidal solution on all boundaries of a small high-resolution domain.
- It gives clear control over amplitude, direction, and phasing.
- It lets us test sensitivity to:
  - tide strength;
  - incoming direction;
  - spring/neap modulation;
  - extra constituents such as K1.
- It keeps the boundary forcing interpretable when diagnosing momentum budgets.

Important caveat:

- It should be described as an **idealized process experiment**, not as a direct TPXO-forced hindcast.

## Current HR test setup

Run family:

- `runBC_HR_analyticM2_south_A001_Nreturn`

Architecture:

- copied from `runBC_telescopic_sponge_test`;
- no prescribed SSH/elevation;
- south normal M2 velocity with `A0 = 0.01 m/s`;
- 30-cell cosine taper near south corners, physically comparable to the 10-cell LR taper;
- north uniform compensating return flow;
- west prescribed velocity exactly zero;
- Stevens/radiative W/S/N;
- 60-cell strong momentum sponge;
- boundary/reflection diagnostics only, no momentum-budget diagnostics yet.

Generated forcing diagnostics:

- grid: `660 x 720 x 60`
- south max/RMS amplitude: `0.0100 / 0.00970 m/s`
- north return max/RMS amplitude: `0.01848 / 0.01390 m/s`
- west max/RMS amplitude: `0.0 / 0.0 m/s`
- transport residual `Q_W + Q_S - Q_N`: `2.3e-10 m3/s`

Run/analyze commands:

```bash
mpiexec -n 16 ./mitgcmuv
python3 check_boundary_reflections.py
python3 modal_characteristic_decomposition.py
```

What to check after running:

- ETAN should remain small and unbiased.
- South should show a clear outgoing internal-tide signal.
- West should be less contaminated than TPXO/SSH cases.
- North return should not dominate the modal diagnostics.
- Compare HR against `runBC_LR_analyticM2_south_A001_Nreturn`, but expect differences because the HR grid and bathymetry resolve more internal-wave structure.

## First HR analytical test result

Run family:

- `runBC_HR_analyticM2_south_A001_Nreturn`

Status:

- The run completed normally.
- STDERR only showed the harmless retired `EXACT_CONSERV` warning.
- The boundary diagnostic scripts produced:
  - `boundary_reflection_summary.csv`
  - `modal_characteristic_summary.csv`
  - `boundary_reflection_summary.png`

Final ETAN snapshot:

| run | ETAN min m | ETAN max m | ETAN mean m | ETAN RMS m |
|---|---:|---:|---:|---:|
| HR analytical A001 N return | -0.0372 | 0.00819 | -0.00674 | 0.0122 |
| LR analytical A001 N return | -0.0392 | 0.00065 | -0.0115 | 0.0159 |

The final-snapshot ETAN RMS is not exploding in HR. In fact, by this crude snapshot metric HR is slightly cleaner than the LR baseline. This suggests the HR issue is not a catastrophic barotropic mass-balance problem.

Boundary/reflection diagnostic comparison:

| run | boundary | offset cells | distance km | bulk outward flux MW | progressive index |
|---|---|---:|---:|---:|---:|
| HR analytical A001 N return | W | 65 | 24.41 | -2.03 | -0.093 |
| HR analytical A001 N return | W | 70 | 25.28 | -2.54 | -0.113 |
| HR analytical A001 N return | W | 80 | 27.01 | -3.55 | -0.149 |
| HR analytical A001 N return | W | 100 | 30.49 | -4.74 | -0.184 |
| HR analytical A001 N return | S | 65 | 26.06 | 10.37 | 0.489 |
| HR analytical A001 N return | S | 70 | 26.98 | 10.55 | 0.490 |
| HR analytical A001 N return | S | 80 | 28.84 | 10.84 | 0.489 |
| HR analytical A001 N return | S | 100 | 32.55 | 11.21 | 0.527 |
| HR analytical A001 N return | N | 65 | 25.27 | 1.70 | 0.338 |
| HR analytical A001 N return | N | 70 | 26.20 | 1.75 | 0.342 |
| HR analytical A001 N return | N | 80 | 28.05 | 1.88 | 0.346 |
| HR analytical A001 N return | N | 100 | 31.77 | 1.91 | 0.356 |

Mode-1 characteristic comparison:

| run | boundary | offset cells | mode-1 outward MW | mode-1 inward MW | incoming energy fraction |
|---|---|---:|---:|---:|---:|
| HR analytical A001 N return | W | 65 | 4.16 | 5.63 | 0.575 |
| HR analytical A001 N return | W | 70 | 4.21 | 5.92 | 0.585 |
| HR analytical A001 N return | W | 80 | 4.20 | 6.43 | 0.605 |
| HR analytical A001 N return | W | 100 | 4.39 | 7.51 | 0.631 |
| HR analytical A001 N return | S | 65 | 9.53 | 1.60 | 0.144 |
| HR analytical A001 N return | S | 70 | 9.68 | 1.53 | 0.137 |
| HR analytical A001 N return | S | 80 | 10.04 | 1.39 | 0.122 |
| HR analytical A001 N return | S | 100 | 9.74 | 1.00 | 0.093 |
| HR analytical A001 N return | N | 65 | 0.85 | 0.36 | 0.296 |
| HR analytical A001 N return | N | 70 | 0.88 | 0.36 | 0.289 |
| HR analytical A001 N return | N | 80 | 0.95 | 0.35 | 0.267 |
| HR analytical A001 N return | N | 100 | 0.91 | 0.31 | 0.253 |

Interpretation:

- The south boundary looks good. It has a stable positive bulk outward flux of about `10-11 MW`, and mode-1 is clearly outward dominated. The incoming mode-1 fraction decreases with distance from the boundary, from `0.144` to `0.093`.
- The north boundary is modest and also outward dominated. The return flow is visible but not dominating the modal diagnostics.
- The west boundary remains the main problem. Bulk flux is negative at all tested offsets and becomes more negative farther inside the domain. Mode-1 incoming fraction increases from `0.575` to `0.631`.
- Because west velocity forcing is exactly zero in this experiment, the west signal is probably not caused by imposed west tidal forcing. It is more likely due to internally generated wave energy interacting with the west open boundary/sponge/bathymetry, or a modal decomposition artifact where the local wave field is not well represented by the simple section-normal modal split.

Decision from this HR test:

- The analytical south-forcing architecture is promising and much cleaner than the TPXO+SSH family.
- It is not quite ready to declare the west boundary solved.
- Before turning on expensive momentum diagnostics, the next useful check is to diagnose why the west boundary is inward dominated in HR:
  - inspect spatial maps/sections of M2 pressure/velocity near the west buffer;
  - compare offsets farther inside/outside the sponge if available;
  - try a modest west sponge/relaxation adjustment;
  - or test a wider/stronger west-only sponge while keeping the south forcing unchanged.

### Follow-up HR gentle-sponge test prepared

Run family:

- `runBC_HR_analyticM2_south_A001_Nreturn_spongeGentle`

Purpose:

- Keep the analytical tide and HR grid unchanged.
- Test whether the west inward signal is partly caused by a too-strong/too-abrupt sponge.

Only intended change relative to `runBC_HR_analyticM2_south_A001_Nreturn`:

| parameter | baseline strong sponge | gentle-sponge test |
|---|---:|---:|
| `spongeThickness` | 60 cells | 60 cells |
| `Urelaxobcsbound` | `1397.25 s` | `11178 s` |
| `Urelaxobcsinner` | `44712 s` | `178848 s` |
| `Vrelaxobcsbound` | `1397.25 s` | `11178 s` |
| `Vrelaxobcsinner` | `44712 s` | `178848 s` |

Notes:

- MITgcm/OBCS exposes per-side sponge switches, but the relaxation timescale values are shared across sponge boundaries in the standard namelist. A truly west-only timescale would require code modification.
- `useLinearSponge` remains `.FALSE.`, preserving the smoother blended sponge target used by the existing OBCS sponge code.
- If west inward energy weakens, the baseline strong sponge was probably part of the reflection problem.
- If west inward energy is unchanged or worse, the next suspect is west-side bathymetry / oblique internal-wave structure / diagnostic geometry rather than sponge strength alone.

Result:

- The gentle-sponge run completed normally.
- Final ETAN RMS was very similar/slightly cleaner than the baseline:
  - strong sponge: `0.01218 m`
  - gentle sponge: `0.01182 m`
- West boundary diagnostics got worse:
  - W offset 65 bulk flux changed from `-2.03 MW` to `-3.17 MW`;
  - W offset 100 bulk flux changed from `-4.74 MW` to `-5.44 MW`;
  - W mode-1 incoming fraction at offset 65 changed from `0.575` to `0.617`;
  - W mode-1 incoming fraction at offset 100 changed from `0.631` to `0.667`.
- South became slightly cleaner by the modal diagnostic:
  - S offset 65 mode-1 incoming fraction changed from `0.144` to `0.101`.
- North remained modest.

Interpretation:

- The west problem is not caused by the strong sponge being too abrupt.
- Weakening the sponge lets more inward west signal remain.
- The original strong-sponge HR analytical run is still the better baseline.
- The next targeted test should either increase effective west absorption or isolate the west sponge by disabling S/N sponge while keeping strong W damping.

### Follow-up HR west-only strong-sponge test completed

Run family:

- `runBC_HR_analyticM2_south_A001_Nreturn_WspongeOnly`

Purpose:

- Keep the analytical tide, grid, diagnostics, and strong sponge timescales unchanged.
- Activate the momentum sponge only on the west boundary.
- Leave south and north with Stevens/radiative behavior and the same prescribed tide/return, but no local momentum sponge.

Only intended sponge-side-switch change relative to `runBC_HR_analyticM2_south_A001_Nreturn`:

| parameter | baseline strong sponge | west-only sponge test |
|---|---|---|
| `OBCSsponge_W` | `.TRUE.` | `.TRUE.` |
| `OBCSsponge_S` | `.TRUE.` | `.FALSE.` |
| `OBCSsponge_N` | `.TRUE.` | `.FALSE.` |
| `OBCSsponge_E` | `.FALSE.` | `.FALSE.` |

The sponge timescales remain:

- `spongeThickness = 60`
- `Urelaxobcsbound = 1397.25 s`
- `Urelaxobcsinner = 44712 s`
- `Vrelaxobcsbound = 1397.25 s`
- `Vrelaxobcsinner = 44712 s`

Result:

| run | ETAN RMS m | W flux65 MW | W flux100 MW | W m1 65 out/in MW frac | W m1 100 out/in MW frac | S m1 65 out/in MW frac | N m1 65 out/in MW frac |
|---|---:|---:|---:|---:|---:|---:|---:|
| HR strong all-sponge | 0.0030 | -2.03 | -4.74 | 4.16/5.63 (0.575) | 4.39/7.51 (0.631) | 9.53/1.60 (0.144) | 0.85/0.36 (0.296) |
| HR gentle all-sponge | 0.0027 | -3.17 | -5.44 | 3.68/5.94 (0.617) | 3.76/7.54 (0.667) | 9.35/1.06 (0.101) | 0.72/0.29 (0.288) |
| HR W-sponge-only | 0.0029 | -2.15 | -4.70 | 4.05/6.04 (0.599) | 4.19/7.64 (0.646) | 9.69/1.13 (0.104) | 0.83/0.40 (0.325) |

Interpretation:

- The W-sponge-only run kept ETAN clean:
  - final ETAN RMS was about `0.0029 m`.
- South remained a strong outward/progressive internal-tide source:
  - S offset 65 mode-1 outward/inward was `9.69 / 1.13 MW`;
  - incoming fraction was `0.104`.
- North remained modest:
  - N offset 65 mode-1 outward/inward was `0.83 / 0.40 MW`.
- West did **not** improve relative to the original strong all-sponge case:
  - W offset 65 bulk flux was `-2.15 MW`, very close to the original `-2.03 MW`;
  - W offset 100 bulk flux was `-4.70 MW`, very close to the original `-4.74 MW`;
  - W offset 65 mode-1 incoming fraction increased slightly from `0.575` to `0.599`;
  - W offset 100 mode-1 incoming fraction increased slightly from `0.631` to `0.646`.

Conclusion:

- The west inward signal is not mainly caused by damping on the south/north boundaries.
- Changing which boundaries receive sponge is not enough to fix the west diagnostic.
- The original strong all-sponge HR analytical run remains at least as good as W-sponge-only, and slightly cleaner at the west boundary.
- The next useful direction is probably not more side-switch sponge tuning; it is either:
  - change the return-flow design, e.g. a carefully scaled north+west return;
  - reduce `A0`;
  - move the west diagnostic section farther from the open boundary / physical domain transition;
  - or accept that west is receiving outgoing/recirculated internal-tide energy from the domain geometry and avoid using it as the primary “reflection score.”

## Recommendation before moving to momentum diagnostics

1. Run `runBC_HR_analyticM2_south_A001_Nreturn`.
2. Execute:

```bash
python3 check_boundary_reflections.py
python3 modal_characteristic_decomposition.py
```

3. If HR diagnostics look qualitatively consistent with the LR analytical baseline, promote this setup to the production momentum-budget run.
4. Only then turn on the heavier vector-invariant momentum diagnostics.

If the HR test shows unexpected west or north contamination, the next most useful sensitivity tests would be:

- lower analytical amplitude, e.g. `A0 = 0.0075 m/s`;
- slightly longer/stronger sponge;
- north+west return at a slightly higher amplitude than the LR `A0 = 0.01` NW-return case;
- boundary-isolation runs with only one boundary active at a time.

## Practical conclusions to remember

- Do not prescribe full TPXO SSH on all boundaries for this small-domain internal-tide experiment.
- No-SSH velocity forcing is much cleaner.
- Sponge helps, but it cannot fix over-constrained boundary data.
- Telescopic buffer zones make sense, but forcing consistency is the controlling issue.
- The best current controlled experiment is analytical south M2 forcing with balanced north return.
- The diagnostic scripts are proxies, not exact reflection calculators; use them comparatively.
- For the momentum budget, boundary cleanliness matters more than tidal realism at this stage.
