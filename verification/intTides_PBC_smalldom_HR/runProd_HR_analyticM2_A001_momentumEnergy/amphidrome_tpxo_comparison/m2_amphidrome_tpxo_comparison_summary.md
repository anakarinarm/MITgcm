# M2 amphidrome TPXO comparison

Run directory: `/Users/karina/Tools/MITgcm/verification/intTides_PBC_smalldom_HR/runProd_HR_analyticM2_A001_momentumEnergy`
Diagnostics directory: `/Volumes/LaCie/MITgcm_outputs/intTides_PBC_smalldom_HR/runProd_HR_analyticM2_A001_momentumEnergy`
TPXO elevation: `/Users/karina/Research/Canyons/data/tides/h_tpxo9.v5a.nc`
TPXO velocity: `/Users/karina/Research/Canyons/data/tides/u_tpxo9.v5a.nc`
Constituent: M2
Samples used after spin-up: 40
Time range used: 134136.0 to 352107.0 s
Sponge cells cropped from W/S/N: 60
Minimum depth used for minimum-|η| amphidrome search: 50 m

Important caveat: this production run uses analytical velocity forcing with no prescribed TPXO SSH,
so absolute SSH amplitudes are not expected to match TPXO.  The most useful comparison here is
whether the model and TPXO share a low-|η| point with amphidromic phase winding, and whether the
phase geometry near the canyon is similar.

## Domain-wide model minus TPXO summary

- median |η| difference: -0.47531 m
- median absolute |η| difference: 0.47531 m
- median absolute phase difference: 0.834997 rad
- median speed-amplitude difference: -0.006688 m/s
- median absolute speed-amplitude difference: 0.0073053 m/s

## Amphidrome locations

- Model minimum |η|: lon 242.90492249, lat 31.41444588
- TPXO minimum |η|: lon 243.17488098, lat 31.96869469
- Model-min to TPXO-min distance: 66.596 km
- Model-min to canyon distance: 52.514 km
- TPXO-min to canyon distance: 21.584 km

## Point metrics

### model: model_min_eta_amp

- lon/lat: 242.90492249, 31.41444588
- depth: 1398.832 m
- |η|: 1.18547e-05 m
- η phase: 0.931596 rad
- barotropic/current speed amplitude: 0.0132394 m/s
- phase winding r=2/4/6/8: 1.776e-15, 8.882e-16, -2.665e-15, 6.283 rad

### model: tpxo_min_eta_amp

- lon/lat: 243.17488098, 31.96869469
- depth: 50.800 m
- |η|: 0.0156236 m
- η phase: -1.60784 rad
- barotropic/current speed amplitude: 0.00470856 m/s
- phase winding r=2/4/6/8: 1.776e-15, nan, nan, nan rad

### model: canyon_marker

- lon/lat: 243.25016785, 31.78505516
- depth: 281.100 m
- |η|: 0.0117311 m
- η phase: -1.44287 rad
- barotropic/current speed amplitude: 0.00588115 m/s
- phase winding r=2/4/6/8: 1.332e-15, -4.441e-16, 4.441e-16, -1.776e-15 rad

### tpxo: model_min_eta_amp

- lon/lat: 242.90492249, 31.41444588
- depth: 1398.832 m
- |η|: 0.482815 m
- η phase: -2.3647 rad
- barotropic/current speed amplitude: 0.0227018 m/s
- phase winding r=2/4/6/8: 0, 0, 0, 0 rad

### tpxo: tpxo_min_eta_amp

- lon/lat: 243.17488098, 31.96869469
- depth: 50.800 m
- |η|: 0.0925503 m
- η phase: -2.38617 rad
- barotropic/current speed amplitude: 0.0021356 m/s
- phase winding r=2/4/6/8: 0, nan, nan, nan rad

### tpxo: canyon_marker

- lon/lat: 243.25016785, 31.78505516
- depth: 281.100 m
- |η|: 0.493051 m
- η phase: -2.38117 rad
- barotropic/current speed amplitude: 0.0355191 m/s
- phase winding r=2/4/6/8: 0, 0, 0, 0 rad

## How to read this

- A true amphidromic point should combine very small |η| with phase winding near ±2π on rings around the point.
- If the model minimum |η| is close to the TPXO minimum |η| and has similar winding, the amphidrome is likely inherited from the regional barotropic tide.
- If the model minimum is displaced, much stronger/weaker, or TPXO has no comparable winding, the amphidrome is more likely produced by our boundary/geometry/sponge choices.
