# Boundary-condition validation: telescopic 60-cell sponge

This is the sponge member of a controlled telescopic-grid pair. Its no-sponge
control is `../runBC_telescopic_test`. The original-grid no-sponge and weak
30-cell-sponge tests remain unchanged in `../runBC_test` and
`../runBC_sponge_test`.

## Design

The complete original 600 x 600 physical domain is retained bit-for-bit at
Python indices `j=60:660, i=60:660` (MITgcm indices `j=61:660, i=61:660`).
The grid adds 60 cells west and 60 cells at both the south and north. There is
no eastern buffer because the east side is a closed land boundary. The new
grid is 660 x 720 x 60.

Cell widths grow geometrically by 2.5% per cell moving from the old domain to
the outer boundary. The transition cell has the same width as its old-domain
neighbor. The buffer widths are about 24 km zonally and 25 km meridionally;
the precise metric distance varies weakly with latitude. Adjacent cell-size
ratios never exceed 1.025.

The bathymetry is a normal extrusion of the old boundary values. Consequently
the buffer contains no new cross-boundary slope, canyon, or seamount. At the
southwest and northwest buffer corners, the corresponding old-domain corner
depth is extended. The horizontally uniform initial T/S profiles and the
constant boundary T/S records are extended identically.

The M2 forcing is not resampled from TPXO at the displaced numerical boundary.
The old complex harmonics are retained over their original coordinate spans.
On the new southwest and northwest corner segments they taper from the old
corner value to zero using a raised cosine, avoiding duplicated full-strength
forcing on two sides of each added corner. A minimum uniform correction over
wet faces is subsequently added, so the corrected outer-corner amplitude is
about 0.000376 m/s rather than exactly zero. The correction restores zero
complex transport and is computed from MITgcm's
initialized stretched-grid `hFacW`, `hFacS`, `DXG`, `DYG`, and `DRF`, rather
than extrapolated old-grid areas. The correction is 0.000375877 m/s and the
remaining complex transport is 2.4e-10 m3/s.

## Sponge

The momentum-only sponge fills the 60 added cells. Its nominal boundary
relaxation time is M2/16 = 2794.5 s and its inner nominal time is 2 M2 =
89424 s. With `useLinearSponge=.FALSE.`, the relaxation target also blends
into the instantaneous interior velocity and therefore the effective damping
goes to zero at the inner edge. Tracer sponge switches remain disabled because
the earlier tracer-sponge test became unstable; Stevens handles boundary T/S.

## Generated inputs

Preparation is deliberately two-stage because exact transport balancing needs
the grid metrics and partial cells produced by MITgcm. For a new grid, first
write the geometry and state inputs:

```bash
python prepare_telescopic_buffer.py --grid-only
```

Then perform a zero-step initialization. Finally rerun without `--grid-only` to
write the exactly balanced tidal files:

```bash
python prepare_telescopic_buffer.py
```

The files are uniquely tagged in `../input`; baseline inputs are not replaced.

## Compile, initialize, run, and analyze

The changed horizontal size requires the separate `../code_telescopic` and
`../build_telescopic` executable. It retains the current 4 x 4 MPI layout,
using 55 x 60 cells per tile and three tiles per process in each direction.

```bash
mkdir -p ../build_telescopic
bash ../code_telescopic/compila.sh
ln -sf ../build_telescopic/mitgcmuv mitgcmuv
mpiexec -n 16 ./mitgcmuv
python check_boundary_reflections.py
```

The checker samples 5, 10, 20, and 40 old-grid cells inward from the inner
sponge edge (outer-boundary offsets 65, 70, 80, and 100). Its reported distance
is the cumulative MITgcm metric distance, not cell count times local spacing.
As before, the pressure/velocity result is an approximate bulk reflection
proxy; a vertical-mode decomposition remains preferable for final assessment.

## Initialization verification

The dedicated executable compiled successfully and a zero-step initialization
on 16 MPI ranks ended normally. MITgcm read `Nx=660`, `Ny=720`, the complete
boundary files, and the 60-cell sponge without errors. A direct comparison of
the generated MITgcm grid confirms zero float32 difference in XC and YC over
the embedded 600 x 600 domain, zero float64 bathymetry difference there, and a
maximum adjacent-cell width ratio of exactly 1.025 in the buffers. `endTime` was
then restored to eight M2 periods for the actual test.

## Corrected completed test result

The corrected eight-M2 telescopic-sponge run completed normally. The post-ramp
domain-mean ETAN M2 amplitude is 0.000873 m, a large improvement over the
invalid first telescopic setup (0.008976 m) and smaller than the old-grid
30-cell sponge test (0.001935 m).

Relative to the no-sponge telescopic control, the 60-cell sponge modestly
improves the south and north progressive-wave indices but does not fix the
west side. In old-domain coordinates, the 5--40-cell sections give westward
outward flux of -3.9 to -9.0 MW and negative progressive-wave index, so the
bulk reflection proxy remains undefined there. South-side reflection-amplitude
proxy improves from about 0.49--0.53 without sponge to 0.46--0.51 with sponge.
North-side proxy improves only slightly, from about 0.73--0.74 to 0.71--0.72.

Interpretation: the balanced forcing and buffer removed the artificial ETAN
growth, and the sponge helps a little, but the limiting problem is now the west
boundary. Further improvement likely requires changing the west forcing/radiation
strategy or diagnosing vertical modes there; simply strengthening this bulk
sponge is unlikely to be the cleanest fix.

## Modal characteristic diagnostic

`modal_characteristic_decomposition.py` was added to split the fitted M2 signal
into approximate hydrostatic vertical modes and outward/inward characteristics.
It uses the linear-EOS initial T/S profile because `RhoRef` is constant in this
configuration. The estimated WKB long-wave speeds are c1--c5 = 3.397, 1.699,
1.132, 0.849, and 0.679 m/s.

For the completed telescopic-sponge run, the west boundary is inward-dominated
in modal characteristic energy. Summing modes 1--5 gives, at old-domain offsets
5 and 40 cells, respectively:

- west: outward/inward/net = 17.2/21.8/-4.6 MW and 14.3/25.4/-11.2 MW
- south: outward/inward/net = 31.0/7.9/23.0 MW and 32.3/6.0/26.4 MW
- north: outward/inward/net = 2.9/1.5/1.5 MW and 3.7/1.7/2.0 MW

The west signal is mostly mode-1 inward excess; mode 1 alone is about
5.7 MW outward and 9.0 MW inward at the 5-cell old-domain section. The
decomposition is approximate, but it confirms that the west issue is not just
an artifact of the bulk pressure-work scalar.

## Invalid first telescopic run retained in the history

The first completed telescopic-sponge run must not be used to evaluate the
sponge. Its tide-preparation step balanced transport using extrapolated
old-grid face areas and extended full-strength corner forcing. Recalculation
with the actual initialized grid found a residual M2 transport of 10311.8 m3/s.
Consistently, domain-mean ETAN amplitude rose to 0.008976 m, westward flux was
inward (-5.9 to -12.9 MW), southward flux reached 28--31 MW, and northward flux
was about 2.2 MW. The west progressive-wave index was negative, so its
reflection proxy was undefined. The old summary is retained as
`invalid_unbalanced_boundary_reflection_summary.csv`; the large model output
was removed before the corrected pair was prepared.
