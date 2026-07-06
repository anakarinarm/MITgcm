#!/usr/bin/env python3
"""Compute a 2-D barotropic-to-baroclinic conversion map.

This script is written for the production diagnostics in this run directory:

    output/diag_state_3d.*  with UVELMASS, VVELMASS, PHIHYD
    output/diag_state_2d.*  with PHIBOT

The diagnostic computed here is a common topographic-conversion estimate:

    C = -rho0 * < phi_bc_bot' * (U_bt . grad(H))' >

where

    phi_bc_bot = PHIBOT - vertical_mean(PHIHYD)

is a bottom baroclinic pressure-potential proxy, H is positive ocean depth,
U_bt is the depth-mean horizontal velocity interpolated to tracer-cell
centers, and primes denote removal of the analysis-period mean.

With this convention, positive C means barotropic-tide energy is converted
into baroclinic/internal-tide energy.  Units are W/m^2.  Domain-integrated
conversion is computed by integrating over wet-cell area.

Notes:

- This is a first-pass, model-diagnostic map.  It is intentionally explicit
  about the sign convention because conversion-sign conventions differ across
  papers and codes.
- Earlier versions of this script used the opposite sign.  Known generation
  sites appeared as negative there; the default is now flipped so known
  generators are positive.  Use --legacy-sign to reproduce the old convention.
- PHIBOT/PHIHYD are pressure potentials (p/rho0).  Multiplying by rho0 gives
  pressure units in the final conversion.
- The calculation streams one time slice at a time from MDS output, so it does
  not load the whole production run into memory.
"""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


@dataclass
class Meta:
    nx: int
    ny: int
    nz: int
    nrecords: int
    dtype: np.dtype
    fields: list[str]
    time: float | None


def parse_meta(path: Path) -> Meta:
    text = path.read_text()
    dim_match = re.search(r"dimList\s*=\s*\[(.*?)\];", text, re.S)
    if not dim_match:
        raise ValueError(f"Cannot parse dimList from {path}")
    numbers = [int(x) for x in re.findall(r"[-+]?\d+", dim_match.group(1))]
    dims = numbers[0::3]
    if len(dims) == 2:
        nx, ny = dims
        nz = 1
    elif len(dims) == 3:
        nx, ny, nz = dims
    else:
        raise ValueError(f"Unexpected dimensions {dims} in {path}")

    precision = re.search(r"dataprec\s*=\s*\[\s*'([^']+)'", text).group(1)
    dtype = np.dtype(">f4" if precision.strip() == "float32" else ">f8")

    nrec_match = re.search(r"nrecords\s*=\s*\[\s*(\d+)", text)
    nrecords = int(nrec_match.group(1)) if nrec_match else 1

    field_match = re.search(r"fldList\s*=\s*\{(.*?)\};", text, re.S)
    fields = re.findall(r"'([^']+)'", field_match.group(1)) if field_match else []
    fields = [field.strip() for field in fields]

    time_match = re.search(r"timeInterval\s*=\s*\[\s*([^\]]+)", text)
    time = float(time_match.group(1)) if time_match else None
    return Meta(nx, ny, nz, nrecords, dtype, fields, time)


def memmap_mds(meta_path: Path) -> tuple[np.memmap, Meta]:
    meta = parse_meta(meta_path)
    data_path = meta_path.with_suffix(".data")
    expected = meta.nrecords * meta.nz * meta.ny * meta.nx
    actual = data_path.stat().st_size // meta.dtype.itemsize
    if actual != expected:
        raise ValueError(f"{data_path}: expected {expected} values, got {actual}")
    array = np.memmap(
        data_path,
        dtype=meta.dtype,
        mode="r",
        shape=(meta.nrecords, meta.nz, meta.ny, meta.nx),
    )
    return array, meta


def read_grid(run_dir: Path, name: str) -> np.ndarray:
    array, meta = memmap_mds(run_dir / f"{name}.meta")
    if meta.nrecords != 1:
        raise ValueError(f"{name}: expected one record, got {meta.nrecords}")
    result = np.asarray(array[0], dtype=np.float64)
    return result if meta.nz > 1 else result[0]


def central_gradient_metric(field: np.ndarray, dx: np.ndarray, dy: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return d(field)/dx, d(field)/dy on C centers using local metric distances."""
    dfdx = np.full_like(field, np.nan, dtype=np.float64)
    dfdy = np.full_like(field, np.nan, dtype=np.float64)

    # Distance from center i-1 to center i+1 is the sum of adjacent half-cell
    # distances.  This works for the telescopic grid without assuming uniform dx.
    denom_x = dx[:, 1:-1] + 0.5 * (dx[:, :-2] + dx[:, 2:])
    denom_y = dy[1:-1, :] + 0.5 * (dy[:-2, :] + dy[2:, :])
    dfdx[:, 1:-1] = (field[:, 2:] - field[:, :-2]) / denom_x
    dfdy[1:-1, :] = (field[2:, :] - field[:-2, :]) / denom_y

    # One-sided edges.
    dfdx[:, 0] = (field[:, 1] - field[:, 0]) / (0.5 * (dx[:, 0] + dx[:, 1]))
    dfdx[:, -1] = (field[:, -1] - field[:, -2]) / (0.5 * (dx[:, -1] + dx[:, -2]))
    dfdy[0, :] = (field[1, :] - field[0, :]) / (0.5 * (dy[0, :] + dy[1, :]))
    dfdy[-1, :] = (field[-1, :] - field[-2, :]) / (0.5 * (dy[-1, :] + dy[-2, :]))
    return dfdx, dfdy


def depth_mean_velocity_to_c(
    u_mass: np.ndarray,
    v_mass: np.ndarray,
    drf: np.ndarray,
    hfac_w: np.ndarray,
    hfac_s: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Depth-average U/V mass velocity and interpolate from faces to C centers."""
    zw = hfac_w * drf[:, None, None]
    zs = hfac_s * drf[:, None, None]
    hw = np.sum(zw, axis=0)
    hs = np.sum(zs, axis=0)

    u_face = np.divide(
        np.sum(u_mass * zw, axis=0),
        hw,
        out=np.zeros_like(hw, dtype=np.float64),
        where=hw > 0,
    )
    v_face = np.divide(
        np.sum(v_mass * zs, axis=0),
        hs,
        out=np.zeros_like(hs, dtype=np.float64),
        where=hs > 0,
    )

    # U is on west/east faces; V is on south/north faces.  Average neighboring
    # faces to tracer-cell centers.  Use nearest available face at the domain edge.
    u_c = np.empty_like(u_face)
    v_c = np.empty_like(v_face)
    u_c[:, :-1] = 0.5 * (u_face[:, :-1] + u_face[:, 1:])
    u_c[:, -1] = u_face[:, -1]
    v_c[:-1, :] = 0.5 * (v_face[:-1, :] + v_face[1:, :])
    v_c[-1, :] = v_face[-1, :]
    return u_c, v_c


def vertical_mean_phi(phi: np.ndarray, drf: np.ndarray, hfac_c: np.ndarray) -> np.ndarray:
    weights = hfac_c * drf[:, None, None]
    depth = np.sum(weights, axis=0)
    return np.divide(
        np.sum(phi * weights, axis=0),
        depth,
        out=np.zeros_like(depth, dtype=np.float64),
        where=depth > 0,
    )


def finite_or_nan(array: np.ndarray, wet: np.ndarray) -> np.ndarray:
    out = np.asarray(array, dtype=np.float64).copy()
    out[~wet] = np.nan
    out[~np.isfinite(out)] = np.nan
    return out


def write_netcdf(
    path: Path,
    xc: np.ndarray,
    yc: np.ndarray,
    wet: np.ndarray,
    physical: np.ndarray,
    conversion: np.ndarray,
    phi_bc_rms: np.ndarray,
    wbt_rms: np.ndarray,
    attrs: dict[str, float | str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with Dataset(path, "w") as ds:
        ny, nx = conversion.shape
        ds.createDimension("y", ny)
        ds.createDimension("x", nx)

        xvar = ds.createVariable("XC", "f8", ("y", "x"), zlib=True)
        yvar = ds.createVariable("YC", "f8", ("y", "x"), zlib=True)
        mask = ds.createVariable("wet_mask", "i1", ("y", "x"), zlib=True)
        pmask = ds.createVariable("physical_domain_mask", "i1", ("y", "x"), zlib=True)
        conv = ds.createVariable(
            "conversion_bt_to_bc", "f8", ("y", "x"), zlib=True, complevel=4, fill_value=np.nan
        )
        pos = ds.createVariable(
            "conversion_positive_part", "f8", ("y", "x"), zlib=True, complevel=4, fill_value=np.nan
        )
        neg = ds.createVariable(
            "conversion_negative_part", "f8", ("y", "x"), zlib=True, complevel=4, fill_value=np.nan
        )
        prms = ds.createVariable(
            "phi_bc_bottom_rms", "f8", ("y", "x"), zlib=True, complevel=4, fill_value=np.nan
        )
        wrms = ds.createVariable(
            "w_bt_topographic_rms", "f8", ("y", "x"), zlib=True, complevel=4, fill_value=np.nan
        )

        xvar[:] = xc
        yvar[:] = yc
        mask[:] = wet.astype(np.int8)
        pmask[:] = physical.astype(np.int8)
        conv[:] = conversion
        pos[:] = np.where(conversion > 0, conversion, 0.0)
        neg[:] = np.where(conversion < 0, conversion, 0.0)
        prms[:] = phi_bc_rms
        wrms[:] = wbt_rms

        xvar.units = "degrees_east"
        yvar.units = "degrees_north"
        conv.units = "W m-2"
        conv.long_name = "barotropic-to-baroclinic conversion"
        conv.positive = "barotropic tide to baroclinic/internal tide"
        pos.units = "W m-2"
        neg.units = "W m-2"
        prms.units = "m2 s-2"
        wrms.units = "m s-1"
        ds.title = "Barotropic-to-baroclinic conversion map"
        for key, value in attrs.items():
            setattr(ds, key, value)


def write_png(path: Path, xc: np.ndarray, yc: np.ndarray, conversion: np.ndarray) -> None:
    # Keep matplotlib cache out of the user's home directory on systems where
    # ~/.matplotlib is not writable.
    import os

    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    path.parent.mkdir(parents=True, exist_ok=True)
    data = conversion * 1e3  # mW/m^2
    vmax = np.nanpercentile(np.abs(data), 99.0)
    if not np.isfinite(vmax) or vmax == 0:
        vmax = 1.0
    fig, ax = plt.subplots(figsize=(8, 7), constrained_layout=True)
    pcm = ax.pcolormesh(xc, yc, data, shading="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    cb = fig.colorbar(pcm, ax=ax, shrink=0.85)
    cb.set_label("BT → BC conversion (mW m$^{-2}$)")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title("Barotropic-to-baroclinic conversion, analysis-period mean")
    ax.set_aspect("equal", adjustable="box")
    fig.savefig(path, dpi=200)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--output-dir", type=Path, default=None,
                        help="Output directory. Default: RUN_DIR/output/postprocessing")
    parser.add_argument("--rho0", type=float, default=1026.0)
    parser.add_argument("--start-time", type=float, default=None,
                        help="Optional minimum diagnostic time in seconds.")
    parser.add_argument("--end-time", type=float, default=None,
                        help="Optional maximum diagnostic time in seconds.")
    parser.add_argument("--prefix-3d", default="diag_state_3d")
    parser.add_argument("--prefix-2d", default="diag_state_2d")
    parser.add_argument("--buffer-cells", type=int, default=60,
                        help="Cells excluded from west/south/north buffers for physical-domain integrals.")
    parser.add_argument("--legacy-sign", action="store_true",
                        help="Use the old sign: +rho0*covariance instead of -rho0*covariance.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_dir = args.run_dir.resolve()
    diag_dir = run_dir / "output"
    out_dir = args.output_dir or (diag_dir / "postprocessing")
    out_dir.mkdir(parents=True, exist_ok=True)

    files_3d = sorted(diag_dir.glob(f"{args.prefix_3d}.[0-9]*.meta"))
    files_2d = sorted(diag_dir.glob(f"{args.prefix_2d}.[0-9]*.meta"))
    if not files_3d or not files_2d:
        raise FileNotFoundError(f"No {args.prefix_3d}/{args.prefix_2d} MDS diagnostics found in {diag_dir}")
    if len(files_3d) != len(files_2d):
        raise ValueError("3-D and 2-D diagnostic file counts differ")

    drf = read_grid(run_dir, "DRF").reshape(-1)
    hfac_c = read_grid(run_dir, "hFacC")
    hfac_w = read_grid(run_dir, "hFacW")
    hfac_s = read_grid(run_dir, "hFacS")
    depth = read_grid(run_dir, "Depth")
    dxc = read_grid(run_dir, "DXC")
    dyc = read_grid(run_dir, "DYC")
    rac = read_grid(run_dir, "RAC")
    xc = read_grid(run_dir, "XC")
    yc = read_grid(run_dir, "YC")

    wet = (hfac_c[0] > 0) & (depth > 0)
    jj, ii = np.indices(depth.shape)
    # The HR telescopic setup preserves the original 600x600 physical domain at
    # i=61:660, j=61:660 in MITgcm 1-based indexing.  In Python 0-based
    # indexing this is i>=60 and 60<=j<660.  There is no east buffer.
    physical = (
        wet
        & (ii >= args.buffer_cells)
        & (jj >= args.buffer_cells)
        & (jj < depth.shape[0] - args.buffer_cells)
    )
    dHdx, dHdy = central_gradient_metric(depth, dxc, dyc)
    dHdx[~wet] = np.nan
    dHdy[~wet] = np.nan

    n = 0
    sum_phi = np.zeros_like(depth, dtype=np.float64)
    sum_w = np.zeros_like(depth, dtype=np.float64)
    sum_phi_w = np.zeros_like(depth, dtype=np.float64)
    sum_phi2 = np.zeros_like(depth, dtype=np.float64)
    sum_w2 = np.zeros_like(depth, dtype=np.float64)
    times: list[float] = []

    for file_3d, file_2d in zip(files_3d, files_2d):
        array_3d, meta_3d = memmap_mds(file_3d)
        array_2d, meta_2d = memmap_mds(file_2d)
        if meta_3d.time is None or meta_2d.time is None:
            raise ValueError("Diagnostic metadata missing timeInterval")
        if abs(meta_3d.time - meta_2d.time) > 0.1:
            raise ValueError(f"3-D and 2-D diagnostics are not synchronous: {file_3d}, {file_2d}")
        if args.start_time is not None and meta_3d.time < args.start_time:
            continue
        if args.end_time is not None and meta_3d.time > args.end_time:
            continue

        fields_3d = meta_3d.fields
        fields_2d = meta_2d.fields
        for needed in ("UVELMASS", "VVELMASS", "PHIHYD"):
            if needed not in fields_3d:
                raise ValueError(f"{needed} missing from {file_3d}")
        if "PHIBOT" not in fields_2d:
            raise ValueError(f"PHIBOT missing from {file_2d}")

        u_mass = np.asarray(array_3d[fields_3d.index("UVELMASS")], dtype=np.float64)
        v_mass = np.asarray(array_3d[fields_3d.index("VVELMASS")], dtype=np.float64)
        phi = np.asarray(array_3d[fields_3d.index("PHIHYD")], dtype=np.float64)
        phibot = np.asarray(array_2d[fields_2d.index("PHIBOT"), 0], dtype=np.float64)

        u_bt, v_bt = depth_mean_velocity_to_c(u_mass, v_mass, drf, hfac_w, hfac_s)
        # Topographic forcing proxy.  With positive ocean depth H, Ubt.grad(H)
        # has the opposite sign of the physical bottom vertical velocity in
        # z-up coordinates.  We apply the sign convention below when converting
        # covariance to BT->BC conversion.
        w_topo = u_bt * dHdx + v_bt * dHdy
        phi_bc_bot = phibot - vertical_mean_phi(phi, drf, hfac_c)

        w_topo = finite_or_nan(w_topo, wet)
        phi_bc_bot = finite_or_nan(phi_bc_bot, wet)

        valid = wet & np.isfinite(w_topo) & np.isfinite(phi_bc_bot)
        sum_phi[valid] += phi_bc_bot[valid]
        sum_w[valid] += w_topo[valid]
        sum_phi_w[valid] += phi_bc_bot[valid] * w_topo[valid]
        sum_phi2[valid] += phi_bc_bot[valid] ** 2
        sum_w2[valid] += w_topo[valid] ** 2
        times.append(meta_3d.time)
        n += 1
        print(f"processed {n:03d}: t={meta_3d.time:.1f} s {file_3d.name}", flush=True)

    if n == 0:
        raise ValueError("No diagnostics selected by the requested time range")

    mean_phi = sum_phi / n
    mean_w = sum_w / n
    covariance = sum_phi_w / n - mean_phi * mean_w
    sign = 1.0 if args.legacy_sign else -1.0
    conversion = sign * args.rho0 * covariance
    phi_rms = np.sqrt(np.maximum(sum_phi2 / n - mean_phi**2, 0.0))
    w_rms = np.sqrt(np.maximum(sum_w2 / n - mean_w**2, 0.0))

    conversion = finite_or_nan(conversion, wet)
    phi_rms = finite_or_nan(phi_rms, wet)
    w_rms = finite_or_nan(w_rms, wet)
    area_full = np.where(wet, rac, 0.0)
    area_phys = np.where(physical, rac, 0.0)

    total_w = float(np.nansum(conversion * area_full))
    positive_w = float(np.nansum(np.where(conversion > 0, conversion, 0.0) * area_full))
    negative_w = float(np.nansum(np.where(conversion < 0, conversion, 0.0) * area_full))
    total_phys_w = float(np.nansum(conversion * area_phys))
    positive_phys_w = float(np.nansum(np.where(conversion > 0, conversion, 0.0) * area_phys))
    negative_phys_w = float(np.nansum(np.where(conversion < 0, conversion, 0.0) * area_phys))

    nc_path = out_dir / "barotropic_to_baroclinic_conversion.nc"
    png_path = out_dir / "barotropic_to_baroclinic_conversion.png"
    csv_path = out_dir / "barotropic_to_baroclinic_conversion_integrals.csv"

    attrs = {
        "rho0_kg_m3": args.rho0,
        "n_snapshots": n,
        "time_start_s": min(times),
        "time_end_s": max(times),
        "formula": (
            "conversion = "
            + ("+" if args.legacy_sign else "-")
            + "rho0 * mean((PHIBOT - vertical_mean(PHIHYD))' * (Ubt*dHdx + Vbt*dHdy)')"
        ),
        "positive_sign": "positive means barotropic tide to baroclinic/internal tide",
        "legacy_sign": int(args.legacy_sign),
        "total_conversion_W": total_w,
        "positive_conversion_W": positive_w,
        "negative_conversion_W": negative_w,
        "physical_domain_total_conversion_W": total_phys_w,
        "physical_domain_positive_conversion_W": positive_phys_w,
        "physical_domain_negative_conversion_W": negative_phys_w,
        "physical_domain_definition": f"wet cells with i>={args.buffer_cells}, j>={args.buffer_cells}, j<ny-{args.buffer_cells}",
        "source_run_dir": str(run_dir),
    }
    write_netcdf(nc_path, xc, yc, wet, physical, conversion, phi_rms, w_rms, attrs)
    write_png(png_path, xc, yc, conversion)
    with csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["quantity", "value_W"])
        writer.writerow(["full_domain_total_conversion", f"{total_w:.16e}"])
        writer.writerow(["full_domain_positive_conversion", f"{positive_w:.16e}"])
        writer.writerow(["full_domain_negative_conversion", f"{negative_w:.16e}"])
        writer.writerow(["physical_domain_total_conversion", f"{total_phys_w:.16e}"])
        writer.writerow(["physical_domain_positive_conversion", f"{positive_phys_w:.16e}"])
        writer.writerow(["physical_domain_negative_conversion", f"{negative_phys_w:.16e}"])

    print("\nWrote:")
    print(f"  {nc_path}")
    print(f"  {png_path}")
    print(f"  {csv_path}")
    print("\nIntegrated conversion:")
    print("  full computational domain:")
    print(f"    total    = {total_w:.6e} W")
    print(f"    positive = {positive_w:.6e} W")
    print(f"    negative = {negative_w:.6e} W")
    print("  physical domain, excluding west/south/north 60-cell buffers:")
    print(f"    total    = {total_phys_w:.6e} W")
    print(f"    positive = {positive_phys_w:.6e} W")
    print(f"    negative = {negative_phys_w:.6e} W")


if __name__ == "__main__":
    main()
