# Read mds files from MITgcm (binary output) using xmitgcm and save as NetCDF files.
import sys
from pathlib import Path
import xarray as xr
import xmitgcm as xgcm

if len(sys.argv) not in (2, 6):
    print(f"Usage: {sys.argv[0]} <output_netcdf_path> [ini_iter end_iter step dt]")
    sys.exit(1)

output_netcdf = Path(sys.argv[1])
data_root = output_netcdf.parent

if not data_root.exists():
    raise FileNotFoundError(f"Data directory not found: {data_root}")

if len(sys.argv) == 6:
    ini_iter = int(sys.argv[2])
    end_iter = int(sys.argv[3])
    step = int(sys.argv[4])
    dt = int(sys.argv[5])
else:
    raise ValueError("Iteration range (ini_iter, end_iter, step) and time step (dt) must be provided as command-line arguments.")

iters = list(range(ini_iter, end_iter + 1, step))

print(f"Creating NetCDF output in: {output_netcdf}")
print(f"Using MITgcm dataset root: {data_root}")
print(f"Using iteration range: {iters}")

vars = ['dynVars']
ds = xgcm.open_mdsdataset(str(data_root), iters=iters, delta_t=dt, geometry='sphericalpolar', prefix=vars)
ds.to_netcdf(path=str(output_netcdf), mode='w', format='NETCDF4')

hydro_filepath = data_root / 'hydroVarsGlob.nc'
vars = ['hydroVars']
ds = xgcm.open_mdsdataset(str(data_root), iters=iters, delta_t=dt, geometry='sphericalpolar', prefix=vars)
ds.to_netcdf(path=str(hydro_filepath), mode='w', format='NETCDF4')

eta_filepath = data_root / 'etaGlob.nc'
vars = ['etan']
ds = xgcm.open_mdsdataset(str(data_root), iters=iters, delta_t=dt, geometry='sphericalpolar', prefix=vars)
ds.to_netcdf(path=str(eta_filepath), mode='w', format='NETCDF4')

phihyd_filepath = data_root / 'phihydGlob.nc'
vars = ['phihyd']
ds = xgcm.open_mdsdataset(str(data_root), iters=iters, delta_t=dt, geometry='sphericalpolar', prefix=vars)
ds.to_netcdf(path=str(phihyd_filepath), mode='w', format='NETCDF4')

phibot_filepath = data_root / 'phibotGlob.nc'
vars = ['phibot']
ds = xgcm.open_mdsdataset(str(data_root), iters=iters, delta_t=dt, geometry='sphericalpolar', prefix=vars)
ds.to_netcdf(path=str(phibot_filepath), mode='w', format='NETCDF4')

vars2 = [
    'Depth', 'DRC', 'DRF', 'DXC', 'DXG', 'DXF', 'DXV', 'DYC', 'DYF', 'DYG', 'DYU',
    'hFacC', 'hFacW', 'hFacS', 'maskInC', 'maskInS', 'maskInW', 'RAC', 'RAS', 'RAW', 'RAZ',
    'RC', 'RF', 'XC', 'XG', 'YC', 'YG'
]
ds = xgcm.open_mdsdataset(str(data_root), delta_t=dt, geometry='sphericalpolar', prefix=vars2)
grid_filepath = data_root / 'gridGlob.nc'
ds.to_netcdf(path=str(grid_filepath), mode='w', format='NETCDF4')

print(f"Wrote dynVars output: {output_netcdf}")
print(f"Wrote grid output: {grid_filepath}")

