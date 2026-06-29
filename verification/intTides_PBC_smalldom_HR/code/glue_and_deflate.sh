
#!/bin/bash
# Activate envirnment before running this script!
# conda activate mitgcm-env-py311

# Exit if any command fails
set -e

# Path to this script and the Python script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_NAME="make_netcdf.py"

# Input and output NetCDF file used for both Python generation and deflation
INPUT_NETCDF="$SCRIPT_DIR/../runTSini_M2K1/dynVarsGlob.nc"
OUTPUT_NETCDF="$INPUT_NETCDF"

# Iteration range for dynVars extraction
INI_ITER=621
END_ITER=96000
STEP=621
DT=9 # time step in seconds

# Run the Python script and pass the target NetCDF path and iter settings
echo "Running Python script: $SCRIPT_NAME"
python "$SCRIPT_DIR/$SCRIPT_NAME" "$INPUT_NETCDF" "$INI_ITER" "$END_ITER" "$STEP" "$DT"

# Deflate NetCDF file
echo "Deflating NetCDF file: $INPUT_NETCDF"
ncks -O -4 -L 4 --cnk_dmn time,1 "$INPUT_NETCDF" "$OUTPUT_NETCDF"
echo "Done. Compressed NetCDF saved as: $OUTPUT_NETCDF"

# Deflate NetCDF file
echo "Deflating NetCDF file: gridGlob.nc"
ncks -O -4 -L 4 --cnk_dmn time,1 "$SCRIPT_DIR/../runTSini_M2K1/gridGlob.nc" "$SCRIPT_DIR/../runTSini_M2K1/gridGlob.nc"
echo "Done. Compressed NetCDF saved as: $SCRIPT_DIR/../runTSini_M2K1/gridGlob.nc"

# Deflate NetCDF file
echo "Deflating NetCDF file: hydroVarsGlob.nc"
ncks -O -4 -L 4 --cnk_dmn time,1 "$SCRIPT_DIR/../runTSini_M2K1/hydroVarsGlob.nc" "$SCRIPT_DIR/../runTSini_M2K1/hydroVarsGlob.nc"
echo "Done. Compressed NetCDF saved as: $SCRIPT_DIR/../runTSini_M2K1/hydroVarsGlob.nc"

# Deflate NetCDF file
echo "Deflating NetCDF file: etaGlob.nc"
ncks -O -4 -L 4 --cnk_dmn time,1 "$SCRIPT_DIR/../runTSini_M2K1/etaGlob.nc" "$SCRIPT_DIR/../runTSini_M2K1/etaGlob.nc"
echo "Done. Compressed NetCDF saved as: $SCRIPT_DIR/../runTSini_M2K1/etaGlob.nc"

# Deflate NetCDF file
echo "Deflating NetCDF file: phihydGlob.nc"
ncks -O -4 -L 4 --cnk_dmn time,1 "$SCRIPT_DIR/../runTSini_M2K1/phihydGlob.nc" "$SCRIPT_DIR/../runTSini_M2K1/phihydGlob.nc"
echo "Done. Compressed NetCDF saved as: $SCRIPT_DIR/../runTSini_M2K1/phihydGlob.nc"

# Deflate NetCDF file
echo "Deflating NetCDF file: phibotGlob.nc"
ncks -O -4 -L 4 --cnk_dmn time,1 "$SCRIPT_DIR/../runTSini_M2K1/phibotGlob.nc" "$SCRIPT_DIR/../runTSini_M2K1/phibotGlob.nc"
echo "Done. Compressed NetCDF saved as: $SCRIPT_DIR/../runTSini_M2K1/phibotGlob.nc"
