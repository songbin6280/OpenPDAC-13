#!/bin/sh

# =============================================================================
# Script 4: Post-Processing
# -----------------------------------------------------------------------------
# This script runs a series of Python-based post-processing tasks to
# generate plots, maps, and other visualizations from the simulation results.
#
# It performs the following steps:
# 1. Activates the dedicated 'OpenPDACconda' environment.
# 2. Cleans up previous post-processing results.
# 3. Generates isosurface plots (plotIso.py).
# 4. Generates ballistic trajectory maps (plotBallistics.py).
# 5. Creates final raster maps from the flow data (createMaps.py).
# 6. Deactivates the Conda environment.
#
# Usage: ./04_run_postprocessing.sh
# =============================================================================

# Change to the script's directory for robust execution
cd "${0%/*}" || exit 1

# Source the OpenFOAM functions for running applications
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"

# =========================================================================
# Activate Conda environment for all Python scripts
# =========================================================================
echo "--> Activating Conda environment: OpenPDACconda"

# Source Conda's shell functions to make 'conda activate' available in scripts
# This is the most robust method and works regardless of installation path.
source "$(conda info --base)/etc/profile.d/conda.sh"

# Activate the specific environment
conda activate OpenPDACconda
# =========================================================================


# --- TASK 1: GENERATE BALLISTIC IMPACT RASTER MAPS ---

echo "--> Cleaning up previous ballistic raster maps..."
rm -rf postProcessing/raster_maps/ballistics

echo "--> Generating ballistic trajectory plots with plotBallistics.py..."
#python3 plotBallistics.py > log.plotBallistics


# --- TASK 2: CREATE RASTER MAPS OF THE FLOW ---

echo "--> Cleaning up previous flow raster maps..."
rm -rf postProcessing/raster_maps/flow

echo "--> Creating final raster maps with createMaps.py..."
# This script also appears to process data in parallel. Adjust --np as needed.
#python3 createMaps.py -np 10 > log.createMaps

# --- TASK 3: GENERATE ISOSURFACE PLOTS ---

echo "--> Cleaning up previous isosurface plot frames..."
rm -rf postProcessing/frames_png

echo "--> Generating isosurface plots with plotIso.py..."
# This script may process data in parallel. Ensure the number of processors
# (--np) matches a relevant setting for your case or script logic.
python3 plotIso.py --np 10 --fr 10 --ps 50 > log.plotiso


# =========================================================================
# Deactivate the Conda environment
# =========================================================================
conda deactivate
echo "--> Conda environment deactivated."
# =========================================================================


# -----------------------------------------------------------------------------
echo
echo "Post-processing script finished successfully."
# =============================================================================
