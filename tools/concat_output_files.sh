#! /bin/sh

module load intel/2018.1 nco

#ROOT_WORK_DIR=$(readlink -f $PWD/..)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
ROOT_WORK_DIR=$(readlink -f $SCRIPT_DIR/..)
echo ROOT_WORK_DIR: $ROOT_WORK_DIR

POEM_WORK_DIR=$ROOT_WORK_DIR/POEM
PISM_WORK_DIR=$ROOT_WORK_DIR/PISM
SIM_START_YEAR=012821
SIM_END_YEAR=012911
CPL_TIMESTEP=10

set -x


# concatenate MOM output
cd $POEM_WORK_DIR/history
INPUT_FILES=$(echo `seq -f "%06g0101.ice-monthly.nc" \
    $SIM_START_YEAR $CPL_TIMESTEP $SIM_END_YEAR`)
OUTPUT_FILE=$(echo `printf "%06g-%06g.ice-monthly.nc" \
    $SIM_START_YEAR $SIM_END_YEAR`)
ncrcat --overwrite $INPUT_FILES $OUTPUT_FILE

INPUT_FILES=$(echo `seq -f "%06g0101.ice-yearly.nc" \
    $SIM_START_YEAR $CPL_TIMESTEP $SIM_END_YEAR`)
OUTPUT_FILE=$(echo `printf "%06g-%06g.ice-yearly.nc" \
    $SIM_START_YEAR $SIM_END_YEAR`)
ncrcat --overwrite $INPUT_FILES $OUTPUT_FILE

INPUT_FILES=$(echo `seq -f "%06g0101.ice-decadal.nc" \
    $SIM_START_YEAR $CPL_TIMESTEP $SIM_END_YEAR`)
OUTPUT_FILE=$(echo `printf "%06g-%06g.ice-decadal.nc" \
    $SIM_START_YEAR $SIM_END_YEAR`)
ncrcat --overwrite $INPUT_FILES $OUTPUT_FILE

INPUT_FILES=$(echo `seq -f "%06g0101.ocean-scalar.nc" \
    $SIM_START_YEAR $CPL_TIMESTEP $SIM_END_YEAR`)
OUTPUT_FILE=$(echo `printf "%06g-%06g.ocean-scalar.nc" \
    $SIM_START_YEAR $SIM_END_YEAR`)
ncrcat --overwrite $INPUT_FILES $OUTPUT_FILE

INPUT_FILES=$(echo `seq -f "%06g0101.ocean-monthly.nc" \
    $SIM_START_YEAR $CPL_TIMESTEP $SIM_END_YEAR`)
OUTPUT_FILE=$(echo `printf "%06g-%06g.ocean-monthly.nc" \
    $SIM_START_YEAR $SIM_END_YEAR`)
ncrcat --overwrite $INPUT_FILES $OUTPUT_FILE

INPUT_FILES=$(echo `seq -f "%06g0101.ocean-yearly.nc" \
    $SIM_START_YEAR $CPL_TIMESTEP $SIM_END_YEAR`)
OUTPUT_FILE=$(echo `printf "%06g-%06g.ocean-yearly.nc" \
    $SIM_START_YEAR $SIM_END_YEAR`)
ncrcat --overwrite $INPUT_FILES $OUTPUT_FILE

INPUT_FILES=$(echo `seq -f "%06g0101.ocean-decadal.nc" \
    $SIM_START_YEAR $CPL_TIMESTEP $SIM_END_YEAR`)
OUTPUT_FILE=$(echo `printf "%06g-%06g.ocean-decadal.nc" \
    $SIM_START_YEAR $SIM_END_YEAR`)
ncrcat --overwrite $INPUT_FILES $OUTPUT_FILE
cd $ROOT_WORK_DIR

# concatenate PISM output
cd $PISM_WORK_DIR/results
ncrcat --overwrite $(echo `seq -f "%06g.pism_extra.nc" $SIM_START_YEAR $CPL_TIMESTEP $SIM_END_YEAR`) \
    $SIM_START_YEAR-$SIM_END_YEAR.pism_extra.nc
ncrcat --overwrite $(echo `seq -f "%06g.pism_snap.nc" $SIM_START_YEAR $CPL_TIMESTEP $SIM_END_YEAR`) \
    $SIM_START_YEAR-$SIM_END_YEAR.pism_snap.nc
ncrcat --overwrite $(echo `seq -f "%06g.pism_ts.nc" $SIM_START_YEAR $CPL_TIMESTEP $SIM_END_YEAR`) \
    $SIM_START_YEAR-$SIM_END_YEAR.pism_ts.nc
cd $ROOT_WORK_DIR

# concatenate inter-model-processing files
cd $ROOT_WORK_DIR/x_MOM-to-PISM
ncrcat --overwrite $(echo `seq -f "%06g0101.processed_MOM.nc" $SIM_START_YEAR $CPL_TIMESTEP $SIM_END_YEAR`) \
    $SIM_START_YEAR-$SIM_END_YEAR.processed_MOM.nc
ncrcat --overwrite $(echo `seq -f "%06g0101.processed_MOM.anomaly.nc" $SIM_START_YEAR $CPL_TIMESTEP $SIM_END_YEAR`) \
    $SIM_START_YEAR-$SIM_END_YEAR.processed_MOM.anomaly.nc
ncrcat --overwrite $(echo `seq -f "%06g0101.PISM_input.nc" $SIM_START_YEAR $CPL_TIMESTEP $SIM_END_YEAR`) \
    $SIM_START_YEAR-$SIM_END_YEAR.PISM_input.nc
cd $ROOT_WORK_DIR/x_PISM-to-MOM
#ncrcat --overwrite $(echo `seq -f "%06g.basin_shelf_depth.nc" $SIM_START_YEAR $CPL_TIMESTEP $SIM_END_YEAR`) \
#    $SIM_START_YEAR-$SIM_END_YEAR.basin_shelf_depth.nc
ncrcat --overwrite $(echo `seq -f "%06g.fluxes.nc" $SIM_START_YEAR $CPL_TIMESTEP $SIM_END_YEAR`) \
    $SIM_START_YEAR-$SIM_END_YEAR.fluxes.nc
cd $ROOT_WORK_DIR

