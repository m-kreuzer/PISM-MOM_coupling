#! /bin/bash

set +x
module load intel/2018.1 nco

#ROOT_WORK_DIR=$(readlink -f $PWD/..)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
ROOT_WORK_DIR=$(readlink -f $SCRIPT_DIR/..)
echo ROOT_WORK_DIR: $ROOT_WORK_DIR

POEM_WORK_DIR=$ROOT_WORK_DIR/POEM
PISM_WORK_DIR=$ROOT_WORK_DIR/PISM
#SIM_EXT_DIR=$ROOT_WORK_DIR/../MOM5_PISM_16km_gmd-2020-230_run21_ext/
SIM_EXT_DIR=$ROOT_WORK_DIR/../MOM5_PISM_16km_gmd-2020-230_run21_ext2/
## ext
#SIM_START_YEAR=013100
#SIM_END_YEAR=014350
#SIM_EXT_START_YEAR=014360
#SIM_EXT_END_YEAR=015610
# ext2
SIM_START_YEAR=013100
SIM_END_YEAR=015610
SIM_EXT_START_YEAR=015620
SIM_EXT_END_YEAR=16420
#CPL_TIMESTEP=10

set -x

OCN_FILES="ice-monthly.nc ice-yearly.nc ice-decadal.nc ocean-scalar.nc ocean-monthly.nc ocean-yearly.nc ocean-decadal.nc"
ICE_FILES="pism_extra.nc pism_snap.nc pism_ts.nc"
OTI_FILES="processed_MOM.nc processed_MOM.anomaly.nc PISM_input.nc"
ITO_FILES="basin_shelf_depth.nc fluxes.nc"

## concatenate MOM output
cd $POEM_WORK_DIR/history
for OF in $OCN_FILES; do
    INPUT_FILE1=$(echo `printf "%06g-%06g.$OF" $SIM_START_YEAR $SIM_END_YEAR`)
    INPUT_FILE2=$(echo $SIM_EXT_DIR/POEM/history/`printf "%06g-%06g.$OF" \
        $SIM_EXT_START_YEAR $SIM_EXT_END_YEAR`)
    OUTPUT_FILE=$(echo `printf "%06g-%06g.$OF" \
        $SIM_START_YEAR $SIM_EXT_END_YEAR`)
    ncrcat --overwrite $INPUT_FILE1 $INPUT_FILE2 $OUTPUT_FILE
done
cd $ROOT_WORK_DIR

# concatenate PISM output
cd $PISM_WORK_DIR/results
for IF in $ICE_FILES; do
    INPUT_FILE1=$(echo `printf "%06g-%06g.$IF" $SIM_START_YEAR $SIM_END_YEAR`)
    INPUT_FILE2=$(echo $SIM_EXT_DIR/PISM/results/`printf "%06g-%06g.$IF" \
        $SIM_EXT_START_YEAR $SIM_EXT_END_YEAR`)
    OUTPUT_FILE=$(echo `printf "%06g-%06g.$IF" \
        $SIM_START_YEAR $SIM_EXT_END_YEAR`)
    ncrcat --overwrite $INPUT_FILE1 $INPUT_FILE2 $OUTPUT_FILE
done
cd $ROOT_WORK_DIR

# concatenate inter-model-processing files
cd $ROOT_WORK_DIR/x_MOM-to-PISM
for F in $OTI_FILES; do
    INPUT_FILE1=$(echo `printf "%06g-%06g.$F" $SIM_START_YEAR $SIM_END_YEAR`)
    INPUT_FILE2=$(echo $SIM_EXT_DIR/x_MOM-to-PISM/`printf "%06g-%06g.$F" \
        $SIM_EXT_START_YEAR $SIM_EXT_END_YEAR`)
    OUTPUT_FILE=$(echo `printf "%06g-%06g.$F" \
        $SIM_START_YEAR $SIM_EXT_END_YEAR`)
    ncrcat --overwrite $INPUT_FILE1 $INPUT_FILE2 $OUTPUT_FILE
done

cd $ROOT_WORK_DIR/x_PISM-to-MOM
for F in $ITO_FILES; do
    INPUT_FILE1=$(echo `printf "%06g-%06g.$F" $SIM_START_YEAR $SIM_END_YEAR`)
    INPUT_FILE2=$(echo $SIM_EXT_DIR/x_PISM-to-MOM/`printf "%06g-%06g.$F" \
        $SIM_EXT_START_YEAR $SIM_EXT_END_YEAR`)
    OUTPUT_FILE=$(echo `printf "%06g-%06g.$F" \
        $SIM_START_YEAR $SIM_EXT_END_YEAR`)
    ncrcat --overwrite $INPUT_FILE1 $INPUT_FILE2 $OUTPUT_FILE
done
cd $ROOT_WORK_DIR
