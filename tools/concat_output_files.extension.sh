#! /bin/bash

set +x
module load intel/2018.1 nco

#ROOT_WORK_DIR=$(readlink -f $PWD/..)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
ROOT_WORK_DIR=$(readlink -f $SCRIPT_DIR/..)
echo ROOT_WORK_DIR: $ROOT_WORK_DIR

POEM_WORK_DIR=$ROOT_WORK_DIR/POEM
PISM_WORK_DIR=$ROOT_WORK_DIR/PISM
SIM_EXT_DIR=$ROOT_WORK_DIR/../MOM5_PISM1.0hash_16km_1pctCO2ext_CCSM4_run01.RESTART_017625/

SIM_START_YEAR=017235
SIM_END_YEAR=017625
SIM_EXT_START_YEAR=017635
SIM_EXT_END_YEAR=021225

set -x

OCN_FILES="ice-monthly.nc ice-yearly.nc ice-decadal.nc ice-decadal_min.nc ice-decadal_max.nc ocean-scalar.nc ocean-monthly.nc ocean-yearly.nc ocean-decadal.nc ocean-decadal_min.nc ocean-decadal_max.nc"
ICE_FILES="pism_extra.nc pism_snap.nc pism_ts.nc"
OTI_FILES="processed_MOM.nc sealevel.processed_MOM.anomaly.nc sealevel.PISM_input.nc tracer.processed_MOM.anomaly.nc tracer.PISM_input.nc"
ITO_FILES="basal_melt_input_depth.nc pico_input_depth.nc fluxes.nc runoff_reference.nc"

## concatenate MOM output
cd $POEM_WORK_DIR/history
for OF in $OCN_FILES; do
    INPUT_FILE1=$(echo `printf "%06g-%06g.$OF" $SIM_START_YEAR $SIM_END_YEAR`)
    INPUT_FILE2=$(echo $SIM_EXT_DIR/POEM/history/`printf "%06g-%06g.$OF" \
        $SIM_EXT_START_YEAR $SIM_EXT_END_YEAR`)
    OUTPUT_FILE=$(echo `printf "%06g-%06g.$OF" \
        $SIM_START_YEAR $SIM_EXT_END_YEAR`)
    echo ncrcat --overwrite $INPUT_FILE1 $INPUT_FILE2 $OUTPUT_FILE
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
    echo ncrcat --overwrite $INPUT_FILE1 $INPUT_FILE2 $OUTPUT_FILE
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
    echo ncrcat --overwrite $INPUT_FILE1 $INPUT_FILE2 $OUTPUT_FILE
    ncrcat --overwrite $INPUT_FILE1 $INPUT_FILE2 $OUTPUT_FILE
done

cd $ROOT_WORK_DIR/x_PISM-to-MOM
for F in $ITO_FILES; do
    INPUT_FILE1=$(echo `printf "%06g-%06g.$F" $SIM_START_YEAR $SIM_END_YEAR`)
    INPUT_FILE2=$(echo $SIM_EXT_DIR/x_PISM-to-MOM/`printf "%06g-%06g.$F" \
        $SIM_EXT_START_YEAR $SIM_EXT_END_YEAR`)
    OUTPUT_FILE=$(echo `printf "%06g-%06g.$F" \
        $SIM_START_YEAR $SIM_EXT_END_YEAR`)
    echo ncrcat --overwrite $INPUT_FILE1 $INPUT_FILE2 $OUTPUT_FILE
    ncrcat --overwrite $INPUT_FILE1 $INPUT_FILE2 $OUTPUT_FILE
done
cd $ROOT_WORK_DIR
