#! /bin/bash
echo " > script to delete intermediate output files"
echo


################################# USER SETTINGS ################################
SIM_START_YEAR=013100   # first output timestamp
SIM_END_YEAR=017090     # last output timestamp
CPL_TIMESTEP=10
KEEP_FREQ_YEAR=500  # don't remove restart files for every n years
################################################################################

### check if command line argument -s is given
#   -> simulates script functionality without removing files
SIMULATE=false

while getopts s flag
do
    case "${flag}" in
        s) SIMULATE=true ;;
        *) echo 'Error in command line parsing' >&2
           exit 1
    esac
done

if [ $SIMULATE = true ]; then
    echo "----------------------------------------------------------------"
    echo ">>> simulating script functionality -> no files are removed! <<<"
    echo "----------------------------------------------------------------"
    echo
    echo
fi

#module load intel/2018.1 nco

### set paths
#ROOT_WORK_DIR=$(readlink -f $PWD/..)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
ROOT_WORK_DIR=$(readlink -f $SCRIPT_DIR/..)
echo ROOT_WORK_DIR: $ROOT_WORK_DIR

POEM_WORK_DIR=$ROOT_WORK_DIR/POEM
PISM_WORK_DIR=$ROOT_WORK_DIR/PISM


#set -x

# -----------------------------------------------------------------------------
cd $PISM_WORK_DIR/results

echo
echo "-------------- PISM ------------"
echo " remove PISM restart files, only keep a few"
PISM_KEEP_FILES=$(echo `seq -f '%06g.pism_out.nc' \
    $(echo "$SIM_START_YEAR-$CPL_TIMESTEP" | bc -l) \
    $KEEP_FREQ_YEAR \
    $SIM_END_YEAR`)
echo files not to remove: $PISM_KEEP_FILES
if [ $SIMULATE = false ]; then
    mkdir keep
    mv $PISM_KEEP_FILES keep/
    rm *.pism_out.nc
    mv keep/* .
    rm -d keep
fi

echo " remove PISM extra and time series output"
PISM_TS_CAT=$(echo `printf "%06g-%06g.pism_ts.nc" \
    $SIM_START_YEAR $SIM_END_YEAR`)
if [ -f "$PISM_TS_CAT" ]; then
    echo " Concated file $PISM_TS_CAT exists. Deleting individual files now."
    if [ $SIMULATE = false ]; then
        find . -regextype posix-egrep -regex '\./[0-9]{6}+\.pism_ts\.nc' -exec rm {} \;
    fi
fi

PISM_EXTRA_CAT=$(echo `printf "%06g-%06g.pism_extra.nc" \
    $SIM_START_YEAR $SIM_END_YEAR`)
if [ -f "$PISM_EXTRA_CAT" ]; then
    echo " Concated file $PISM_EXTRA_CAT exists. Deleting individual files now." 
    if [ $SIMULATE = false ]; then
        find . -regextype posix-egrep -regex '\./[0-9]{6}+\.pism_extra\.nc' -exec rm {} \;
    fi
fi


# -----------------------------------------------------------------------------
cd $ROOT_WORK_DIR/POEM/history
echo
echo "--------- POEM --------"
echo " remove POEM intermediate restart files, only keep a few"
POEM_KEEP_FILES=$(echo `seq -f '%06g0101.RESTART.tar.bz2' \
    $(echo "$SIM_START_YEAR-$CPL_TIMESTEP" | bc -l) \
    $KEEP_FREQ_YEAR \
    $SIM_END_YEAR`)
echo files not to remove: $POEM_KEEP_FILES
if [ $SIMULATE = false ]; then
    mkdir keep
    mv $POEM_KEEP_FILES keep/
    rm *.RESTART.tar.bz2
    mv keep/* .
    rm -d keep
fi


# -----------------------------------------------------------------------------
cd $ROOT_WORK_DIR/x_MOM-to-PISM
echo
echo "--------- x_MOM-to-PISM --------"
echo " remove regridded MOM output, only keep a few"
REGRID_KEEP_FILES=$(echo `seq -f '%06g0101.regrid.MOM-to-PISM.bil.cdo.nc' \
    $(echo "$SIM_START_YEAR-$CPL_TIMESTEP" | bc -l) \
    $KEEP_FREQ_YEAR \
    $SIM_END_YEAR`)
echo files not to remove: $REGRID_KEEP_FILES
if [ $SIMULATE = false ]; then
    mkdir keep
    mv $REGRID_KEEP_FILES keep/
    rm *.regrid.MOM-to-PISM.bil.cdo.nc
    mv keep/* .
    rm -d keep
fi
