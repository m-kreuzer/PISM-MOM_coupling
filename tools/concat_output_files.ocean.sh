#!/bin/bash

### concatenating ocean output files from multiple consecutive runs
#
#   When MOM runs are restarted in between (e.g. due to crashes), the resulting
#   output is spread over multiple files with different output time stamps. In
#   order to make analysis easier, this script is concatenating the output into
#   single files.
#   Be aware that some files might have corrupted fields of static data (like
#   cell area, ...) if the model did not finish as planned but crashed earlier.
#   When concatenating, the corrupted fields might be copied to the
#   concatenated output. To correct this, use the script
#   'repair_concat_files_after_crash.ocean.sh`
#
#   Execute this script in $MOM_WORK_DIR/history

module load intel/2018.1 nco

TIMESTAMPS_IN="0190250201 0206250201 0212250201"
TIMESTAMP_OUT="0212250201"

FILE_NAMES="ice-decadal ice-decadal_max ice-decadal_min ocean-decadal ocean-decadal_min ocean-decadal_max ocean-scalar ocean-yearly"

CMD="ncrcat"

mkdir intermediate_files

for f in $FILE_NAMES; do
    FILE_LIST=''
    # put all files to be concated in a list
    for t in $TIMESTAMPS_IN; do
        FILE_LIST=$FILE_LIST' '$t.$f.nc
    done

    # move files to backup dir
    echo mv $FILE_LIST intermediate_files
    mv $FILE_LIST intermediate_files
    cd intermediate_files

    # append output file name to list
    FILE_LIST=$FILE_LIST' '../$TIMESTAMP_OUT.$f.nc
    set -x
    $CMD $FILE_LIST
    set +x

    cd ..

done
