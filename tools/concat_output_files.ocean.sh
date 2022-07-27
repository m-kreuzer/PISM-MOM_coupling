#!/bin/bash
# run script with 'bash concat_files.sh'

module load intel/2018.1 nco

TIMESTAMPS_IN="0141900101 0151900101 0171900101"
TIMESTAMP_OUT="0171900101"

FILE_NAMES="ice-decadal ocean-decadal ocean-scalar ocean-yearly"

CMD="ncrcat"


for f in $FILE_NAMES; do
    FILE_LIST=''
    # put all files to be concated in a list
    for t in $TIMESTAMPS_IN; do
        FILE_LIST=$FILE_LIST' '$t.$f.nc
    done
    # append output file name to list
    FILE_LIST=$FILE_LIST' '$TIMESTAMP_OUT.$f.nc.cat
    set -x
    $CMD $FILE_LIST
    set +x

    # rename concated file
    mv $TIMESTAMP_OUT.$f.nc $TIMESTAMP_OUT.$f.nc.bak
    mv $TIMESTAMP_OUT.$f.nc.cat $TIMESTAMP_OUT.$f.nc
done
