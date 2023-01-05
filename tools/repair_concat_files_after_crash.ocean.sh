#!/bin/bash

### repairing corrupted static fields in concatenated ocean output
#
#   The script 'concat_output_files.ocean.sh' concatenates MOM output that is
#   spread over multiple files (due to model restarts and possible crashes). It
#   is possible that corrupted static fields (like cell area, ...) are present
#   in the concatenated files.
#   To correct for this, this script copies the requested fields from a
#   correctly terminated job output (usually the last time slice) and adds them
#   to the concatenated output.

module load intel/2018.1 nco

TIMESTAMP_IN="0212250201"
TIMESTAMP_OUT="0212250201"

FILE_NAMES_OCN="ocean-decadal ocean-yearly"
VARIABLES_OCN="-v area_t,area_u,geolat_c,geolat_t,geolon_c,geolon_t,ht,hu,rescale_rho0_mask"

FILE_NAMES_ICE="ice-decadal"
VARIABLES_ICE="-v CELL_AREA,COSROT,GEOLON,GEOLAT,SINROT"

CMD="ncks -A "


for f in $FILE_NAMES_OCN; do
    set -x
    $CMD $VARIABLES_OCN intermediate_files/$TIMESTAMP_IN.$f.nc $TIMESTAMP_OUT.$f.nc
    set +x
done

for f in $FILE_NAMES_ICE; do
    set -x
    $CMD $VARIABLES_ICE intermediate_files/$TIMESTAMP_IN.$f.nc $TIMESTAMP_OUT.$f.nc
    set +x
done
