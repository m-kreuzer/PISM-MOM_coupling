#! /bin/bash

#  Copyright (C) 2019-2021 PISM-MOM_coupling authors, see AUTHORS file
#
#  This file is part of PISM-MOM_coupling
#
#  PISM-MOM_coupling is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  PISM-MOM_coupling is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with PISM-MOM_coupling.  If not, see <https://www.gnu.org/licenses/>.

# -----------------------------------------------------------------------------
# standalone version of process_mom_to_pism() and process_atm_to_pism()
# routines of PISM-MOM coupling framework
#
# This script calculates the PISM/PICO's input variables (temp, salt, eta) from
# ocean standalone output as well as atmospheric anomalies (precipitation,
# surface temperature) from atmospheric output.
# Use case: compare PISM/PICO input variables in standalone and coupled mode;
# or force standalone ice sheet model with pre-computed climate anomalies.
# The script requires some information for regridding and reference files,
# which are typically taken from coupled runs. They can be specified in
# process_poem_to_pism.settings
# -----------------------------------------------------------------------------

#SBATCH --account=poem
#SBATCH --job-name=process_poem_to_pism
#SBATCH --tasks=1
#SBATCH --mem=60000     # unit: MB
#SBATCH --qos=priority
##SBATCH --qos=short
#SBATCH --time=0:30:00
#SBATCH --output=process_poem_to_pism.%j.out
#SBATCH --error=process_poem_to_pism.%j.err
#SBATCH --mail-type=END,FAIL,REQUEUE,STAGE_OUT,TIME_LIMIT_90,TIME_LIMIT

TIME_START_SCRIPT=$(date +%s.%N)

# ---------------------------- function definitions ----------------------------
concat_poem_output(){
    echo
    echo
    echo " >> concat scattered outputs into single file"
    echo
    TIME_START_CONCAT_POEM_OUTPUT=$(date +%s.%N)

    cd $POEM_OUTPUT_DIR

    # get list of file name bases and check quantity
    bases=`ls -1 *.nc | sed -r "s/([[:digit:]]+)\.([^[:digit:]]*\.nc)/\2/" | sort -u`
    n_bases=$(echo $bases | awk -F' ' '{print NF; exit}')
    # get list of timestamps and check quantity
    timestamps=`ls -1 *.nc | sed -r "s/([[:digit:]]+)\.([^[:digit:]]*\.nc)/\1/" | sort -u`
    n_timestamps=$(echo $timestamps | awk -F' ' '{print NF; exit}')

    echo bases:      $bases
    echo timestamps: $timestamps
    echo

    if [ "$n_bases" -ge 1 ] && [ $n_timestamps -ge 2 ]; then
        echo "  > concatenating POEM output files along time dimension"

        mkdir -p intermediate_files

        for b in $bases; do
            file_list=''
            for t in $timestamps; do
                file_list="$file_list $t.$b"
            done
            echo mv $file_list intermediate_files
            mv $file_list intermediate_files
            cd intermediate_files
            echo ncrcat $file_list ../$t.$b
            ncrcat $file_list ../$t.$b
            cd ..
        done
        POEM_TIME_END=$t
    else
        echo "  > no concatenating required"
        if [ $n_timestamps -eq 1 ]; then
            POEM_TIME_END=$timestamps
        else
            for t in $timestamps; do
                POEM_TIME_END=$timestamps
            done
        fi
    fi

    cd $ROOT_WORK_DIR

    TIME_END_CONCAT_POEM_OUTPUT=$(date +%s.%N)
    TIME_CONCAT_POEM_OUTPUT=$(echo "$TIME_END_CONCAT_POEM_OUTPUT - $TIME_START_CONCAT_POEM_OUTPUT" | bc | awk '{printf "%f", $0}')
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
prepare_ocean_tracer_anomaly_reference_file(){
    echo
    echo
    echo " >> prepare ocean tracer anomaly reference file"
    echo
    TIME_START_OCN_TRACER_ANOMALY_REF=$(date +%s.%N)

    echo COMPUTE_OCEAN_TRACER_ANOMALY_FILE=$COMPUTE_OCEAN_TRACER_ANOMALY_FILE
    case $COMPUTE_OCEAN_TRACER_ANOMALY_FILE in
        true)
            echo "  > computing ocean tracer anomaly reference state from $CALC_OCN_TRACER_ANOMALY_INPUT_PATH"
            PWD_SAVE=$PWD
            cd $CALC_OCN_TRACER_ANOMALY_INPUT_PATH
            echo "   > average model output (possibly spread across multile files)"
            INPUT_FILES=$(echo `seq -f "$CALC_OCN_TRACER_ANOMALY_NAME_FORMAT_IN" \
                                       "$CALC_OCN_TRACER_ANOMALY_YR_START"  \
                                       "$CALC_OCN_TRACER_ANOMALY_YR_STEP"   \
                                       "$CALC_OCN_TRACER_ANOMALY_YR_END" `)
            AVG_OUT_FILE=$(echo `printf "$CALC_OCN_TRACER_ANOMALY_NAME_FORMAT_OUT" \
                                        "$CALC_OCN_TRACER_ANOMALY_YR_START"   \
                                        "$CALC_OCN_TRACER_ANOMALY_YR_END" `)
            AVG_OUT_PATH=$CALC_OCN_TRACER_ANOMALY_OUTPUT_PATH/$AVG_OUT_FILE
            ncra --overwrite -v temp,salt $INPUT_FILES $AVG_OUT_PATH
            RESULT=$?
            return_check $RESULT "ncra-operation"


            echo "   > bilinear regridding averaged MOM output to PISM grid"

            ## extract temp, salt for regridding
            #FILE_EXTENSION="${AVG_OUT_PATH##*.}"
            #FILE_PATH_NO_EXTENSION="${AVG_OUT_PATH%.*}"
            #REGRID_IN=$FILE_PATH_NO_EXTENSION.sparse.$FILE_EXTENSION
            #ncks -O -v temp,salt $AVG_OUT_PATH tmp.nc
            #ncks -C -O -x -v geolat_t,geolon_t tmp.nc $REGRID_IN
            #rm tmp.nc

            REGRID_IN=$AVG_OUT_PATH
            ncks -C -O -x -v geolat_t,geolon_t $REGRID_IN tmp.nc
            RESULT=$?
            return_check $RESULT "ncks-operation"
            mv tmp.nc $REGRID_IN

            REGRID_OUT=$CALC_OCN_TRACER_ANOMALY_OUTPUT_PATH/ocean_tracer_mean.regrid.MOM-to-PISM.bil.cdo.nc
            # remove global attribute _NCProperties from POEM output which causes cdo to crash
            ncatted -O -a _NCProperties,global,d,, $REGRID_IN
            #export REMAP_EXTRAPOLATE=off
            cdo -b F64 -f nc4c remap,$PISM_PRE_OUT_FILE,$WEIGHTS_OCN_PATH $REGRID_IN $REGRID_OUT
            RESULT=$?
            return_check $RESULT "cdo_regridding"

            OUTFILE=$(echo `printf "%06g-%06g.tracer_mean.processed_MOM.nc" \
                                   "$CALC_OCN_TRACER_ANOMALY_YR_START"      \
                                   "$CALC_OCN_TRACER_ANOMALY_YR_END" `)
            OCEAN_TRACER_ANOMALY_REFERENCE_FILE=$CALC_OCN_TRACER_ANOMALY_OUTPUT_PATH/$OUTFILE
            echo "  > regriddedMOM-to-PISM_processing script"
            $SCRIPT_DIR/regriddedMOM-to-PISM_processing.py          \
                --input     $REGRID_OUT                             \
                --basins    $PISM_BASIN_FILE                        \
                --edges     $PISM_EDGES_FILE                        \
                --fill      temp salt                               \
                --depth     $PICO_INPUT_DEPTH_FILE                  \
                --output    $OCEAN_TRACER_ANOMALY_REFERENCE_FILE    \
                --verbose
            RESULT=$?
            return_check $RESULT "regriddedMOM-to-PISM_processing.py"

            echo OCEAN_TRACER_ANOMALY_REFERENCE_FILE=$OCEAN_TRACER_ANOMALY_REFERENCE_FILE
            cd $PWD_SAVE

        ;;

        false)
            echo "  > using ocean prescribed tracer anomaly reference file"
            echo OCEAN_TRACER_ANOMALY_REFERENCE_FILE=$OCEAN_TRACER_ANOMALY_REFERENCE_FILE

            # copy reference file to own evaluation directory and update path
            cp -av $OCEAN_TRACER_ANOMALY_REFERENCE_FILE $CALC_OCN_TRACER_ANOMALY_OUTPUT_PATH/
            RESULT=$?
            return_check $RESULT "cp"

            FILE_NAME="$(basename "${OCEAN_TRACER_ANOMALY_REFERENCE_FILE}")"
            OCEAN_TRACER_ANOMALY_REFERENCE_FILE=$CALC_OCN_TRACER_ANOMALY_OUTPUT_PATH/$FILE_NAME
            echo OCEAN_TRACER_ANOMALY_REFERENCE_FILE=$OCEAN_TRACER_ANOMALY_REFERENCE_FILE

        ;;

        *)
            echo COMPUTE_OCEAN_TRACER_ANOMALY_FILE is supposed to be either "true" or "false"
            echo given value: COMPUTE_OCEAN_TRACER_ANOMALY_FILE=$COMPUTE_OCEAN_TRACER_ANOMALY_FILE
            exit -1
        ;;
    esac

    # copy temp,salt fields, so that basin mean fields are added by ncbo calls below
    ncap2 -O -s 'theta_ocean_basin_mean=theta_ocean;salinity_ocean_basin_mean=salinity_ocean' \
        $OCEAN_TRACER_ANOMALY_REFERENCE_FILE $OCEAN_TRACER_ANOMALY_REFERENCE_FILE.tmp
    RESULT=$?
    return_check $RESULT "ncap2-operation"
    mv $OCEAN_TRACER_ANOMALY_REFERENCE_FILE.tmp $OCEAN_TRACER_ANOMALY_REFERENCE_FILE

    # removing time dimension from ocean tracer anomaly reference file
    #  -> needed for clean broadcasting with ncbo
    FILE_EXTENSION="${OCEAN_TRACER_ANOMALY_REFERENCE_FILE##*.}"
    FILE_PATH_NO_EXTENSION="${OCEAN_TRACER_ANOMALY_REFERENCE_FILE%.*}"
    OCEAN_TRACER_ANOMALY_REFERENCE_FILE_NOTIME=$FILE_PATH_NO_EXTENSION.no-time.$FILE_EXTENSION

    echo OCEAN_TRACER_ANOMALY_REFERENCE_FILE_NOTIME=$OCEAN_TRACER_ANOMALY_REFERENCE_FILE_NOTIME
    ncwa -O -a time $OCEAN_TRACER_ANOMALY_REFERENCE_FILE $OCEAN_TRACER_ANOMALY_REFERENCE_FILE_NOTIME
    RESULT=$?
    return_check $RESULT "ncwa-operation"


    TIME_END_OCN_TRACER_ANOMALY_REF=$(date +%s.%N)
    TIME_OCN_TRACER_ANOMALY_REF=$(echo "$TIME_END_OCN_TRACER_ANOMALY_REF - $TIME_START_OCN_TRACER_ANOMALY_REF" | bc | awk '{printf "%f", $0}')

}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
prepare_ocean_sealevel_anomaly_reference_file(){
    echo
    echo
    echo " >> prepare ocean sealevel anomaly reference file"
    echo
    TIME_START_OCN_SL_ANOMALY_REF=$(date +%s.%N)

    echo COMPUTE_OCEAN_SEALEVEL_ANOMALY_FILE=$COMPUTE_OCEAN_SEALEVEL_ANOMALY_FILE
    case $COMPUTE_OCEAN_SEALEVEL_ANOMALY_FILE in
        true)
            echo "  > computing ocean sealevel anomaly reference state from $CALC_OCN_SEALEVEL_ANOMALY_INPUT_PATH"
            PWD_SAVE=$PWD
            cd $CALC_OCN_SEALEVEL_ANOMALY_INPUT_PATH
            echo "   > average model output (possibly spread across multile files)"
            INPUT_FILES=$(echo `seq -f "$CALC_OCN_SEALEVEL_ANOMALY_NAME_FORMAT_IN"  \
                                       "$CALC_OCN_SEALEVEL_ANOMALY_YR_START"        \
                                       "$CALC_OCN_SEALEVEL_ANOMALY_YR_STEP"         \
                                       "$CALC_OCN_SEALEVEL_ANOMALY_YR_END" `)
            AVG_OUT_FILE=$(echo `printf "$CALC_OCN_SEALEVEL_ANOMALY_NAME_FORMAT_OUT" \
                                        "$CALC_OCN_SEALEVEL_ANOMALY_YR_START"        \
                                        "$CALC_OCN_SEALEVEL_ANOMALY_YR_END" `)
            AVG_OUT_PATH=$CALC_OCN_SEALEVEL_ANOMALY_INPUT_PATH/$AVG_OUT_FILE
            ncra --overwrite -v eta_t $INPUT_FILES $AVG_OUT_PATH
            RESULT=$?
            return_check $RESULT "ncra-operation"


            echo "   > bilinear regridding averaged MOM output to PISM grid"

            ## extract temp, salt for regridding
            #FILE_EXTENSION="${AVG_OUT_PATH##*.}"
            #FILE_PATH_NO_EXTENSION="${AVG_OUT_PATH%.*}"
            #REGRID_IN=$FILE_PATH_NO_EXTENSION.sparse.$FILE_EXTENSION
            #ncks -O -v temp,salt $AVG_OUT_PATH tmp.nc
            #ncks -C -O -x -v geolat_t,geolon_t tmp.nc $REGRID_IN
            #rm tmp.nc

            REGRID_IN=$AVG_OUT_PATH
            ncks -C -O -x -v geolat_t,geolon_t $REGRID_IN tmp.nc
            RESULT=$?
            return_check $RESULT "ncks-operation"
            mv tmp.nc $REGRID_IN

            REGRID_OUT=$CALC_OCN_SEALEVEL_ANOMALY_INPUT_PATH/ocean_sealevel_mean.regrid.MOM-to-PISM.bil.cdo.nc
            # remove global attribute _NCProperties from POEM output which causes cdo to crash
            ncatted -O -a _NCProperties,global,d,, $REGRID_IN
            cdo -b F64 -f nc4c remap,$PISM_PRE_OUT_FILE,$WEIGHTS_OCN_PATH $REGRID_IN $REGRID_OUT
            RESULT=$?
            return_check $RESULT "cdo_regridding"

            OUTFILE=$(echo `printf "%06g-%06g.sealevel_mean.processed_MOM.nc" \
                                   "$CALC_OCN_SEALEVEL_ANOMALY_YR_START"      \
                                   "$CALC_OCN_SEALEVEL_ANOMALY_YR_END" `)
            OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE=$CALC_OCN_SEALEVEL_ANOMALY_INPUT_PATH/$OUTFILE
            echo "  > regriddedMOM-to-PISM_processing script"
            $SCRIPT_DIR/regriddedMOM-to-PISM_processing.py          \
                --input     $REGRID_OUT                             \
                --basins    $PISM_BASIN_FILE                        \
                --edges     $PISM_EDGES_FILE                        \
                --fill      eta_t                                   \
                --depth     $PICO_INPUT_DEPTH_FILE                  \
                --output    $OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE  \
                --verbose
            RESULT=$?
            return_check $RESULT "regriddedMOM-to-PISM_processing.py"

            echo OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE=$OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE
            cd $PWD_SAVE

        ;;

        false)
            echo "  > using ocean prescribed sealevel anomaly reference file"
            echo OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE=$OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE

            # copy reference file to own evaluation directory and update path
            cp -av $OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE $CALC_OCN_SEALEVEL_ANOMALY_OUTPUT_PATH/
            RESULT=$?
            return_check $RESULT "cp"

            FILE_NAME="$(basename "${OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE}")"
            OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE=$CALC_OCN_SEALEVEL_ANOMALY_OUTPUT_PATH/$FILE_NAME
            echo OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE=$OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE

        ;;

        *)
            echo COMPUTE_OCEAN_SEALEVEL_ANOMALY_FILE is supposed to be either "true" or "false"
            echo given value: COMPUTE_OCEAN_SEALEVEL_ANOMALY_FILE=$COMPUTE_OCEAN_SEALEVEL_ANOMALY_FILE
            exit -1
        ;;
    esac

    # copy temp,salt fields, so that basin mean fields are added by ncbo calls below
    ncap2 -O -s 'delta_SL_basin_mean=delta_SL' \
        $OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE $OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE.tmp
    RESULT=$?
    return_check $RESULT "ncap2-operation"
    mv $OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE.tmp $OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE

    # removing time dimension from ocean tracer anomaly reference file
    #  -> needed for clean broadcasting with ncbo
    FILE_EXTENSION="${OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE##*.}"
    FILE_PATH_NO_EXTENSION="${OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE%.*}"
    OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE_NOTIME=$FILE_PATH_NO_EXTENSION.no-time.$FILE_EXTENSION

    echo OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE_NOTIME=$OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE_NOTIME
    ncwa -O -a time $OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE $OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE_NOTIME
    RESULT=$?
    return_check $RESULT "ncwa-operation"

    TIME_END_OCN_SL_ANOMALY_REF=$(date +%s.%N)
    TIME_OCN_SL_ANOMALY_REF=$(echo "$TIME_END_OCN_SL_ANOMALY_REF - $TIME_START_OCN_SL_ANOMALY_REF" | bc | awk '{printf "%f", $0}')
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
prepare_atmos_anomaly_reference_file(){
    echo
    echo
    echo " >> prepare atmosphere anomaly reference file"
    echo
    TIME_START_ATM_ANOMALY_REF=$(date +%s.%N)

    echo COMPUTE_OCEAN_TRACER_ANOMALY_FILE=$COMPUTE_OCEAN_TRACER_ANOMALY_FILE
    case $COMPUTE_OCEAN_TRACER_ANOMALY_FILE in
        true)
            echo "  > computing ocean tracer anomaly reference state from $CALC_ATM_ANOMALY_INPUT_PATH"
            PWD_SAVE=$PWD
            cd $CALC_ATM_ANOMALY_INPUT_PATH
            echo "   > average model output (possibly spread across multile files)"
            INPUT_FILES=$(echo `seq -f "$CALC_ATM_ANOMALY_NAME_FORMAT_IN" \
                                       "$CALC_ATM_ANOMALY_YR_START"  \
                                       "$CALC_ATM_ANOMALY_YR_STEP"   \
                                       "$CALC_ATM_ANOMALY_YR_END" `)
            AVG_OUT_FILE=$(echo `printf "$CALC_ATM_ANOMALY_NAME_FORMAT_OUT" \
                                        "$CALC_ATM_ANOMALY_YR_START"   \
                                        "$CALC_ATM_ANOMALY_YR_END" `)
            AVG_OUT_PATH=$CALC_ATM_ANOMALY_OUTPUT_PATH/$AVG_OUT_FILE
            ncra --overwrite -v precip,t_surf $INPUT_FILES $AVG_OUT_PATH
            RESULT=$?
            return_check $RESULT "ncra-operation"


            echo "   > bilinear regridding averaged ATM output to PISM grid"

            ## extract temp, salt for regridding
            #FILE_EXTENSION="${AVG_OUT_PATH##*.}"
            #FILE_PATH_NO_EXTENSION="${AVG_OUT_PATH%.*}"
            #REGRID_IN=$FILE_PATH_NO_EXTENSION.sparse.$FILE_EXTENSION
            #ncks -O -v temp,salt $AVG_OUT_PATH tmp.nc
            #ncks -C -O -x -v geolat_t,geolon_t tmp.nc $REGRID_IN
            #rm tmp.nc

            #REGRID_IN=$AVG_OUT_PATH
            #ncks -C -O -x -v geolat_t,geolon_t $REGRID_IN tmp.nc
            #RESULT=$?
            #return_check $RESULT "ncks-operation"
            #mv tmp.nc $REGRID_IN

            REGRID_OUT=$CALC_ATM_ANOMALY_OUTPUT_PATH/atmos_reference_mean.regrid.ATM-to-PISM.bil.cdo.nc
            # remove global attribute _NCProperties from POEM output which causes cdo to crash
            #ncatted -O -a _NCProperties,global,d,, $REGRID_IN
            #export REMAP_EXTRAPOLATE=on
            cdo -b F64 -f nc4c remap,$PISM_PRE_OUT_FILE,$WEIGHTS_ATM_PATH $AVG_OUT_PATH $REGRID_OUT
            RESULT=$?
            return_check $RESULT "cdo_regridding"
            ATMOS_ANOMALY_REFERENCE_FILE=$REGRID_OUT
            cd $PWD_SAVE

        ;;

        false)
            echo "  > using ocean prescribed tracer anomaly reference file"
            echo ATMOS_ANOMALY_REFERENCE_FILE=$ATMOS_ANOMALY_REFERENCE_FILE

            # copy reference file to own evaluation directory and update path
            cp -av $ATMOS_ANOMALY_REFERENCE_FILE $CALC_ATM_ANOMALY_OUTPUT_PATH/
            RESULT=$?
            return_check $RESULT "cp"

            FILE_NAME="$(basename "${ATMOS_ANOMALY_REFERENCE_FILE}")"
            ATMOS_ANOMALY_REFERENCE_FILE=$CALC_ATM_ANOMALY_OUTPUT_PATH/$FILE_NAME
            echo ATMOS_ANOMALY_REFERENCE_FILE=$ATMOS_ANOMALY_REFERENCE_FILE
        ;;

        *)
            echo COMPUTE_ATMOS_ANOMALY_FILE is supposed to be either "true" or "false"
            echo given value: COMPUTE_ATMOS_ANOMALY_FILE=$COMPUTE_ATMOS_ANOMALY_FILE
            exit -1
        ;;
    esac

    # removing time dimension from atmosphere anomaly reference file
    #  -> needed for clean broadcasting with ncbo
    FILE_EXTENSION="${ATMOS_ANOMALY_REFERENCE_FILE##*.}"
    FILE_PATH_NO_EXTENSION="${ATMOS_ANOMALY_REFERENCE_FILE%.*}"
    ATMOS_ANOMALY_REFERENCE_FILE_NOTIME=$FILE_PATH_NO_EXTENSION.no-time.$FILE_EXTENSION

    echo ATMOS_ANOMALY_REFERENCE_FILE_NOTIME=$ATMOS_ANOMALY_REFERENCE_FILE_NOTIME
    ncwa -O -a time $ATMOS_ANOMALY_REFERENCE_FILE $ATMOS_ANOMALY_REFERENCE_FILE_NOTIME
    RESULT=$?
    return_check $RESULT "ncwa-operation"


    TIME_END_ATM_ANOMALY_REF=$(date +%s.%N)
    TIME_ATM_ANOMALY_REF=$(echo "$TIME_END_ATM_ANOMALY_REF - $TIME_START_ATM_ANOMALY_REF" | bc | awk '{printf "%f", $0}')

}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
process_mom_to_pism(){
    echo
    echo
    echo " >> process MOM output for input in PISM"
    echo
    TIME_START_MOM_to_PISM_PROCESS=$(date +%s.%N)

    echo "  > bilinear regridding MOM output to PISM grid"

   OCEAN_OUT=$POEM_OUTPUT_DIR/$POEM_TIME_END.$OCN_OUTPUT_FILE_BASE.nc

    ### average yearly ocean output over bins of 10 years (decadal resolution) or not
    case $AVERAGE_OCN_OUTPUT_YEARLY_TO_DECADAL in
        true)
            echo "decrease temporal resolution of ocean output from yearly to decadal"
            echo OCN_OUTPUT_FILE_BASE: $OCN_OUTPUT_FILE_BASE

            ncra -O --mro -d time,,,10,10 -v temp,salt,eta_t $OCEAN_OUT tmp.nc
            RESULT=$?
            return_check $RESULT "ncra command"

            OCEAN_OUT_MOD=$X_MOM_TO_PISM_PATH/$POEM_TIME_END.$OCN_OUTPUT_FILE_BASE_DECADAL.nc
            ncap2 -O -s 'defdim("bnds",2); time_bounds=make_bounds(time,$bnds,"time_bounds");' \
                tmp.nc $OCEAN_OUT_MOD
            RESULT=$?
            return_check $RESULT "ncap2 command"
            rm tmp.nc 
        ;;

        false)
            echo "don't change temporal resolution of ocean output"
            echo OCN_OUTPUT_FILE_BASE: $OCN_OUTPUT_FILE_BASE

            OCEAN_OUT_MOD=$X_MOM_TO_PISM_PATH/$POEM_TIME_END.$OCN_OUTPUT_FILE_BASE.nc
            ncks -O -v temp,salt,eta_t,time_bounds $OCEAN_OUT $OCEAN_OUT_MOD
            RESULT=$?
            return_check $RESULT "ncks command"
        ;;

        *)
            echo AVERAGE_OCN_OUTPUT_YEARLY_TO_DECADAL is supposed to be either "true" or "false"
            echo given value: AVERAGE_OCN_OUTPUT_YEARLY_TO_DECADAL=$AVERAGE_OCN_OUTPUT_YEARLY_TO_DECADAL
            exit -1
        ;;
    esac

    ## extract temp, salt, eta_t for regridding
    #OCEAN_OUT=$POEM_OUTPUT_DIR/$POEM_TIME_END.$OCN_OUTPUT_FILE_BASE.nc
    #REGRID_IN=$X_MOM_TO_PISM_PATH/$POEM_TIME_END.$OCN_OUTPUT_FILE_BASE.sparse.nc
    #ncks -O -v temp,salt,eta_t,time_bounds $OCEAN_OUT tmp.nc
    #ncks -C -O -x -v geolat_t,geolon_t tmp.nc $REGRID_IN
    #rm tmp.nc

    REGRID_OUT=$X_MOM_TO_PISM_PATH/$POEM_TIME_END.regrid.MOM-to-PISM.bil.cdo.nc
    # remove global attribute _NCProperties from POEM output which causes cdo to crash
    ncatted -O -a _NCProperties,global,d,, $REGRID_IN
    #export REMAP_EXTRAPOLATE=off
    #cdo -b F64 -f nc4c remap,$PISM_PRE_OUT_FILE,$WEIGHTS_OCN_PATH $REGRID_IN $REGRID_OUT
    cdo -b F64 -f nc4c remap,$PISM_PRE_OUT_FILE,$WEIGHTS_OCN_PATH $OCEAN_OUT_MOD $REGRID_OUT
    RESULT=$?
    return_check $RESULT "cdo_regridding"

    echo "  > regriddedMOM-to-PISM_processing script"
    $SCRIPT_DIR/regriddedMOM-to-PISM_processing.py                                   \
        --input     $REGRID_OUT                                                      \
        --basins    $PISM_BASIN_FILE                                                 \
        --edges     $PISM_EDGES_FILE                                                 \
        --fill      temp salt eta_t                                                  \
        --depth     $PICO_INPUT_DEPTH_FILE                                           \
        --output    $X_MOM_TO_PISM_PATH/$POEM_TIME_END.processed_MOM.nc     \
        --verbose
    RESULT=$?
    return_check $RESULT "regriddedMOM-to-PISM_processing.py"


    ### ocean tracer anomaly approach: use regular PISM ocean forcing (e.g.
    #       Schmidtko data) but add difference of
    #       regriddedMOM-to-PISM_processing.py output from current coupling
    #       time step to reference file

    # difference from current ocean output to reference (temp, salt)
    ncbo -O -v theta_ocean,salinity_ocean,theta_ocean_basin_mean,salinity_ocean_basin_mean \
        --op_typ=subtract \
        $X_MOM_TO_PISM_PATH/$POEM_TIME_END.processed_MOM.nc \
        $OCEAN_TRACER_ANOMALY_REFERENCE_FILE_NOTIME \
        $X_MOM_TO_PISM_PATH/$POEM_TIME_END.tracer.processed_MOM.anomaly.nc
    RESULT=$?
    return_check $RESULT "ncbo"
    # copy time and basins variables to PISM forcing file (omitted by ncbo add operation above)
    ncks -A -v time,time_bounds,basins \
        $X_MOM_TO_PISM_PATH/$POEM_TIME_END.processed_MOM.nc \
        $X_MOM_TO_PISM_PATH/$POEM_TIME_END.tracer.processed_MOM.anomaly.nc
    RESULT=$?
    return_check $RESULT "ncks"

    # add difference to given ocean forcing file (temp, salt)
    cp -av $PISM_OCN_FORCING $X_MOM_TO_PISM_PATH
    RESULT=$?
    return_check $RESULT "cp"

    PISM_OCN_FORCING_FILE=$(basename "$PISM_OCN_FORCING")
    PISM_OCN_FORCING_PATH=$X_MOM_TO_PISM_PATH/$PISM_OCN_FORCING_FILE
    FILE_EXTENSION="${PISM_OCN_FORCING_PATH##*.}"
    FILE_PATH_NO_EXTENSION="${PISM_OCN_FORCING_PATH%.*}"
    PISM_OCN_FORCING_MOD_PATH=$FILE_PATH_NO_EXTENSION.mod.$FILE_EXTENSION
    #ncap2 -O -s 'theta_ocean_basin_mean=theta_ocean;salinity_ocean_basin_mean=salinity_ocean' \
    #    $ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$PISM_OCN_FORCING_FILE \
    #    $ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$PISM_OCN_FORCING_FILE.mod
    #ncbo -O --op_typ=add \
    #    $ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$POEM_TIME_END.tracer.processed_MOM.anomaly.nc \
    #    $ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$PISM_OCN_FORCING_FILE.mod \
    #    $ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$POEM_TIME_END.tracer.PISM_input.nc
    ncap2 -O -s 'theta_ocean_basin_mean=theta_ocean;salinity_ocean_basin_mean=salinity_ocean' \
        $PISM_OCN_FORCING_PATH      \
        $PISM_OCN_FORCING_MOD_PATH
    RESULT=$?
    return_check $RESULT "ncap2"
    ncbo -O --op_typ=add \
        $X_MOM_TO_PISM_PATH/$POEM_TIME_END.tracer.processed_MOM.anomaly.nc \
        $PISM_OCN_FORCING_MOD_PATH  \
        $X_MOM_TO_PISM_PATH/$POEM_TIME_END.tracer.PISM_input.nc
    RESULT=$?
    return_check $RESULT "ncbo"
    # copy time and basins variables to PISM forcing file (omitted by ncbo add operation above)
    ncks -A -v time,time_bounds,basins \
        $X_MOM_TO_PISM_PATH/$POEM_TIME_END.processed_MOM.nc \
        $X_MOM_TO_PISM_PATH/$POEM_TIME_END.tracer.PISM_input.nc
    RESULT=$?
    return_check $RESULT "ncks"




    if [ $DO_OCEAN_SEALEVEL_ANOMALY = true ]; then
        ### ocean sealevel anomaly approach: force PISM with sea level change
        #       from ocean model. Calculated through difference of current local
        #       sea level (output of regriddedMOM-to-PISM_processing.py) to
        #       reference file.
        echo OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE=$OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE
        echo OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE_NOTIME=$OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE_NOTIME

        # difference from current ocean output to reference (temp, salt)
        ncbo -O -v delta_SL,delta_SL_basin_mean --op_typ=subtract \
            $X_MOM_TO_PISM_PATH/$POEM_TIME_END.processed_MOM.nc \
            $OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE_NOTIME \
            $X_MOM_TO_PISM_PATH/$POEM_TIME_END.sealevel.processed_MOM.anomaly.nc
        RESULT=$?
        return_check $RESULT "ncbo"
        # copy time variable to PISM forcing file (omitted by ncbo add operation above)
        ncks -A -v time \
            $X_MOM_TO_PISM_PATH/$POEM_TIME_END.processed_MOM.nc \
            $X_MOM_TO_PISM_PATH/$POEM_TIME_END.sealevel.processed_MOM.anomaly.nc
        RESULT=$?
        return_check $RESULT "ncks"

        # set NaNs to 0
        ncatted -O -a _FillValue,delta_SL,m,f,0 \
            $X_MOM_TO_PISM_PATH/$POEM_TIME_END.sealevel.processed_MOM.anomaly.nc
        RESULT=$?
        return_check $RESULT "ncatted"
        ncatted -O -a _FillValue,delta_SL,d,, \
            $X_MOM_TO_PISM_PATH/$POEM_TIME_END.sealevel.processed_MOM.anomaly.nc
        RESULT=$?
        return_check $RESULT "ncatted"

        # create symlink from PISM input file to sealevel anomaly file
        PWD_SAVE=$PWD
        cd $X_MOM_TO_PISM_PATH
        ln -sf $POEM_TIME_END.sealevel.processed_MOM.anomaly.nc $POEM_TIME_END.sealevel.PISM_input.nc
    fi
    cd $PWD_SAVE

    TIME_END_MOM_to_PISM_PROCESS=$(date +%s.%N)
    TIME_MOM_to_PISM_PROCESS=$(echo "$TIME_END_MOM_to_PISM_PROCESS - $TIME_START_MOM_to_PISM_PROCESS" | bc | awk '{printf "%f", $0}')
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
process_atm_to_pism() {
    echo
    echo
    echo " >> process ATM output for input in PISM"
    echo
    TIME_START_ATM_to_PISM_PROCESS=$(date +%s.%N)

    echo "  > bilinear regridding ATM output to PISM grid"

    ATMOS_OUT=$POEM_OUTPUT_DIR/$POEM_TIME_END.$ATM_OUTPUT_FILE_BASE.nc

    ### average yearly atmospheric output over bins of 10 years (decadal resolution) or not
    case $AVERAGE_ATM_OUTPUT_YEARLY_TO_DECADAL in
        true)
            echo "decrease temporal resolution of atmospheric output from yearly to decadal"
            echo ATM_OUTPUT_FILE_BASE: $ATM_OUTPUT_FILE_BASE

            ncra -O --mro -d time,,,10,10 -v t_surf,precip $ATMOS_OUT tmp.nc
            RESULT=$?
            return_check $RESULT "ncra command"

            ATMOS_OUT_MOD=$X_ATM_TO_PISM_PATH/$POEM_TIME_END.$ATM_OUTPUT_FILE_BASE_DECADAL.nc
            ncap2 -O -s 'defdim("bnds",2); time_bounds=make_bounds(time,$bnds,"time_bounds");' \
                tmp.nc $ATMOS_OUT_MOD
            RESULT=$?
            return_check $RESULT "ncap2 command"
            rm tmp.nc
        ;;

        false)
            echo "don't change temporal resolution of atmospheric output"
            echo ATM_OUTPUT_FILE_BASE: $ATM_OUTPUT_FILE_BASE

            ATMOS_OUT_MOD=$X_ATM_TO_PISM_PATH/$POEM_TIME_END.$ATM_OUTPUT_FILE_BASE.nc
            ncks -O -v precip,t_surf,time_bounds $ATMOS_OUT $ATMOS_OUT_MOD
            RESULT=$?
            return_check $RESULT "ncks command"
        ;;

        *)
            echo AVERAGE_ATM_OUTPUT_YEARLY_TO_DECADAL is supposed to be either "true" or "false"
            echo given value: AVERAGE_ATM_OUTPUT_YEARLY_TO_DECADAL=$AVERAGE_ATM_OUTPUT_YEARLY_TO_DECADAL
            exit -1
        ;;
    esac

    # average atmospheric output over coupling time step before regridding
    #REGRID_IN=$POEM_TIME_END.atmos-yearly.mean.nc
    #ncra -O -v t_surf,precip -F -d time,1,,1 $ATMOS_OUT $REGRID_IN
    REGRID_OUT=$X_ATM_TO_PISM_PATH/$POEM_TIME_END.regrid.ATM-to-PISM.bil.cdo.nc
    ## remove global attribute _NCProperties from POEM output which causes cdo to crash
    #ncatted -O -a _NCProperties,global,d,, $REGRID_IN
    #export REMAP_EXTRAPOLATE=on
    cdo -b F64 -f nc4c remap,$PISM_PRE_OUT_FILE,$WEIGHTS_ATM_PATH $ATMOS_OUT_MOD $REGRID_OUT
    RESULT=$?
    return_check $RESULT "cdo_regridding"

    ### atmosphere anomaly approach: use regular PISM atmosphere forcing (e.g.
    #       RACMO data) but add difference of atmosphere output from current coupling
    #       time step to reference file

    # difference from current atmosphere output to reference
    ncbo -O -v precip,t_surf --op_typ=subtract \
        $REGRID_OUT \
        $ATMOS_ANOMALY_REFERENCE_FILE_NOTIME \
        $X_ATM_TO_PISM_PATH/$POEM_TIME_END.atmos_anomaly.absolute.nc
    RESULT=$?
    return_check $RESULT "ncbo"

    # copy time variable as this seems to be omitted by ncbo
    ncks -A -v time $REGRID_OUT \
        $X_ATM_TO_PISM_PATH/$POEM_TIME_END.atmos_anomaly.absolute.nc
    RESULT=$?
    return_check $RESULT "ncks"

    ### apply relative precipitation anomaly from POEM to PISM baseline forcing
    ncrename -O -v precip,precip_anomaly_POEM -v t_surf,air_temp_anomaly \
        $X_ATM_TO_PISM_PATH/$POEM_TIME_END.atmos_anomaly.absolute.nc \
        tmp1.nc
    RESULT=$?
    return_check $RESULT "ncrename"

    ncrename -O -v precip,precip_reference_POEM \
        $ATMOS_ANOMALY_REFERENCE_FILE_NOTIME \
        tmp2.nc
    RESULT=$?
    return_check $RESULT "ncrename"

    ncks -A -v precip_reference_POEM tmp2.nc tmp1.nc
    RESULT=$?
    return_check $RESULT "ncks"

    ncap2 -O -s precip_anomaly_relative_POEM=precip_anomaly_POEM/precip_reference_POEM \
        tmp1.nc tmp3.nc
    RESULT=$?
    return_check $RESULT "ncap2"

    # preprocess pism atm file (otherwise follow up operations fail)
    ncwa -O -a time -v precipitation \
        $PISM_ATM_FORCING \
        tmp4.nc
    RESULT=$?
    return_check $RESULT "ncwa"

    ncrename -O -v precipitation,precip_pism_base tmp4.nc tmp5.nc
    RESULT=$?
    return_check $RESULT "ncrename"

    ncks -A -v precip_pism_base tmp5.nc tmp3.nc
    RESULT=$?
    return_check $RESULT "ncks"

    ncap2 -O -s precipitation_anomaly=precip_pism_base*precip_anomaly_relative_POEM \
        tmp3.nc $X_ATM_TO_PISM_PATH/$POEM_TIME_END.atmos_anomaly.PISM_input.nc
    RESULT=$?
    return_check $RESULT "ncap2"

    rm tmp*.nc


    ### some sanity checks
    VMAX=$(ncmax air_temp_anomaly $X_ATM_TO_PISM_PATH/$POEM_TIME_END.atmos_anomaly.PISM_input.nc)
    if (( $(echo "$VMAX > 100" |bc -l) ))
    then
        echo "[sanity check] Air temperature anomaly in file "\
            "$POEM_TIME_END.atmos_anomaly.PISM_input.nc is greater than 100 degC. "\
            "Exiting as this seems abnormal."
        exit_script 4
    fi

    VMIN=$(ncmin air_temp_anomaly $X_ATM_TO_PISM_PATH/$POEM_TIME_END.atmos_anomaly.PISM_input.nc)
    if (( $(echo "$VMIN < -100" |bc -l) ))
    then
        echo "[sanity check] Air temperature anomaly in file "\
            "$POEM_TIME_END.atmos_anomaly.PISM_input.nc is smaller than -100 degC. "\
            "Exiting as this seems abnormal."
        exit_script 4
    fi

    VMAX=$(ncmax precipitation_anomaly $X_ATM_TO_PISM_PATH/$POEM_TIME_END.atmos_anomaly.PISM_input.nc)
    if (( $(echo "$VMAX > 100000" |bc -l) ))
    then
        echo "[sanity check] Precipitation anomaly in file "\
            "$POEM_TIME_END.atmos_anomaly.PISM_input.nc is greater than 100000 kg/m2/year. "\
            "Exiting as this seems abnormal."
        exit_script 4
    fi

    VMIN=$(ncmin precipitation_anomaly $X_ATM_TO_PISM_PATH/$POEM_TIME_END.atmos_anomaly.PISM_input.nc)
    if (( $(echo "$VMIN < -100000" |bc -l) ))
    then
        echo "[sanity check] Precipitation anomaly in file "\
            "$POEM_TIME_END.atmos_anomaly.PISM_input.nc is smaller than -100000 kg/m2/year. "\
            "Exiting as this seems abnormal."
        exit_script 4
    fi


    TIME_END_ATM_to_PISM_PROCESS=$(date +%s.%N)
    TIME_ATM_to_PISM_PROCESS=$(echo "$TIME_END_ATM_to_PISM_PROCESS - $TIME_START_ATM_to_PISM_PROCESS" | bc | awk '{printf "%f", $0}')
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

calc_basin_contshelf_mean(){
    echo
    echo
    echo " >> calculate continental shelf mean values for each basin"
    echo

    TIME_START_CALC_BASIN_CONTSHELF_MEAN=$(date +%s.%N)


    echo "  > ocean output (processed to PISM grid)"
    echo
    OCEAN_OUTPUT_FILE=$X_MOM_TO_PISM_PATH/$POEM_TIME_END.processed_MOM.nc
    FILE_EXTENSION="${OCEAN_OUTPUT_FILE##*.}"
    FILE_PATH_NO_EXTENSION="${OCEAN_OUTPUT_FILE%.*}"
    OCEAN_OUTPUT_FILE_BASIN_MEAN=$FILE_PATH_NO_EXTENSION.basin_mean.$FILE_EXTENSION

    $SCRIPT_DIR/calc_basin_contshelf_mean.py                        \
        --input     $OCEAN_OUTPUT_FILE                              \
        --basins    $PISM_BASIN_FILE                                \
        --fill      theta_ocean salinity_ocean delta_SL             \
        --output    $OCEAN_OUTPUT_FILE_BASIN_MEAN                   \
        --verbose
    RESULT=$?
    return_check $RESULT "calc_basin_contshelf_mean"

    echo
    echo "  > tracer reference (baseline for computing anomalies)"
    echo
    FILE_EXTENSION="${OCEAN_TRACER_ANOMALY_REFERENCE_FILE##*.}"
    FILE_PATH_NO_EXTENSION="${OCEAN_TRACER_ANOMALY_REFERENCE_FILE%.*}"
    OCEAN_TRACER_ANOMALY_REFERENCE_FILE_BASIN_MEAN=$FILE_PATH_NO_EXTENSION.basin_mean.$FILE_EXTENSION

    $SCRIPT_DIR/calc_basin_contshelf_mean.py                        \
        --input     $OCEAN_TRACER_ANOMALY_REFERENCE_FILE            \
        --basins    $PISM_BASIN_FILE                                \
        --fill      theta_ocean salinity_ocean                      \
        --output    $OCEAN_TRACER_ANOMALY_REFERENCE_FILE_BASIN_MEAN \
        --verbose
    RESULT=$?
    return_check $RESULT "calc_basin_contshelf_mean"

    echo
    echo "  > tracer anomaly (diff of ocean output to reference)"
    echo
    OCEAN_OUTPUT_TRACER_ANOMALY_FILE=$X_MOM_TO_PISM_PATH/$POEM_TIME_END.tracer.processed_MOM.anomaly.nc
    FILE_EXTENSION="${OCEAN_OUTPUT_TRACER_ANOMALY_FILE##*.}"
    FILE_PATH_NO_EXTENSION="${OCEAN_OUTPUT_TRACER_ANOMALY_FILE%.*}"
    OCEAN_TRACER_ANOMALY_FILE_BASIN_MEAN=$FILE_PATH_NO_EXTENSION.basin_mean.$FILE_EXTENSION

    $SCRIPT_DIR/calc_basin_contshelf_mean.py                        \
        --input     $OCEAN_OUTPUT_TRACER_ANOMALY_FILE  \
        --basins    $PISM_BASIN_FILE                                \
        --fill      theta_ocean salinity_ocean                      \
        --output    $OCEAN_TRACER_ANOMALY_FILE_BASIN_MEAN           \
        --verbose
    RESULT=$?
    return_check $RESULT "calc_basin_contshelf_mean"

    echo
    echo "  > tracer PISM baseline forcing"
    echo
    #PISM_OCN_FORCING_MOD_PATH=$FILE_PATH_NO_EXTENSION.mod.$FILE_EXTENSION
    FILE_EXTENSION="${PISM_OCN_FORCING_MOD_PATH##*.}"
    FILE_PATH_NO_EXTENSION="${PISM_OCN_FORCING_MOD_PATH%.*}"
    PISM_OCN_FORCING_MOD_BASIN_MEAN_PATH=$FILE_PATH_NO_EXTENSION.basin_mean.$FILE_EXTENSION

    $SCRIPT_DIR/calc_basin_contshelf_mean.py                        \
        --input     $PISM_OCN_FORCING_MOD_PATH                      \
        --basins    $PISM_BASIN_FILE                                \
        --fill      theta_ocean salinity_ocean                      \
        --output    $PISM_OCN_FORCING_MOD_BASIN_MEAN_PATH           \
        --verbose
    RESULT=$?
    return_check $RESULT "calc_basin_contshelf_mean"

    echo
    echo "  > tracer PISM input (anomalies + PISM baseline forcing)"
    echo
    PISM_INPUT_TRACER_FILE=$X_MOM_TO_PISM_PATH/$POEM_TIME_END.tracer.PISM_input.nc
    FILE_EXTENSION="${PISM_INPUT_TRACER_FILE##*.}"
    FILE_PATH_NO_EXTENSION="${PISM_INPUT_TRACER_FILE%.*}"
    PISM_INPUT_TRACER_FILE_BASIN_MEAN=$FILE_PATH_NO_EXTENSION.basin_mean.$FILE_EXTENSION

    $SCRIPT_DIR/calc_basin_contshelf_mean.py                        \
        --input     $PISM_INPUT_TRACER_FILE                         \
        --basins    $PISM_BASIN_FILE                                \
        --fill      theta_ocean salinity_ocean                      \
        --output    $PISM_INPUT_TRACER_FILE_BASIN_MEAN              \
        --verbose
    RESULT=$?
    return_check $RESULT "calc_basin_contshelf_mean"


    if [ $DO_OCEAN_SEALEVEL_ANOMALY = true ]; then
        echo
        echo "  > sea level reference (baseline for computing anomalies)"
        echo
        FILE_EXTENSION="${OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE##*.}"
        FILE_PATH_NO_EXTENSION="${OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE%.*}"
        OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE_BASIN_MEAN=$FILE_PATH_NO_EXTENSION.basin_mean.$FILE_EXTENSION

        $SCRIPT_DIR/calc_basin_contshelf_mean.py                            \
            --input     $OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE              \
            --basins    $PISM_BASIN_FILE                                    \
            --fill      delta_SL                                            \
            --output    $OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE_BASIN_MEAN   \
            --verbose
        RESULT=$?
        return_check $RESULT "calc_basin_contshelf_mean"

        echo
        echo "  > sea level anomaly (diff of ocean output to reference)"
        echo
        OCEAN_OUTPUT_SEALEVEL_ANOMALY_FILE=$X_MOM_TO_PISM_PATH/$POEM_TIME_END.sealevel.PISM_input.nc
        FILE_EXTENSION="${OCEAN_OUTPUT_SEALEVEL_ANOMALY_FILE##*.}"
        FILE_PATH_NO_EXTENSION="${OCEAN_OUTPUT_SEALEVEL_ANOMALY_FILE%.*}"
        OCEAN_SEALEVEL_ANOMALY_FILE_BASIN_MEAN=$FILE_PATH_NO_EXTENSION.basin_mean.$FILE_EXTENSION

        $SCRIPT_DIR/calc_basin_contshelf_mean.py                    \
            --input     $OCEAN_OUTPUT_SEALEVEL_ANOMALY_FILE         \
            --basins    $PISM_BASIN_FILE                            \
            --fill      delta_SL                                    \
            --output    $OCEAN_SEALEVEL_ANOMALY_FILE_BASIN_MEAN     \
            --verbose
        RESULT=$?
        return_check $RESULT "calc_basin_contshelf_mean"
    fi

    TIME_END_CALC_BASIN_CONTSHELF_MEAN=$(date +%s.%N)
    TIME_CALC_BASIN_CONTSHELF_MEAN=$(echo "$TIME_END_CALC_BASIN_CONTSHELF_MEAN - $TIME_START_CALC_BASIN_CONTSHELF_MEAN" | bc | awk '{printf "%f", $0}')


    #echo "  > copying basin continental shelf mean files to ocean history directory"

    #OCN_TO_ICE_OUTPUT_DIR=$POEM_OUTPUT_DIR/x_MOM-to-PISM/
    #mkdir -p $OCN_TO_ICE_OUTPUT_DIR

    #cp -av $OCEAN_OUTPUT_FILE_BASIN_MEAN $OCN_TO_ICE_OUTPUT_DIR
    #cp -av $OCEAN_TRACER_ANOMALY_REFERENCE_FILE_BASIN_MEAN $OCN_TO_ICE_OUTPUT_DIR
    #cp -av $OCEAN_TRACER_ANOMALY_FILE_BASIN_MEAN $OCN_TO_ICE_OUTPUT_DIR
    #cp -av $PISM_OCN_FORCING_MOD_BASIN_MEAN_PATH $OCN_TO_ICE_OUTPUT_DIR
    #cp -av $PISM_INPUT_TRACER_FILE_BASIN_MEAN $OCN_TO_ICE_OUTPUT_DIR
    #cp -av $OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE_BASIN_MEAN $OCN_TO_ICE_OUTPUT_DIR
    #cp -av $OCEAN_SEALEVEL_ANOMALY_FILE_BASIN_MEAN $OCN_TO_ICE_OUTPUT_DIR


    echo
    echo "  > cleaning up (move basin mean files)"
    echo

    OCN_TO_ICE_OUTPUT_DIR=$X_MOM_TO_PISM_PATH/
    mkdir -p $OCN_TO_ICE_OUTPUT_DIR/basin_mean
    #mkdir -p $OCN_TO_ICE_OUTPUT_DIR/tmp

    mv \
    $OCEAN_OUTPUT_FILE_BASIN_MEAN \
    $OCEAN_TRACER_ANOMALY_REFERENCE_FILE_BASIN_MEAN \
    $OCEAN_TRACER_ANOMALY_FILE_BASIN_MEAN \
    $PISM_OCN_FORCING_MOD_BASIN_MEAN_PATH \
    $PISM_INPUT_TRACER_FILE_BASIN_MEAN \
    $OCN_TO_ICE_OUTPUT_DIR/basin_mean

    if [ $DO_OCEAN_SEALEVEL_ANOMALY = true ]; then
        mv \
        $OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE_BASIN_MEAN \
        $OCEAN_SEALEVEL_ANOMALY_FILE_BASIN_MEAN \
        $OCN_TO_ICE_OUTPUT_DIR/basin_mean
    fi

    #mv $OCN_TO_ICE_OUTPUT_DIR/*.nc $OCN_TO_ICE_OUTPUT_DIR/tmp
    #mv $OCN_TO_ICE_OUTPUT_DIR/keep/* $OCN_TO_ICE_OUTPUT_DIR
    #rm -d $OCN_TO_ICE_OUTPUT_DIR/keep



}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function ncmax { ncap2 -O -C -v -s "foo=${1}.max();print(foo)" ${2} ~/tmp/foo.nc | cut -f 3- -d ' ' ; }
function ncmin { ncap2 -O -C -v -s "foo=${1}.min();print(foo)" ${2} ~/tmp/foo.nc | cut -f 3- -d ' ' ; }

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
return_check() {
    # $1 : RETURN_VAL
    # $2 : FUNC_STRING
    # usage: return_check $RETURN_VAL "func xy"

    if [ 0 = "$1" ] ; then
        echo "'"$2"'" ended successfully
    else
        echo "'"$2"'" exited with code $1
        exit_script $1
        #cd $ROOT_WORK_DIR
        #rm -f .workdir_locked_by_batchjob
        #exit $1
    fi
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
exit_script() {
    # $1 : RETURN_VAL
    # usage: exit_script $RESULT
    cd $ROOT_WORK_DIR
    rm -f .workdir_locked_by_batchjob
    exit $1
}

################################ start script ##################################

# load settings
source ./process_poem_to_pism.settings

export LANG=C
export LC_NUMERIC=en_US.UTF-8

# --- load packages ---
module purge
module load poem/2017.4

# load anaconda python environment for pre & inter-model processing scripts
module load anaconda
source activate py3_netcdf_xarray

# create folder to store intermediate and final data
mkdir -p $X_MOM_TO_PISM_PATH/
mkdir -p $X_ATM_TO_PISM_PATH/

#if [ $COMPUTE_OCEAN_TRACER_ANOMALY_FILE = true ]; then
mkdir -p $CALC_OCN_TRACER_ANOMALY_OUTPUT_PATH
#fi
#if [ $COMPUTE_OCEAN_SEALEVEL_ANOMALY_FILE = true ]; then
if [ $DO_OCEAN_SEALEVEL_ANOMALY = true ]; then
    mkdir -p $CALC_OCN_SEALEVEL_ANOMALY_OUTPUT_PATH
fi
#if [ $COMPUTE_ATMOS_ANOMALY_FILE = true ]; then
mkdir -p $CALC_ATM_ANOMALY_OUTPUT_PATH
#fi


# copy required scripts from coupling source directory
cp -av $CPL_SRC_DIR/inter-model-processing/regriddedMOM-to-PISM_processing.py $SCRIPT_DIR
cp -av $CPL_SRC_DIR/post-processing/calc_basin_contshelf_mean.py $SCRIPT_DIR

echo Hard Limits:
ulimit -aH
# unset some limits
#ulimit -c unlimited             # unlimit maximum size of core files created
ulimit -s unlimited             # unlimit maximum stack size
ulimit -d unlimited             # unlimit maximum size of a process's data segment
#ulimit -m unlimited             # unlimit maximum resident set size
ulimit -m `ulimit -H -m`
#ulimit -v unlimited             # unlimit maximum amount of virtual memory for shell
#ulimit -f unlimited             # unlimit maximum size of files written
ulimit -t unlimited             # unlimit maximum amount of cpu time (sec)

echo Soft limits:
ulimit -a
echo

set -x


# check whether sealevel anomaly should be calculated
case $DO_OCEAN_SEALEVEL_ANOMALY in
    true)
        echo DO_OCEAN_SEALEVEL_ANOMALY=$DO_OCEAN_SEALEVEL_ANOMALY
    ;;
    false)
        echo DO_OCEAN_SEALEVEL_ANOMALY=$DO_OCEAN_SEALEVEL_ANOMALY
    ;;
    *)
        echo DO_OCEAN_SEALEVEL_ANOMALY is supposed to be either "true" or "false"
        echo given value: DO_OCEAN_SEALEVEL_ANOMALY=$DO_OCEAN_SEALEVEL_ANOMALY
        exit -1
    ;;
esac

concat_poem_output
prepare_ocean_tracer_anomaly_reference_file
prepare_atmos_anomaly_reference_file
if [ $DO_OCEAN_SEALEVEL_ANOMALY = true ]; then
    prepare_ocean_sealevel_anomaly_reference_file
fi
process_mom_to_pism
process_atm_to_pism
calc_basin_contshelf_mean


### --------------------------- runtime statistics -----------------------------
set +x

print_stat() {
    # $1 : subroutine name
    # $2 : subroutine time (sec)
    # $3 : total script time (sec)

    format="  %-30s %10.2f \t %6.2f\n"
    TIME_PERCENT=$(echo "$2 / $3 * 100" | bc -l)
    printf "$format" $1 $2 $TIME_PERCENT
}

TIME_END_SCRIPT=$(date +%s.%N)
TIME_SCRIPT=$(echo "$TIME_END_SCRIPT - $TIME_START_SCRIPT" | bc | awk '{printf "%f", $0}')


str="-------------"
sep_str=$str$str$str$str$str
header="\n  %-30s %10s %11s \n"
width=54

echo
echo
echo
echo
echo " ---------------- runtime statistics ----------------- "
echo
echo " model time (yrs)          " $SIMULATION_TIME
echo " coupling timestep (yrs)   " $CPL_TIMESTEP
echo
echo
printf "$header" "ROUTINE" "TIME(s)" "ratio(%)"
printf "%$width.${width}s\n" " $sep_str"
print_stat total $TIME_SCRIPT $TIME_SCRIPT
print_stat concat_poem_output $TIME_CONCAT_POEM_OUTPUT $TIME_SCRIPT
print_stat ocn_tracer_anomaly_ref $TIME_OCN_TRACER_ANOMALY_REF $TIME_SCRIPT
print_stat ocn_sealevel_anomaly_ref $TIME_OCN_SL_ANOMALY_REF $TIME_SCRIPT
print_stat process_pism_to_mom $TIME_MOM_to_PISM_PROCESS $TIME_SCRIPT
print_stat process_atm_to_mom $TIME_ATM_to_PISM_PROCESS $TIME_SCRIPT
print_stat calc_basin_contshelf_mean $TIME_CALC_BASIN_CONTSHELF_MEAN $TIME_SCRIPT
printf "%$width.${width}s\n" " $sep_str"

echo
echo


set -x

exit_script $RESULT
