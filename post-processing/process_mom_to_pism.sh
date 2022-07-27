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
# standalone version of process_mom_to_pism() routine of PISM-MOM coupling
# framework
#
# This script calculates the PISM/PICO's input variables (temp, salt, eta) from
# ocean standalone output. Use case: compare PISM/PICO input variables in
# standalone and coupled mode.  
# The script requires some information for regridding and reference files,
# which are typically taken from coupled runs. They can be specified in
# process_mom_to_pism.settings
# -----------------------------------------------------------------------------

#SBATCH --account=poem
#SBATCH --job-name=process_mom_to_pism
#SBATCH --tasks=1
#SBATCH --mem=60000     # unit: MB
##SBATCH --qos=priority
#SBATCH --qos=short
#SBATCH --time=0:30:00
#SBATCH --output=process_mom_to_pism.%j.out
#SBATCH --error=process_mom_to_pism.%j.err
#SBATCH --mail-type=END,FAIL,REQUEUE,STAGE_OUT,TIME_LIMIT_90,TIME_LIMIT

TIME_START_SCRIPT=$(date +%s.%N)

# ---------------------------- function definitions ----------------------------
concat_mom_output(){
    echo
    echo
    echo " >> concat scattered outputs into single file"
    echo
    TIME_START_CONCAT_MOM_OUTPUT=$(date +%s.%N)
    
    cd $OCN_OUTPUT_DIR
    
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
        echo "  > concatenating MOM output files along time dimension"

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

    TIME_END_CONCAT_MOM_OUTPUT=$(date +%s.%N)
    TIME_CONCAT_MOM_OUTPUT=$(echo "$TIME_END_CONCAT_MOM_OUTPUT - $TIME_START_CONCAT_MOM_OUTPUT" | bc | awk '{printf "%f", $0}')
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
            echo "  > computing ocean tracer anomaly reference state from $CALC_OCN_TRACER_ANOMALY_PATH"
            PWD_SAVE=$PWD
            cd $CALC_OCN_TRACER_ANOMALY_PATH
            echo "   > average model output (possibly spread across multile files)"
            INPUT_FILES=$(echo `seq -f "$CALC_OCN_TRACER_ANOMALY_NAME_FORMAT_IN" \
                                       "$CALC_OCN_TRACER_ANOMALY_YR_START"  \
                                       "$CALC_OCN_TRACER_ANOMALY_YR_STEP"   \
                                       "$CALC_OCN_TRACER_ANOMALY_YR_END" `)
            AVG_OUT_FILE=$(echo `printf "$CALC_OCN_TRACER_ANOMALY_NAME_FORMAT_OUT" \
                                        "$CALC_OCN_TRACER_ANOMALY_YR_START"   \
                                        "$CALC_OCN_TRACER_ANOMALY_YR_END" `)  
            AVG_OUT_PATH=$ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$AVG_OUT_FILE
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

            REGRID_OUT=$ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/ocean_tracer_mean.regrid.MOM-to-PISM.bil.cdo.nc
            # remove global attribute _NCProperties from POEM output which causes cdo to crash
            ncatted -O -a _NCProperties,global,d,, $REGRID_IN
            cdo -b F64 -f nc4c remap,$PISM_PRE_OUT_FILE,$WEIGHTS_PATH $REGRID_IN $REGRID_OUT
            RESULT=$?
            return_check $RESULT "cdo_regridding"
            
            OUTFILE=$(echo `printf "%06g-%06g.tracer_mean.processed_MOM.nc" \
                                   "$CALC_OCN_TRACER_ANOMALY_YR_START"      \
                                   "$CALC_OCN_TRACER_ANOMALY_YR_END" `)  
            OCEAN_TRACER_ANOMALY_REFERENCE_FILE=$ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$OUTFILE
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
            echo "  > computing ocean sealevel anomaly reference state from $CALC_OCN_SEALEVEL_ANOMALY_PATH"
            PWD_SAVE=$PWD
            cd $CALC_OCN_SEALEVEL_ANOMALY_PATH
            echo "   > average model output (possibly spread across multile files)"
            INPUT_FILES=$(echo `seq -f "$CALC_OCN_SEALEVEL_ANOMALY_NAME_FORMAT_IN"  \
                                       "$CALC_OCN_SEALEVEL_ANOMALY_YR_START"        \
                                       "$CALC_OCN_SEALEVEL_ANOMALY_YR_STEP"         \
                                       "$CALC_OCN_SEALEVEL_ANOMALY_YR_END" `)
            AVG_OUT_FILE=$(echo `printf "$CALC_OCN_SEALEVEL_ANOMALY_NAME_FORMAT_OUT" \
                                        "$CALC_OCN_SEALEVEL_ANOMALY_YR_START"        \
                                        "$CALC_OCN_SEALEVEL_ANOMALY_YR_END" `)  
            AVG_OUT_PATH=$ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$AVG_OUT_FILE
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

            REGRID_OUT=$ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/ocean_sealevel_mean.regrid.MOM-to-PISM.bil.cdo.nc
            # remove global attribute _NCProperties from POEM output which causes cdo to crash
            ncatted -O -a _NCProperties,global,d,, $REGRID_IN
            cdo -b F64 -f nc4c remap,$PISM_PRE_OUT_FILE,$WEIGHTS_PATH $REGRID_IN $REGRID_OUT
            RESULT=$?
            return_check $RESULT "cdo_regridding"
            
            OUTFILE=$(echo `printf "%06g-%06g.sealevel_mean.processed_MOM.nc" \
                                   "$CALC_OCN_SEALEVEL_ANOMALY_YR_START"      \
                                   "$CALC_OCN_SEALEVEL_ANOMALY_YR_END" `)  
            OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE=$ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$OUTFILE
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
process_pism_to_mom(){
    echo
    echo
    echo " >> process MOM output for input in PISM"
    echo
    TIME_START_MOM_to_PISM_PROCESS=$(date +%s.%N)
    
    echo "  > bilinear regridding MOM output to PISM grid"

    # extract temp, salt, eta_t for regridding
    OCEAN_OUT=$OCN_OUTPUT_DIR/$POEM_TIME_END.$OCN_OUTPUT_FILE_BASE.nc
    REGRID_IN=$ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$POEM_TIME_END.$OCN_OUTPUT_FILE_BASE.sparse.nc
    ncks -O -v temp,salt,eta_t,time_bounds $OCEAN_OUT tmp.nc
    ncks -C -O -x -v geolat_t,geolon_t tmp.nc $REGRID_IN
    rm tmp.nc

    REGRID_OUT=$ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$POEM_TIME_END.regrid.MOM-to-PISM.bil.cdo.nc
    # remove global attribute _NCProperties from POEM output which causes cdo to crash
    ncatted -O -a _NCProperties,global,d,, $REGRID_IN
    cdo -b F64 -f nc4c remap,$PISM_PRE_OUT_FILE,$WEIGHTS_PATH $REGRID_IN $REGRID_OUT
    RESULT=$?
    return_check $RESULT "cdo_regridding"
    
    echo "  > regriddedMOM-to-PISM_processing script"
    $SCRIPT_DIR/regriddedMOM-to-PISM_processing.py                                   \
        --input     $REGRID_OUT                                                      \
        --basins    $PISM_BASIN_FILE                                                 \
        --edges     $PISM_EDGES_FILE                                                 \
        --fill      temp salt eta_t                                                  \
        --depth     $PICO_INPUT_DEPTH_FILE                                           \
        --output    $ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$POEM_TIME_END.processed_MOM.nc     \
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
        $ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$POEM_TIME_END.processed_MOM.nc \
        $OCEAN_TRACER_ANOMALY_REFERENCE_FILE_NOTIME \
        $ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$POEM_TIME_END.tracer.processed_MOM.anomaly.nc
    RESULT=$?
    return_check $RESULT "ncbo"
    # copy time and basins variables to PISM forcing file (omitted by ncbo add operation above)
    ncks -A -v time,time_bounds,basins \
        $ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$POEM_TIME_END.processed_MOM.nc \
        $ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$POEM_TIME_END.tracer.processed_MOM.anomaly.nc
    RESULT=$?
    return_check $RESULT "ncks"
    
    # add difference to given ocean forcing file (temp, salt)
    cp -av $PISM_OCN_FORCING $ROOT_WORK_DIR/evaluation/x_MOM-to-PISM
    PISM_OCN_FORCING_FILE=$(basename "$PISM_OCN_FORCING")
    PISM_OCN_FORCING_PATH=$ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$PISM_OCN_FORCING_FILE
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
    ncbo -O --op_typ=add \
        $ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$POEM_TIME_END.tracer.processed_MOM.anomaly.nc \
        $PISM_OCN_FORCING_MOD_PATH  \
        $ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$POEM_TIME_END.tracer.PISM_input.nc
    RESULT=$?
    return_check $RESULT "ncbo"
    # copy time and basins variables to PISM forcing file (omitted by ncbo add operation above)
    ncks -A -v time,time_bounds,basins \
        $ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$POEM_TIME_END.processed_MOM.nc \
        $ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$POEM_TIME_END.tracer.PISM_input.nc
    RESULT=$?
    return_check $RESULT "ncks"
    
    
    
    
    ### ocean sealevel anomaly approach: force PISM with sea level change
    #       from ocean model. Calculated through difference of current local 
    #       sea level (output of regriddedMOM-to-PISM_processing.py) to 
    #       reference file.
    echo OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE=$OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE
    echo OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE_NOTIME=$OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE_NOTIME
    
    # difference from current ocean output to reference (temp, salt)
    ncbo -O -v delta_SL,delta_SL_basin_mean --op_typ=subtract \
        $ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$POEM_TIME_END.processed_MOM.nc \
        $OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE_NOTIME \
        $ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$POEM_TIME_END.sealevel.processed_MOM.anomaly.nc
    RESULT=$?
    return_check $RESULT "ncbo"
    # copy time variable to PISM forcing file (omitted by ncbo add operation above)
    ncks -A -v time \
        $ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$POEM_TIME_END.processed_MOM.nc \
        $ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$POEM_TIME_END.sealevel.processed_MOM.anomaly.nc
    RESULT=$?
    return_check $RESULT "ncks"
    
    # set NaNs to 0
    ncatted -O -a _FillValue,delta_SL,m,f,0 \
        $ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$POEM_TIME_END.sealevel.processed_MOM.anomaly.nc
    RESULT=$?
    return_check $RESULT "ncatted"
    ncatted -O -a _FillValue,delta_SL,d,, \
        $ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$POEM_TIME_END.sealevel.processed_MOM.anomaly.nc
    RESULT=$?
    return_check $RESULT "ncatted"
    
    # create symlink from PISM input file to sealevel anomaly file
    PWD_SAVE=$PWD
    cd $ROOT_WORK_DIR/evaluation/x_MOM-to-PISM
    ln -sf $POEM_TIME_END.sealevel.processed_MOM.anomaly.nc $POEM_TIME_END.sealevel.PISM_input.nc
    cd $PWD_SAVE
    
    
    TIME_END_MOM_to_PISM_PROCESS=$(date +%s.%N)
    TIME_MOM_to_PISM_PROCESS=$(echo "$TIME_END_MOM_to_PISM_PROCESS - $TIME_START_MOM_to_PISM_PROCESS" | bc | awk '{printf "%f", $0}')
}


calc_basin_contshelf_mean(){
    echo
    echo
    echo " >> calculate continental shelf mean values for each basin"
    echo

    TIME_START_CALC_BASIN_CONTSHELF_MEAN=$(date +%s.%N)

    echo "  > ocean output (processed to PISM grid)"
    echo
    OCEAN_OUTPUT_FILE=$ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$POEM_TIME_END.processed_MOM.nc
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

    echo "  > tracer anomaly (diff of ocean output to reference)"
    echo
    OCEAN_OUTPUT_TRACER_ANOMALY_FILE=$ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$POEM_TIME_END.tracer.processed_MOM.anomaly.nc
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

    ### tracer PISM input (anomalies + PISM baseline forcing)
    PISM_INPUT_TRACER_FILE=$ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$POEM_TIME_END.tracer.PISM_input.nc
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

    echo "  > sea level anomaly (diff of ocean output to reference)"
    echo
    OCEAN_OUTPUT_SEALEVEL_ANOMALY_FILE=$ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/$POEM_TIME_END.sealevel.PISM_input.nc
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

    TIME_END_CALC_BASIN_CONTSHELF_MEAN=$(date +%s.%N)
    TIME_CALC_BASIN_CONTSHELF_MEAN=$(echo "$TIME_END_CALC_BASIN_CONTSHELF_MEAN - $TIME_START_CALC_BASIN_CONTSHELF_MEAN" | bc | awk '{printf "%f", $0}')


    #echo "  > copying basin continental shelf mean files to ocean history directory"

    #OCN_TO_ICE_OUTPUT_DIR=$OCN_OUTPUT_DIR/x_MOM-to-PISM/
    #mkdir -p $OCN_TO_ICE_OUTPUT_DIR

    #cp -av $OCEAN_OUTPUT_FILE_BASIN_MEAN $OCN_TO_ICE_OUTPUT_DIR
    #cp -av $OCEAN_TRACER_ANOMALY_REFERENCE_FILE_BASIN_MEAN $OCN_TO_ICE_OUTPUT_DIR
    #cp -av $OCEAN_TRACER_ANOMALY_FILE_BASIN_MEAN $OCN_TO_ICE_OUTPUT_DIR
    #cp -av $PISM_OCN_FORCING_MOD_BASIN_MEAN_PATH $OCN_TO_ICE_OUTPUT_DIR
    #cp -av $PISM_INPUT_TRACER_FILE_BASIN_MEAN $OCN_TO_ICE_OUTPUT_DIR
    #cp -av $OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE_BASIN_MEAN $OCN_TO_ICE_OUTPUT_DIR
    #cp -av $OCEAN_SEALEVEL_ANOMALY_FILE_BASIN_MEAN $OCN_TO_ICE_OUTPUT_DIR


    echo "  > cleaning up intermediate files (move to tmp/ dir)"

    OCN_TO_ICE_OUTPUT_DIR=$ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/
    mkdir -p $OCN_TO_ICE_OUTPUT_DIR/keep
    mkdir -p $OCN_TO_ICE_OUTPUT_DIR/tmp

    mv \
    $OCEAN_OUTPUT_FILE_BASIN_MEAN \
    $OCEAN_TRACER_ANOMALY_REFERENCE_FILE_BASIN_MEAN \
    $OCEAN_TRACER_ANOMALY_FILE_BASIN_MEAN \
    $PISM_OCN_FORCING_MOD_BASIN_MEAN_PATH \
    $PISM_INPUT_TRACER_FILE_BASIN_MEAN \
    $OCEAN_SEALEVEL_ANOMALY_REFERENCE_FILE_BASIN_MEAN \
    $OCEAN_SEALEVEL_ANOMALY_FILE_BASIN_MEAN \
    $OCN_TO_ICE_OUTPUT_DIR/keep
    
    mv $OCN_TO_ICE_OUTPUT_DIR/*.nc $OCN_TO_ICE_OUTPUT_DIR/tmp
    mv $OCN_TO_ICE_OUTPUT_DIR/keep/* $OCN_TO_ICE_OUTPUT_DIR
    rm -d $OCN_TO_ICE_OUTPUT_DIR/keep


}
    
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
source ./process_mom_to_pism.settings

export LANG=C
export LC_NUMERIC=en_US.UTF-8

# --- load packages ---
module purge
module load poem/2017.4

# load anaconda python environment for pre & inter-model processing scripts
module load anaconda
source activate py3_netcdf_xarray

# create folder to store intermediate and final data
mkdir -p $ROOT_WORK_DIR/evaluation/x_MOM-to-PISM/


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

concat_mom_output
prepare_ocean_tracer_anomaly_reference_file
prepare_ocean_sealevel_anomaly_reference_file
process_pism_to_mom
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
print_stat concat_mom_output $TIME_CONCAT_MOM_OUTPUT $TIME_SCRIPT
print_stat ocn_tracer_anomaly_ref $TIME_OCN_TRACER_ANOMALY_REF $TIME_SCRIPT
print_stat ocn_sealevel_anomaly_ref $TIME_OCN_SL_ANOMALY_REF $TIME_SCRIPT
print_stat process_pism_to_mom $TIME_MOM_to_PISM_PROCESS $TIME_SCRIPT
print_stat calc_basin_contshelf_mean $TIME_CALC_BASIN_CONTSHELF_MEAN $TIME_SCRIPT
printf "%$width.${width}s\n" " $sep_str"

echo
echo


set -x

exit_script $RESULT
