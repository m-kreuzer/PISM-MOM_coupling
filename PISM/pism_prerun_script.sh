#! /bin/sh

# execution of PISM for a single coupling iteration 
#
# this script is part of the coupled PISM-MOM setup and gets called by the 
# main coupling routine 'run_coupled.sh' in each coupling iteration
#
# all PISM options, parameters and flags are specified here except the 
# input/restart file $PISM_RESTART_FILE taken from in run_coupled.sh


echo SLURM_NTASKS:      $SLURM_NTASKS 
echo ROOT_WORK_DIR:     $ROOT_WORK_DIR
echo PISM_PROJ_DIR:     $PISM_PROJ_DIR 
echo PISM_WORK_DIR:     $PISM_WORK_DIR 
echo PISM_RESTART_FILE: $PISM_RESTART_FILE
echo __TIMESTEP:        $__TIMESTEP  
echo PISM_TIME_BEGIN:   $PISM_TIME_BEGIN
echo PISM_TIME_END:     $PISM_TIME_END 
echo $SLURM_NTASKS $ROOT_WORK_DIR $PISM_PROJ_DIR $PISM_WORK_DIR $PISM_RESTART_FILE $__TIMESTEP $PISM_TIME_BEGIN $PISM_TIME_END 
#echo SLURM_NTASKS $SLURM_NTASKS


srun -n $SLURM_NTASKS $PISM_PROJ_DIR/bin/pismr \
-i $PISM_RESTART_FILE \
-bootstrap -Mx $P_Mx -My $P_Mx -Lz $P_Lz -Lbz $P_Lbz -Mz $P_Mz -Mbz $P_Mbz \
-config_override $PISM_WORK_DIR/initdata/pism_config_override.nc \
-y $__TIMESTEP \
-verbose 2 \
-options_left \
-o_format netcdf4_parallel \
-atmosphere pik \
-atmosphere_pik era_interim \
-atmosphere_pik_temp_file initdata/racmo_wessem_initmip16km_mean1986_2005.nc \
-surface pdd \
-ocean pico \
-ocean_pico_file $PISM_OCEAN_START_FILE \
-o $PISM_WORK_DIR/prerun/prerun.pism_out.nc \
-extra_file $PISM_WORK_DIR/prerun/prerun.pism_extra.nc  \
-extra_times $PISM_TIME_BEGIN:$__TIMESTEP:$PISM_TIME_END \
-extra_vars basal_mass_flux_floating,basins,amount_fluxes,pdd_fluxes \
-save_file $PISM_WORK_DIR/prerun/prerun.pism_snap.nc \
-save_times $PISM_TIME_BEGIN:$__TIMESTEP:$PISM_TIME_END \
    > $PISM_WORK_DIR/prerun/prerun.pism.out 2>&1


RESULT=$?

exit $RESULT
