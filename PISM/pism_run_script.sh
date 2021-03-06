#! /bin/sh

# execution of PISM for a single coupling iteration 
#
# this script is part of the coupled PISM-MOM setup and gets called by the 
# main coupling routine 'run_coupled.sh' in each coupling iteration
#
# all PISM options, parameters and flags are specified here except the 
# input/restart file $PISM_RESTART_FILE taken from in run_coupled.sh


srun -n $SLURM_NTASKS $PISM_PROJ_DIR/bin/pismr \
-i $PISM_RESTART_FILE \
-config_override ./initdata/pism_config_override.nc \
-y $__TIMESTEP \
-verbose 2 \
-options_left \
-o_format netcdf4_parallel \
-atmosphere pik \
-atmosphere_pik era_interim \
-atmosphere_pik_temp_file initdata/racmo_wessem_initmip16km_mean1986_2005.nc \
-surface pdd \
-ocean pico \
-ocean_pico_file $ROOT_WORK_DIR/x_MOM-to-PISM/$POEM_TIME_END.PISM_input.nc  \
-gamma_T 1.0e-5 \
-overturning_coeff 0.5e6 \
-exclude_icerises \
-continental_shelf_depth -2000 \
-bed_def lc \
-hydrology null \
-calving eigen_calving,ocean_kill \
-eigen_calving_K 1.0e17  \
-ocean_kill_file initdata/bedmap2_albmap_racmo_wessem_tillphi_pism_initmip16km.nc \
-pik \
-sia_e 1.0 \
-ssa_e 1.0 \
-ssa_method fd \
-stress_balance ssa+sia \
-sia_flow_law gpbld \
-ssa_flow_law gpbld \
-pseudo_plastic \
-pseudo_plastic_q 0.75 \
-pseudo_plastic_uthreshold 100.0 \
-till_effective_fraction_overburden 0.04 \
-subgl \
-no_subgl_basal_melt \
-o results/$PISM_TIME_END.pism_out.nc \
-o_size big \
-log_view \
-extra_file results/$PISM_TIME_END.pism_extra.nc \
-extra_times $PISM_TIME_BEGIN:$__TIMESTEP:$PISM_TIME_END \
-extra_vars basal_mass_flux_floating,basal_mass_flux_grounded,basins,pico_overturning,pico_salinity_box0,pico_temperature_box0,pico_box_mask,pico_shelf_mask,pico_ice_rise_mask,pico_basal_melt_rate,pico_contshelf_mask,pico_salinity,pico_temperature,pico_T_star,pico_basal_temperature,amount_fluxes,pdd_fluxes,ice_mass,enthalpy \
-save_file results/$PISM_TIME_END.pism_snap.nc \
-save_times $PISM_TIME_BEGIN:$__TIMESTEP:$PISM_TIME_END \
-ts_file results/$PISM_TIME_END.pism_ts.nc \
-ts_times $PISM_TIME_BEGIN:$__TIMESTEP:$PISM_TIME_END \
    > results/$PISM_TIME_END.pism.out 2>&1

RESULT=$?

exit $RESULT
