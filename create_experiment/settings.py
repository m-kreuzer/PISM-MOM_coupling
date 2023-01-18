###
# Set basic settings here,e.g., adjust path to input files, 
# select the run type, ensemble parameters...
# runscripts are created via create_run.py or create_set.py
# more options are set in templates/pism_run.sh.jinja2

import os
import pwd
import collections

import pism_grids

# import settings including path to input, output directories
from pikcluster_settings import *


# ----------------------------- coupling settings ------------------------------

experiment = "MOM5_PISM1.0hash_16km_1pctCO2_CCSM4_run01"
experiment_dir      = os.path.join(working_dir, experiment)


coupling_timestep = 5       # in years, must be greater or equal 1
max_cpl_iteration = 27      # number of coupling iterations


# - - - - - - - - - ice -> ocean fluxes [1st coupling iteration] - - - - - - - -

# in case of no coupled restart:
# specify ice-to-ocean runoff for first coupling iteration as ocean runs before
# ice and default setup for ocean has no Antarctic runoff from land/ice
pism_to_mom_flux_init_file = 'equi_16km_110000yrs.last_1ka.mean.fluxes.nc'
pism_to_mom_flux_init_path = os.path.join('/p/tmp/kreuzer/coupled_PISM_MOM/experiments/pism1.1_equi_16km_100000_plus_run05/output_processed/mom5.2023-01/', pism_to_mom_flux_init_file)

# - - - - - - - - - - - - - - - - - - restart - - - - - - - - - - - - - - - - -
# option to restart a coupled setup from a previous run
#  -> gives the opportunity to add PISM input fluxes to MOM from last coupling 
#     iteration of the previous run
#  -> still requires to specify PISM restart files by hand (see pism_infile 
#     and pism_infile_dir below)
#  -> MOM restart files are automatically copied from $restart_dir/POEM/INPUT/
#     if path exists
coupled_restart = False
restart_dir = "/p/tmp/kreuzer/coupled_PISM_MOM/experiments/coupling_framework_basal_melt_dev04"

### specify files computed by PISM-to-MOM_processing.py from PISM output in
#   last coupling iteration of the previous run

# PISM to MOM fluxes, inserted in MOM's data_table
pism_to_mom_flux_restart_file = '002994.fluxes.nc'
pism_to_mom_flux_restart_path = os.path.join(restart_dir, 'x_PISM-to-MOM', pism_to_mom_flux_restart_file)

# runoff reference, inserted in MOM's data_table
# -> only used when (do_runoff_slc==True AND runoff_reference_surf_accum==True)
runoff_reference_restart_file = '002994.runoff_reference.nc'
runoff_reference_restart_path = os.path.join(restart_dir, 'x_PISM-to-MOM', runoff_reference_restart_file)

# basal melt input depth, inserted in MOM's data_table
# -> only used when (insert_basal_melt_at_depth==True AND use_prescribed_basal_melt_input_depth==False)
basal_melt_input_depth_restart_file = '002994.basal_melt_input_depth.nc'
basal_melt_input_depth_restart_path = os.path.join(restart_dir, 'x_PISM-to-MOM', basal_melt_input_depth_restart_file)

# PICO input depth, used in regriddedMOM-to-PISM.py
# -> only used when (use_prescribed_pico_input_depth==False)
pico_input_depth_restart_file = '002994.pico_input_depth.nc'
pico_input_depth_restart_path = os.path.join(restart_dir, 'x_PISM-to-MOM', pico_input_depth_restart_file)



# - - - - - - - - - - - - - - MOM -> PISM timeseries - - - - - - - - - - - - - - 
# option to pass a time series of ocean forcing to PISM when coupling time step
# is greater than a year. Otherwise the ocean forcing is averaged over the 
# coupling time step before being processed and passed to PISM.
ocean_to_ice_timeseries = False

# - - - - - - - - - - - - - - ice sheet runoff SLC - - - - - - - - - - - - - - -
# In ocean/sea ice only configuration (no atm), the ocean model uses a
# normalisation of surface mass fluxes (precip-evap+river=0) to keep the mass
# of the ocean and sea ice system constant.
# (In MOM5 this is done via the `zero_net_water_coupler` flag in the `ocean_sbc`
# namelist).
# To avoid that changes in the ice sheet mass (and subsequently ice-to-ocean
# runoff) are being balanced out when provided to the ocean via the river routing
# scheme, the ice sheet runoff can be split:
#  (1) reference mass flux to ocean in steady state/ ice sheet equilibrium
#       (e.g. average at end of ice sheet standalone spinup)
#  (2) anomaly/difference to (1), the reference runoff
# While (1) is provided to MOM via the regular river runoff mechanism, (2) is
# given to MOM5 via the new implemented `runoff_slc` data table entry (which
# adds it to the same runoff field, but excludes it from the pme+river=0
# normalisation). With the `do_runoff_slc` flag this split can be enabled. It might
# be reasonable to switch it off during spinup phase of coupled ice/ocean
# simulations. 
# There are two possiblities to choose the ice to ocean reference runoff:
#  (a) compute it from PISM's surface accumulation
#       In ice sheet equilibrium, the surface accumulation flux should balance
#       the mass flux from ice to ocean. This method is recommended especially
#       with transient surface forcing on the ice sheet. To use this option,
#       set `runoff_reference_surf_accum = True`
#  (b) provide a file with the reference ice to ocean runoff This can be
#       computed manually as an average ice to ocean mass flux from an
#       standalone ice spinup/equilibrium state. 
#       WARNING: As the runoff (especially the calving) shows quasi-stochastic
#       behaviour, it is very likely that the computed mean does not exactly
#       represent the ice to ocean mean runoff. This introduces a systematic
#       positive/negative bias when computing the anomalies, which piles up in
#       a constant mass drift of the coupled ice-ocean system.  This drift can
#       be reduced, using a manual correction of the reference while on basis
#       of the coupled mass drift. 
#       This option is used when `runoff_reference_surf_accum = False`. In this
#       case the variables `runoff_reference_file` (one timestamp only!) and
#       `runoff_reference_path` are used. 
do_runoff_slc = True
runoff_reference_surf_accum = True
# use reference file if (do_runoff_slc==True & runoff_reference_surf_accum==False) 
runoff_reference_file = "equi_16km_110000yrs.mean_last_1ka.fluxes.nc"
runoff_reference_path = os.path.join("/p/tmp/kreuzer/coupled_PISM_MOM/experiments/pism1.1_equi_16km_100000_plus_run03/output_processed", runoff_reference_file)


# - - - - - - - - - - - - - basal melt insertion depth - - - - - - - - - - - - -
# MOM5/6 by default offers no possibility to insert freshwater other than via
# river runoff at surface. With a modified version, it is now also possible to
# insert freswater at depth.
insert_basal_melt_at_depth = True
use_prescribed_basal_melt_input_depth = False
prescribed_basal_melt_input_depth_file = 'basin_melt_input_depth.nc'
prescribed_basal_melt_input_depth_path = os.path.join('template_path', prescribed_basal_melt_input_depth_file)

# - - - - - - - - - - - - - - - ocean tracer anomaly - - - - - - - - - - - - - -
do_ocean_tracer_anomaly    = True

# in case of coupled restart with do_ocean_tracer_anomaly:
#     specify ocean tracer anomaly reference file from previous run
use_ocean_tracer_anomaly_from_prev_run = True
ocean_tracer_anomaly_reference_file = "017090-017090.tracer_mean.processed_MOM.nc"
#ocean_tracer_anomaly_reference_path = os.path.join(restart_dir, 'x_MOM-to-PISM', ocean_tracer_anomaly_reference_file)
ocean_tracer_anomaly_reference_path = os.path.join('/p/tmp/kreuzer/coupled_PISM_MOM/experiments/MOM5_standalone_EM3_spinup_PISM_runoff_run13/evaluation/x_MOM-to-PISM/tmp/', ocean_tracer_anomaly_reference_file)
#     or specify MOM output files used for computing ocean tracer anomaly reference state
#       -> used if do_ocean_tracer_anomaly==True and use_ocean_tracer_anomaly_from_prev_run==False
calc_ocn_tracer_anomaly = {}
calc_ocn_tracer_anomaly['path'] = "/p/tmp/kreuzer/coupled_PISM_MOM/experiments/MOM5_standalone_EM3_spinup_PISM_runoff_run13/history"
calc_ocn_tracer_anomaly['yr_start'] = "17090"
calc_ocn_tracer_anomaly['yr_end'] = "17090"
calc_ocn_tracer_anomaly['yr_step'] = "10"
calc_ocn_tracer_anomaly['name_format_in']  = "%06g0101.ocean-decadal.nc"
#calc_ocn_tracer_anomaly['name_format_in']  = "%06g.ocean-decadal.tracer_mean.nc.bak"
calc_ocn_tracer_anomaly['name_format_out'] = "%06g-%06g.ocean-decadal.tracer_mean.nc"

# - - - - - - - - - - - - - - ocean sealevel anomaly - - - - - - - - - - - - - -
do_ocean_sealevel_anomaly    = False

# in case of coupled restart with do_ocean_sealevel_anomaly:
#     specify ocean sealevel anomaly reference file from previous run
use_ocean_sealevel_anomaly_from_prev_run = False
ocean_sealevel_anomaly_reference_file = "017090-017090.sealevel_mean.processed_MOM.nc"
ocean_sealevel_anomaly_reference_path = os.path.join(restart_dir, 'x_MOM-to-PISM', ocean_sealevel_anomaly_reference_file)
#     or specify MOM output files used for computing ocean sealevel anomaly reference state
#       -> used if do_ocean_sealevel_anomaly==True and use_ocean_sealevel_anomaly_from_prev_run==False
calc_ocn_sealevel_anomaly = {}
calc_ocn_sealevel_anomaly['path'] = "/p/tmp/kreuzer/coupled_PISM_MOM/experiments/MOM5_standalone_EM3_spinup_PISM_runoff_run13/history"
calc_ocn_sealevel_anomaly['yr_start'] = "17090"
calc_ocn_sealevel_anomaly['yr_end'] = "17090"
calc_ocn_sealevel_anomaly['yr_step'] = "10"
calc_ocn_sealevel_anomaly['name_format_in']  = "%06g0101.ocean-decadal.nc"
#calc_ocn_sealevel_anomaly['name_format_in']  = "%06g.ocean-decadal.sealevel_mean.nc.bak"
calc_ocn_sealevel_anomaly['name_format_out'] = "%06g-%06g.ocean-decadal.sealevel_mean.nc"

# - - - - - - - - - - - - - prescribed PICO input depth - - - - - - - - - - - -
# option to use static PICO input depth values read from file
#  -> gives possibility to use prescribed and fixed values for PICO input
#  -> depths are used by regriddedMOM-to-PISM_processing.py to select PICO
#     temperatures and salinities from vertical profile on a basin averaged
#     basins
#  -> if option is not used then basin shelf depth is calculated each
#     coupling time step in PISM-to-MOM_processing.py

use_prescribed_pico_input_depth = False
prescribed_pico_input_depth_file = 'basin_shelf_depth_oceanic_gateways.nc'
prescribed_pico_input_depth_path = os.path.join('/p/projects/climber3/kreuzer/POEM/POEM_PISM_coupling_templates/pico_input_depth', prescribed_pico_input_depth_file)



# ------------------------------- POEM settings --------------------------------

poem_exp_dir        = os.path.join(experiment_dir, 'POEM')

# - - - - - - - - - - - - - - - - POEM restart - - - - - - - - - - - - - - - - -
### in case of no coupled restart, select directory to copy POEM restart 
#   files from
poem_restart_files_dir = '/p/tmp/kreuzer/coupled_PISM_MOM/experiments/MOM5_standalone_EM3_spinup_PISM_runoff_run13_plus_1month/RESTART'


# - - - - - - - - - - - - - - - - POEM forcing - - - - - - - - - - - - - - - - - 
### select data_table template
#   no additional forcing:            data_table-dummy_mom5-clim
#   CMIP5 1pctCO2 forcing:            data_table-dummy_mom5-clim+CMIP5-anom
#   CMIP5 1pctCO2 extension forcing:  data_table-dummy_mom5-clim+CMIP5-anom

poem_data_table_dummy = 'data_table-dummy_mom5-clim+CMIP5-anom'
poem_data_table_replace = {'MODEL':            'CCSM4',
                           'SCENARIO':         '1pctCO2',
                           'CLIMEXTENSION':    ''} 
poem_copy_forcing_data = True
poem_forcing_data_source_dir = '/p/projects/pism/kreuzer/coupled_PISM_MOM/forcing/CMIP5_reanalysis/Amon_CCSM4_1pctCO2_r1i1p1_month_anomaly_plus_MOM5_climatology'
poem_forcing_data_source_pattern = '*.nc'
poem_forcing_data_source_path = os.path.join(poem_forcing_data_source_dir, poem_forcing_data_source_pattern)
poem_forcing_data_target_path = os.path.join(experiment_dir, 'POEM','INPUT','CMIP5_forcing')
poem_forcing_time_shift_years = (17090 - 1)

#poem_data_table_replace = {'MODEL':            'CCSM4',
#                           'SCENARIO':         '1pctCO2',
#                           'CLIMEXTENSION':    'ext-clim_',
#                           '\.shift':          ''} 
#poem_copy_forcing_data = True
#poem_forcing_data_source_dir = '/p/projects/pism/kreuzer/coupled_PISM_MOM/forcing/CMIP5_reanalysis/Amon_CCSM4_1pctCO2_r1i1p1_ext-clim_month_anomaly_plus_MOM5_climatology'
#poem_forcing_data_source_pattern = '*.nc'
#poem_forcing_data_source_path = os.path.join(poem_forcing_data_source_dir, poem_forcing_data_source_pattern)
#poem_forcing_data_target_path = os.path.join(experiment_dir, 'POEM','INPUT','CMIP5_forcing')
##poem_forcing_time_shift_years = 17090


# ------------------------------- PISM settings --------------------------------

# select resolution of the run
grid_id = "initmip16km" 
#grid_id = "initmip8km" 
pism_grid = pism_grids.grids[grid_id]

# directories and path definitions
pism_exp_dir        = os.path.join(experiment_dir, 'PISM')
pism_exp_bin_dir    = os.path.join(experiment_dir, 'PISM', 'bin')
pism_exp_bin        = os.path.join(pism_exp_bin_dir, pism_exec)
#pism_exp_bin        = os.path.join(settings.pism_code_dir, 'bin', settings.pism_exec)
pism_sys_bin        = os.path.join(pism_code_dir, 'bin', pism_exec)

# input data
pism_infile_dir = "/p/tmp/kreuzer/coupled_PISM_MOM/experiments/pism1.0_precipscale_hash_equi_16km_run10/output_pism/"
pism_infile = "result_equi_16km_50000.nc"
pism_infile_path = os.path.join(pism_infile_dir, pism_infile)

#pism_atm_data_dir = os.path.join(pism_input_root_dir, "racmo_wessem")
pism_atm_data_dir = pism_input_root_dir
pism_atm_file = "racmo_wessem_"+grid_id+"_mean1986_2005.nc"
#pism_atm_data_dir = os.path.join(pism_input_root_dir, "merged")
#pism_atm_file = "bedmap2_albmap_racmo_wessem_tillphi_pism_"+grid_id+".nc"
pism_atm_data_path = os.path.join(pism_atm_data_dir,pism_atm_file)

### 1pctCO2 atm anomaly forcing
pism_use_atm_anomaly_file = True
pism_atm_anomaly_data_dir = "/p/tmp/kreuzer/coupled_PISM_MOM/experiments/pism1.0_precipscale_hash_q0.625_16km_1pctCO2_CCSM4_yearly_run02/initdata/"
pism_atm_anomaly_file = "pdd_Amon_CCSM4_1pctCO2_r1i1p1_yearly_anomaly_relativeprecip_initmip16km.timeshift.nc"
pism_atm_anomaly_data_path = os.path.join(pism_atm_anomaly_data_dir, pism_atm_anomaly_file)
pism_atm_anomaly_time_shift_years = (17090 - 50001) 

### 1pctCO2ext atm anomaly forcing
#pism_use_atm_anomaly_file = True
#pism_atm_anomaly_data_dir = "/p/tmp/kreuzer/coupled_PISM_MOM/experiments/pism1.0_precipscale_hash_q0.625_16km_1pctCO2ext_CCSM4_yearly_run02/initdata/"
#pism_atm_anomaly_file = "pdd_Amon_CCSM4_1pctCO2_r1i1p1_ext-clim_yearly_anomaly_relativeprecip_initmip16km.timeshift.nc"
#pism_atm_anomaly_data_path = os.path.join(pism_atm_anomaly_data_dir, pism_atm_anomaly_file)
##pism_atm_anomaly_time_shift_years = (17090 - 50001) 

pism_use_atm_lapse_rate_file = True
#pism_atm_lapse_rate_data_dir = os.path.join(pism_input_root_dir, "merged")
pism_atm_lapse_rate_data_dir = pism_input_root_dir
pism_atm_lapse_rate_file = "bedmap2_albmap_racmo_wessem_tillphi_pism_"+grid_id+".nc"
pism_atm_lapse_rate_data_path = os.path.join(pism_atm_lapse_rate_data_dir,pism_atm_lapse_rate_file)

#pism_ocn_data_dir = os.path.join(pism_input_root_dir, "schmidtko")
#pism_ocn_file = "schmidtko_"+grid_id+"_means.nc"
pism_ocn_data_dir = pism_input_root_dir
#pism_ocn_file = "schmidtko_"+grid_id+".nc"
pism_ocn_file = "oceanWarming_schmidtko_0.0K.nc"
pism_ocn_data_path = os.path.join(pism_ocn_data_dir,pism_ocn_file)

#pism_ocnkill_data_dir = os.path.join(pism_input_root_dir, "bedmap2")
#pism_ocnkill_file = "bedmap2_"+grid_id+".nc"
#pism_ocnkill_data_dir = os.path.join(pism_input_root_dir, "merged")
pism_ocnkill_data_dir = pism_input_root_dir
pism_ocnkill_file = "bedmap2_albmap_racmo_wessem_tillphi_pism_"+grid_id+".nc"
pism_ocnkill_data_path = os.path.join(pism_ocnkill_data_dir,pism_ocnkill_file)




# - - - - - - - - - - - - - - - - config file - - - - - - - - - - - - - - -  -

# set pism parameters that apply to all runs (unless part of the ensemble)
pism_config_file = os.path.join(pism_code_dir,"src/pism_config.cdl")

# override parameters that deviate from default.
override_params = collections.OrderedDict([

("atmosphere.precip_exponential_factor_for_temperature", 0.04879),
("atmosphere.precip_exponential_factor_for_temperature_doc", "= ln(1.05); a 5% change of precipitation rate for every one degC of temperature change"),
("atmosphere.precip_exponential_factor_for_temperature_type", "scalar"),
("atmosphere.precip_exponential_factor_for_temperature_units", "Kelvin-1"),

("atmosphere.precip_lapse_scaling","yes"),
("atmosphere.precip_lapse_scale_factor", 0.41),
("atmosphere.precip_lapse_scale_factor_doc", "=8.2 K km-1 * 5 % K-1; atmospheric lapse rate * precipitation change per degree of warming; Scale precipitation according to change in surface elevation"),
("atmosphere.precip_lapse_scale_factor_option", "precip_scale_factor"),
("atmosphere.precip_lapse_scale_factor_type", "scalar"),
("atmosphere.precip_lapse_scale_factor_units", "km-1"),


("basal_yield_stress.model",  "mohr_coulomb"),
("basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden", 0.04),
("basal_resistance.pseudo_plastic.enabled",  "true"),
("basal_resistance.pseudo_plastic.q",  0.625),
("basal_resistance.pseudo_plastic.u_threshold", 100.0),

("bed_deformation.model", "lc"),

("calving.eigen_calving.K", 1.0e17), 
("calving.methods", "eigen_calving,thickness_calving,ocean_kill"), 
("calving.thickness_calving.threshold", 50),

("energy.basal_melt.use_grounded_cell_fraction", "false"),  # -no_sugl_basal_melt

("geometry.part_grid.enabled", "true"), # -pik
("geometry.remove_icebergs", "true"), # -pik
("geometry.grounded_cell_fraction", "true"), # -pik

("hydrology.model", "null"),

("ocean.pico.heat_exchange_coefficent", 1.0e-5),
("ocean.pico.overturning_coefficent", 0.5e6),
("ocean.pico.exclude_ice_rises", "yes"),
("ocean.pico.continental_shelf_depth", -2000),

("stress_balance.calving_front_stress_bc", "true"), # -pik
("stress_balance.model","ssa+sia"),
("stress_balance.sia.enhancement_factor",1.0),
("stress_balance.ssa.enhancement_factor",1.0),
("stress_balance.ssa.method", "fd"),

("surface.pdd.factor_snow", 0.0032967032967033),
("surface.pdd.factor_ice", 0.00879120879120879),
("surface.pdd.std_dev", 5),

])

# - - - - - - - - - - - - - - - command line options - - - - - - - - - - - - - -
pism_general_opt = "-verbose 2 -options_left -o_format netcdf4_parallel"
#pism_atm_opt = "-atmosphere pik -atmosphere_pik_file initdata/"+pism_atm_file+" -surface pdd"
pism_atm_opt = "-atmosphere pik_temp,anomaly,lapse_rate -atmosphere_pik_temp_file initdata/"+pism_atm_file+" -atmosphere_anomaly_file initdata/"+pism_atm_anomaly_file+" -atmosphere_lapse_rate_file initdata/"+pism_atm_lapse_rate_file+" -temp_era_interim -temp_lapse_rate 0.0 -precip_scale_factor .410 -surface pdd"
pism_add_opt = "-ocean_kill_file initdata/"+pism_ocnkill_file



# - - - - - - - - - - - - - - - - - - output - - - - - - - - - - - - - - - - -

# set extra variables
#pism_extra_vars = "mask,thk,velsurf_mag,velbase_mag,flux_mag,tillwat,tauc,pico_overturning,pico_temperature_box0,pico_salinity_box0,tillphi,shelfbmassflux,shelfbtemp,basins,bmelt,bfrict,bfrict,tendency_of_ice_amount,amount_fluxes,ice_mass,pico_contshelf_mask,pdd_fluxes,pdd_rates,climatic_mass_balance,pico_shelf_mask,pico_box_mask"
pism_extra_vars = "mask,thk,velsurf_mag,velbase_mag,flux_mag,tillwat,tauc,pico_overturning,pico_temperature_box0,pico_salinity_box0,tillphi,shelfbmassflux,shelfbtemp,basins,bmelt,bfrict,bfrict,tendency_of_ice_amount,amount_fluxes,ice_mass,pico_contshelf_mask,pdd_fluxes,pdd_rates,climatic_mass_balance,pico_shelf_mask,pico_box_mask,velbar_mag,velbar,surface_melt_rate,topg,dHdt,dbdt,usurf,ice_surface_temp,air_temp_mean_summer"

pism_ts_vars = "ice_volume,ice_volume_glacierized_temperate,ice_volume_glacierized_grounded,ice_volume_glacierized_floating,ice_volume_glacierized_cold,ice_volume_glacierized,ice_mass_glacierized,ice_enthalpy_glacierized,ice_area_glacierized_temperate_base,ice_area_glacierized_grounded,ice_area_glacierized_floating,ice_area_glacierized_cold_base,ice_area_glacierized,tendency_of_ice_volume_glacierized,tendency_of_ice_mass_glacierized,slvol,tendency_of_ice_mass_due_to_discharge,max_hor_vel,max_diffusivity,dt,basal_mass_flux_floating,basal_mass_flux_grounded,tendency_of_ice_mass_due_to_basal_mass_flux,tendency_of_ice_mass_due_to_flow,tendency_of_ice_mass_due_to_surface_mass_flux"

# set pism diagnostic timesteps
pism_diag_extra_timestep = coupling_timestep
pism_diag_snap_timestep = coupling_timestep 
pism_diag_ts_timestep = coupling_timestep 




# ---------------------- no edits below this line needed. ----------------------
project_root = os.path.dirname(os.path.abspath(__file__))
user = pwd.getpwuid(os.getuid()).pw_name
# ------------------------------------------------------------------------------ 




