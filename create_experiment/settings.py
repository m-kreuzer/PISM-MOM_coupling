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

experiment = "MOM5_PISM_16km_gmd-2020-230_run04_post"
experiment_dir      = os.path.join(working_dir, experiment)


coupling_timestep = 10      # in years, must be greater or equal 1
max_cpl_iteration = 10      # number of coupling iterations


# - - - - - - - - - - - - - - - - - - restart - - - - - - - - - - - - - - - - -
# option to restart a coupled setup from a previous run
#  -> gives the opportunity to add PISM input fluxes to MOM from last coupling 
#     iteration of the previous run
#  -> still requires to specify PISM restart files by hand (see pism_infile 
#     and pism_infile_dir below)
#  -> MOM restart files are automatically copied from $restart_dir/POEM/INPUT/
#     if path exists
coupled_restart = True
pism_to_mom_flux_restart_file = '12811.fluxes.nc'
restart_dir = "/p/tmp/kreuzer/coupled_PISM_MOM/experiments/MOM5_PISM_16km_gmd-2020-230_run04"
pism_to_mom_flux_restart_path = os.path.join(restart_dir, 'x_PISM-to-MOM', pism_to_mom_flux_restart_file)

# - - - - - - - - - - - - - - - - ocean anomaly - - - - - - - - - - - - - - -  -
do_ocean_anomaly    = False

# in case of coupled restart with do_ocean_anomaly:
#     specify ocean anomaly reference file from previous run
use_ocean_anomaly_from_prev_run = False
#ocean_anomaly_reference_file = "001234-002345.mean_MOM.nc" 
ocean_anomaly_reference_file = "12811.fluxes.nc" 
ocean_anomaly_reference_path = os.path.join(restart_dir, 'x_PISM-to-MOM', ocean_anomaly_reference_file)
#     or specify MOM output files used for computing ocean anomaly reference state
#       -> used if do_ocean_anomally==True and use_ocean_anomaly_from_prev_run==False
calc_ocn_anomaly = {}
calc_ocn_anomaly['path'] = "/p/tmp/kreuzer/coupled_PISM_MOM/experiments/MOM5_standalone_equi_10000yrs_PISM_runoff_run02_post/history"
calc_ocn_anomaly['yr_start'] = "7811"
calc_ocn_anomaly['yr_end'] = "7811"
calc_ocn_anomaly['yr_step'] = "10"
calc_ocn_anomaly['name_format_in']  = "%06g0101.ocean-yearly.nc"
calc_ocn_anomaly['name_format_out'] = "%06g-%06g.ocean-yearly.mean.nc"

# - - - - - - - - - - - - - - - basin shelf depth - - - - - - - - - - - - - - -
# option to use static basin shelf depth values read from file
#  -> gives possibility to use prescribed shelf depth values
#  -> depths are used by regriddedMOM-to-PISM_processing.py to select basin
#     temperatures and salinity from vertical profile
#  -> if option is not used then basin shelf depth is calculated each coupling
#     time step in PISM-to-MOM_processing.py

use_prescribed_basin_shelf_depth = True
prescribed_basin_shelf_depth_path = '/p/projects/climber3/kreuzer/POEM/POEM_PISM_coupling_templates/basin_shelf_depth/basin_shelf_depth_oceanic_gateways.nc'



# ------------------------------- POEM settings --------------------------------

poem_exp_dir        = os.path.join(experiment_dir, 'POEM')

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
#pism_infile_dir = "/p/tmp/albrecht/pism19/pismOut/equi/equi9000/results"
#pism_infile = "result_equi_16km_100000yrs.nc"
#pism_infile_dir = "/p/tmp/reese/pism_out/pism_025_initmip8km_ismip_merged_schmidtko_woa18_cold1.25_thkgradient_subgl_subglmelt_hmin700_decay7_7a241540/"
#pism_infile = "snapshots_112000.000.nc"
pism_infile_dir = "/p/tmp/kreuzer/coupled_PISM_MOM/experiments/MOM5_PISM_16km_gmd-2020-230_run04/PISM/results/"
pism_infile = "12811.pism_out.nc"
pism_infile_path = os.path.join(pism_infile_dir, pism_infile)

pism_atm_data_dir = os.path.join(pism_input_root_dir, "racmo_wessem")
pism_atm_file = "racmo_wessem_"+grid_id+"_mean1986_2005.nc"
#pism_atm_data_dir = os.path.join(pism_input_root_dir, "merged")
#pism_atm_file = "bedmap2_albmap_racmo_wessem_tillphi_pism_"+grid_id+".nc"
pism_atm_data_path = os.path.join(pism_atm_data_dir,pism_atm_file)

pism_ocn_data_dir = os.path.join(pism_input_root_dir, "schmidtko")
pism_ocn_file = "schmidtko_"+grid_id+"_means.nc"
#pism_ocn_data_dir = "/p/projects/pism/reese/ISMIP6_input_data/ocean_merged_schmidtko_woa18/"
#pism_ocn_file = "ocean_merged_schmidtko_woa18_initmip8km_means_amundsen-1.25.nc"
pism_ocn_data_path = os.path.join(pism_ocn_data_dir,pism_ocn_file)

pism_ocnkill_data_dir = os.path.join(pism_input_root_dir, "bedmap2")
pism_ocnkill_file = "bedmap2_"+grid_id+".nc"
#pism_ocnkill_data_dir = os.path.join(pism_input_root_dir, "merged")
#pism_ocnkill_file = "bedmap2_albmap_racmo_wessem_tillphi_pism_"+grid_id+".nc"
pism_ocnkill_data_path = os.path.join(pism_ocnkill_data_dir,pism_ocnkill_file)




# - - - - - - - - - - - - - - - - config file - - - - - - - - - - - - - - -  -

# set pism parameters that apply to all runs (unless part of the ensemble)
pism_config_file = os.path.join(pism_code_dir,"src/pism_config.cdl")

# override parameters that deviate from default.
override_params = collections.OrderedDict([


#("atmosphere.precip_exponential_factor_for_temperature", 0.04879),
#("atmosphere.precip_exponential_factor_for_temperature_doc", "ln(1.05); a 5% change of precipitation rate for every one degC of temperature change"),
#("atmosphere.precip_exponential_factor_for_temperature_type", "scalar"),
#("atmosphere.precip_exponential_factor_for_temperature_units", "Kelvin-1"),
#("atmosphere.lapse_rate.precipitation_lapse_rate", 0.),
#("atmosphere.lapse_rate.precipitation_lapse_rate_doc", "Elevation lapse rate for the surface mass balance"), 
#("atmosphere.lapse_rate.precipitation_lapse_rate_type", "scalar"),
#("atmosphere.lapse_rate.precipitation_lapse_rate_units", "(m / year) / km"),
#("atmosphere.lapse_rate.precipitation_lapse_rate_option", "precip_lapse_rate"), 
#("bed_deformation.lithosphere_flexural_rigidity", 6.e+24),
#("bed_deformation.mantle_viscosity", 1.e+21),
#("bed_deformation.update_interval", 10.),
#("surface.pdd.factor_snow", 0.0032967032967033),
#("surface.pdd.factor_ice", 0.00879120879120879),
#("surface.pdd.std_dev", 5.),
#("hydrology.tillwat_decay_rate", 1.0),
#("basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden", 0.04),
#("energy.enthalpy.temperate_ice_thermal_conductivity_ratio", 0.1), 
#("stress_balance.sia.max_diffusivity", 100.),



("atmosphere.pik.parameterization", "era_interim"),
("atmosphere.precip_exponential_factor_for_temperature", 0.04879),
("atmosphere.precip_exponential_factor_for_temperature_doc", "ln(1.05); a 5% change of precipitation rate for every one degC of temperature change"),
("atmosphere.precip_exponential_factor_for_temperature_type", "scalar"),
("atmosphere.precip_exponential_factor_for_temperature_units", "Kelvin-1"),
("atmosphere.lapse_rate.precipitation_lapse_rate", 0.),
("atmosphere.lapse_rate.precipitation_lapse_rate_doc", "Elevation lapse rate for the surface mass balance"), 
("atmosphere.lapse_rate.precipitation_lapse_rate_type", "scalar"),
("atmosphere.lapse_rate.precipitation_lapse_rate_units", "(m / year) / km"),
("atmosphere.lapse_rate.precipitation_lapse_rate_option", "precip_lapse_rate"), 

("basal_resistance.pseudo_plastic.enabled","true"),
("basal_resistance.pseudo_plastic.q", 0.75),
("basal_resistance.pseudo_plastic.u_threshold", 100.0),
("basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden", 0.04),

("bed_deformation.lc.elastic_model", "yes"),
("bed_deformation.lithosphere_flexural_rigidity", 6.e+24),
("bed_deformation.mantle_viscosity", 1.e+21),
("bed_deformation.update_interval", 10.),

("calving.methods", "eigen_calving,ocean_kill"), 
("calving.eigen_calving.K", 1.0e17), 

("energy.enthalpy.temperate_ice_thermal_conductivity_ratio", 0.1), 
("energy.basal_melt.use_grounded_cell_fraction", "false"),

("geometry.part_grid.enabled", "true"),
("geometry.remove_icebergs", "true"),
("geometry.grounded_cell_fraction", "true"),

("hydrology.model", "null"),
("hydrology.tillwat_decay_rate", 1.0),

#("ocean.models", "pico"),
("ocean.pico.heat_exchange_coefficent", 1.0e-5),
("ocean.pico.overturning_coefficent", 0.5e6),
("ocean.pico.exclude_ice_rises", "yes"),
("ocean.pico.continental_shelf_depth", -2000),

("stress_balance.calving_front_stress_bc", "true"),
("stress_balance.model","ssa+sia"),
("stress_balance.sia.enhancement_factor",1.0),
("stress_balance.ssa.enhancement_factor",1.0),
("stress_balance.ssa.method", "fd"),
("stress_balance.sia.max_diffusivity", 100.),
("stress_balance.sia.flow_law", "gpbld"),
("stress_balance.ssa.flow_law", "gpbld"),

#("surface.models", "pdd"),
("surface.pdd.factor_snow", 0.0032967032967033),
("surface.pdd.factor_ice", 0.00879120879120879),
("surface.pdd.std_dev", 5.),
#
#
##("basal_yield_stress.mohr_coulomb.topg_to_phi.enabled",  "yes"),
##("basal_yield_stress.mohr_coulomb.topg_to_phi.phi_min", 2.0), # FIXME 5.0
##("basal_yield_stress.mohr_coulomb.topg_to_phi.phi_max", 50.0),
##("basal_yield_stress.mohr_coulomb.topg_to_phi.topg_min", -700.0),
##("basal_yield_stress.mohr_coulomb.topg_to_phi.topg_max", 500.0),
### grounding line interpolations of melting
##("energy.basal_melt.use_grounded_cell_fraction", "true"),
##("energy.basal_melt.no_melting_first_floating_cell", "false"),
##("calving.methods", "eigen_calving,thickness_calving,ocean_kill"), # FIXME 
##("calving.eigen_calving.K", 1e16), # FIXME 1e17
##("calving.thickness_calving.threshold", 50), # 200
##("calving.ocean_kill.file", "initdata/"+pism_ocnkill_file), 
###("calving.float_kill.calve_near_grounding_line", "false"), # FIXME remove? Keep one shelf cell
### the following four options are equivalent to command line option -pik
### if all set to true
##("hydrology.set_tillwat_ocean", "yes"), # use Mattias tillwat fix
#### Include limit for the nomass runs! FIXME nomass only! And for Bedmachine because of convergence errors
##("stress_balance.ssa.fd.max_speed", 20e3),
##("stress_balance.sia.limit_diffusivity", "yes"),
##("stress_balance.sia.max_diffusivity", 10),
##("stress_balance.sia.surface_gradient_method", "GL_thk"),
##("stress_balance.ssa.fd.relative_convergence", "1.e-07"),
###("hydrology.use_const_bmelt", "yes"),
])

# - - - - - - - - - - - - - - - command line options - - - - - - - - - - - - - -
pism_general_opt = "-verbose 2 -options_left -o_format netcdf4_parallel -log_view"
pism_atm_opt = "-atmosphere pik -atmosphere_pik_temp_file initdata/"+pism_atm_file+" -surface pdd"
pism_add_opt = "-ocean_kill_file initdata/"+pism_ocnkill_file



# - - - - - - - - - - - - - - - - - - output - - - - - - - - - - - - - - - - -

# set extra variables
pism_extra_vars = "mask,thk,velsurf_mag,velbase_mag,flux_mag,tillwat,tauc,pico_overturning,pico_temperature_box0,pico_salinity_box0,tillphi,shelfbmassflux,shelfbtemp,basins,bmelt,bfrict,bfrict,tendency_of_ice_amount,amount_fluxes,ice_mass,pico_contshelf_mask"

# set pism diagnostic timesteps
pism_diag_extra_timestep = coupling_timestep
pism_diag_snap_timestep = coupling_timestep 
pism_diag_ts_timestep = coupling_timestep 




# ---------------------- no edits below this line needed. ----------------------
project_root = os.path.dirname(os.path.abspath(__file__))
user = pwd.getpwuid(os.getuid()).pw_name
# ------------------------------------------------------------------------------ 




########### ---- old code below ---- ###################
# select pism code version
#pism_code_version = "pism1.2.1" #"dev" # "pism1.1"

# directories
#pism_experiments_dir = os.path.join(home_dir,"pism_experiments")
#pism_code_dir = os.path.join(home_dir,"pism")
#pism_code_dir = os.path.join(home_dir_mengel,"pism")

#input_data_dir = os.path.join(pism_input_root_dir,"merged")
##input_data_dir_bedmachine = os.path.join(pism_input_root_dir, "BedMachine")

#atm_data_dir = os.path.join(pism_input_root_dir,"merged")
#ocn_data_dir = os.path.join(pism_input_root_dir,"schmidtko")
# FIXME
# 4km, 16km
#ocn_data_dir2 = os.path.join(pism_input_root_dir,"pycmip5/p003_testing")
# 8km
#ocn_data_dir2 = os.path.join(pism_input_root_dir,"pycmip5/p004_8kmprojections")
ocn_data_dir2 = "dummy_path/"

#tillwat_data_dir = os.path.join(pism_input_root_dir, "tillwat")




# select the run steps, possible are "smoothing","nomass", "full_physics", 
#                                       "forcing", "continue"
#   smoothing:      optionally run 1-year sia only to smooth the fields
#   nomass:         no_mass run for temperature evolution (could be run with 
#                   sia or sia+ssa advection, in the latter case you might 
#                   want to limit sia-diffusivity and ssa velocities using 
#                   config_override)
#   full_physics:   run with full physics, select start year and input type 
#                   below, select parametes for ensemble below
#   forcing:        run forcing experiment, select forcing files below, 
#                   FIXME this does not exist yet!
#steps = ["full_physics"]
#steps = ["forcing"]

# FIXME For full_physics, forcing: select the start year and duration
# "full_physics"
#startyear = 1000
#length = 2000
# "forcing"
#startyear = 1850
#length = 450


# Only full_physics: select the init type
# "bootstrapping": the file is bootstrapped, 
# "regrid": bootrstrapping + regridding of certain variables
# "": -i option
#init="regrid" 



# FIXME BEDMAP2 input files, used in nomass and for topography and ice thickness in full_physics (regridding) 
#bootstrapfile = os.path.join(input_data_dir,
#                      "bedmap2_albmap_racmo_wessem_tillphi_pism_"+grid_id+".nc")
# FIXME BEDMACHINE input file
#bootstrapfile = os.path.join(input_data_dir_bedmachine,
#                        "bedmachine_"+grid_id+".nc")


# "nomass" run: regrid only the tillwat variable from a fit with rignot velocities
#regridfile_tillwat = os.path.join(tillwat_data_dir,
#                       "tillwat_initmip16km_adj_rignot_tillwat100.nc")
# infile = "no_mass.nc"

## use the smoothing file created with esia=essa=1
#infile_smoothing = os.path.join(working_dir,"dev_061_initmip16km_testing_small_ensemble_dbeb47a0/",
#                      "no_mass.nc")

# infile = bootstrapfile
#infile_full_physics = os.path.join(working_dir,"picobw_052_initmip4km_testing_tillphi_tw5/no_mass_tillphi_tillwatmod.nc")
#infile_full_physics = os.path.join(store_data_dir,"pism1.1_070_initmip16km_ensemble_bedmap2_nomass_adjtillwat/no_mass_tillwat")
#infile_full_physics = os.path.join(store_data_dir,"pism1.1_070_initmip16km_ensemble_bedmap2_nomass_adjtillwat/no_mass_tillwat100")


# forcing: see below 

#atmfile = "bedmap2_albmap_racmo_wessem_tillphi_pism_"+grid_id+".nc"
#ocean_opts = "-ocean pico -ocean_pico_file $oceanfile"

# ocean_data_dir = ""
# FIXME select ocean file
# oceanfile = os.path.join(pism_input_root_dir, "schmidtko" ,"schmidtko_"+grid_id+"_means_cold.nc")
# oceanfile = os.path.join(ocn_data_dir,"schmidtko_"+grid_id+"_means_amundsen_m0p37.nc")
# 8km
#oceanfile = os.path.join(ocn_data_dir2,"thetao_Omon_GFDL-CM3_historical+rcp85_r1i1p1/schmidtko_anomaly/thetao_Omon_GFDL-CM3_historical+rcp85_r1i1p1_"+grid_id+"_100km_time0.nc") 
# 4km and 16km
#oceanfile = os.path.join(ocn_data_dir2, "thetao_Omon_GFDL-CM3_historical+rcp85_r1i1p1/schmidtko_anomaly/thetao_Omon_GFDL-CM3_historical+rcp85_r1i1p1_"+grid_id+"_100km_time0.nc")
#oceanfile = os.path.join("/gpfs/work/pn69ru/di52cok/pism_input/schmidtko/schmidtko_"+grid_id+"_means_intermediate.nc")

## forcing: ocean data iterables 4km
ocean_data_dir = ocn_data_dir2 #"/gpfs/work/pn69ru/di52cok/pism_input/pycmip5/p003_testing"
its = ["CSIRO-Mk3-6-0_historical+rcp85","GFDL-CM3_historical+rcp85","IPSL-CM5A-LR_historical+rcp85"]

iterables = {}
# FIXME include the ocean file iterables for "forcing" runs: 
#iterables["oceanfile"] = { k : os.path.join(ocean_data_dir,
#   "thetao_Omon_"+k+"_r1i1p1/schmidtko_anomaly/thetao_Omon_"+k+"_r1i1p1_"+grid_id+"_100km.nc")
#   for k in its}

## "full_physics": to create parameter ensemble
param_iterables = {}
# FIXME: include parameters for full_physics ensemble
#param_iterables["stress_balance.sia.enhancement_factor"] = [1.0,2.0]
#param_iterables["stress_balance.ssa.enhancement_factor"] = [1.0,0.4]
param_iterables["basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden"] = [0.01,0.02,0.03,0.04]
param_iterables["basal_resistance.pseudo_plastic.q"] = [0.75,0.5,0.33,0.25,0.1]
#param_iterables["hydrology.tillwat_decay_rate"] = [2,5,8]
#param_iterables["calving.eigen_calving.K"] = [1.0e16, 5.0e16, 1.0e17, 5.0e17, 1.0e18]
# special case topg_to_phi caught by if clause later:
#param_iterables["topg_to_phi"] = [
##[2.,20.,-700.,500.],
##[2.,50.,-500.,0.],
##[2.,20.,-500.,0.],
##[2.,50.,-500.,500.],
##[2.,20.,-500.,500.],
##[2.,30.,-500.,0.],
##[2.,50.,-500.,1000.]
## new:
#[1.0,50.,-700.,500.], 
#[2.0,50.,-700.,500.],
#[3.0,50.,-700.,500.],
#[4.0,50.,-700.,500.]]
# param_iterables["ocean.pico.overturning_coefficent"] = [5e5,1e6]
#param_iterables["ocean.pico.heat_exchange_coefficent"] = [1e-5,2e-5,4e-5]


## FIXME forcing: add also a control run in which the ocean data from before is used (e.g., schmidtko) 
##iterables["oceanfile"].update({"base":oceanfile})
#
## for continue_set.py and create_set_forcing.py:
## This allows to continue a number of runs as specified in runs_to_continue from the full_physics ensemble
## it is also used to create forcing runs (can also be a dummy containing only one run):
#source_ensemble_table = "pism1.1_075_initmip4km_bedmap2_testing_calvthk50_okill_tillphi2_mediumtemps"
#
## For continue_set.py:  a subset of the hashes in ensemble_table, can also be "all".
#runs_to_continue = "data/lists_of_best/pism1.1_075_initmip8km_bedmap2_mini_ensemble_thkcalv50_phimin.txt" #"data/lists_of_best/dev_058_initmip4km_resoensemble5best_20_amundsen_vel_gl.txt"
#
## For "forcing": infile is created in create_set_forcing by using get_infile_to_continue(ehash, year)
## specify the infile(s) for the forcing run  
##runs_for_forcing = "data/lists_of_best/dev_063_initmip16km_testing_small_ensemble_full_physics_forcing.txt"
#runs_for_forcing = "data/lists_of_best/pism1.1_075_initmip4km_bedmap2_testing_calvthk50_okill_tillphi2_intermediatetemps.txt"
#
#
## ensemble hash is inserted between infile_continue[0] and infile_continue[1]
## infile_continue = ["/gpfs/work/pn69ru/di36lav2/pism_store/dev_058_initmip4km_resoensemble5/dev_058_initmip4km_resoensemble5_",
## "snapshots_2300.000.nc"]
#def get_infile_to_continue(ehash, year):
#    pre = os.path.join(store_data_dir, "pism1.1_075_initmip4km_bedmap2_testing_ensemble", pism_code_version+"_071_"+grid_id+"_bedmap2_testing_calvthk50_okill_tillphi2_mediumtemps_")
#    fle = ["snapshots_",".000.nc"]
#    return os.path.join(pre+ehash,fle[0]+str(year)+fle[1])
#




