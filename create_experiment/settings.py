###
# Set basic settings here,e.g., adjust path to input files, 
# select the run type, ensemble parameters...
# runscripts are created via create_run.py or create_set.py
# more options are set in templates/pism_run.sh.jinja2

import os
import pwd
import collections
import pwd

import pism_grids

# import settings including path to input, output directories
from pikcluster_settings import *


# ----------------------------- coupling settings ------------------------------

# ATTENTION: make sure to adjust this, otherwise, files will be overwritten
# a useful approach is to have one number (_061_) for a suite of runs that get a
# common name (_small_ensemble_) and an additional identifies for the current 
# run step (_forcing_) 
# FIXME
#experiment = pism_code_version+"_075_"+grid_id+"_bedmachine_ensemble_amedtem" # no _
experiment = "coupled_run_setup_test"
#experiment_dir      = os.path.join(settings.working_dir, settings.experiment)
#experiment_dir      = os.path.join(working_dir, experiment)


coupling_timestep = 10      # in years, must be greater or equal 1
max_cpl_iteration = 40      # number of coupling iterations


# ------------------------------- POEM settings --------------------------------



# ------------------------------- PISM settings --------------------------------

# select resolution of the run
grid_id = "initmip16km" 
pism_grid = pism_grids.grids[grid_id]



# input data
pism_infile_dir = "/p/tmp/albrecht/pism19/pismOut/equi/equi9000/results"
pism_infile = "result_equi_16km_100000yrs.nc"
pism_infile_path = os.path.join(pism_infile_dir, pism_infile)

pism_atm_data_dir = os.path.join(pism_input_root_dir, "racmo_wessem")
pism_atm_file = "racmo_wessem_"+grid_id+"_mean1986_2005.nc"
pism_atm_data_path = os.path.join(pism_atm_data_dir,pism_atm_file)

pism_ocn_data_dir = os.path.join(pism_input_root_dir, "schmidtko")
pism_ocn_file = "schmidtko_"+grid_id+"_means.nc"
pism_ocn_data_path = os.path.join(pism_ocn_data_dir,pism_ocn_file)

pism_ocnkill_data_dir = os.path.join(pism_input_root_dir, "bedmap2")
pism_ocnkill_file = "bedmap2_"+grid_id+".nc"
pism_ocnkill_data_path = os.path.join(pism_ocnkill_data_dir,pism_ocnkill_file)


# set pism parameters that apply to all runs (unless part of the ensemble)
pism_config_file = os.path.join(pism_code_dir,"src/pism_config.cdl")

# override parameters that deviate from default.
override_params = collections.OrderedDict([
# "ocean.pico.continental_shelf_depth", -2000,
("stress_balance.sia.enhancement_factor",1.0),
("stress_balance.ssa.enhancement_factor",1.0),
("stress_balance.model","ssa+sia"),
("time_stepping.skip.enabled", "yes"),
("basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden", 0.03),
("basal_resistance.pseudo_plastic.q", 0.75),
("basal_yield_stress.mohr_coulomb.topg_to_phi.enabled",  "yes"),
("basal_yield_stress.mohr_coulomb.topg_to_phi.phi_min", 2.0), # FIXME 5.0
("basal_yield_stress.mohr_coulomb.topg_to_phi.phi_max", 50.0),
("basal_yield_stress.mohr_coulomb.topg_to_phi.topg_min", -700.0),
("basal_yield_stress.mohr_coulomb.topg_to_phi.topg_max", 500.0),
("basal_resistance.pseudo_plastic.enabled","true"),
("hydrology.tillwat_decay_rate", 5.0),
# grounding line interpolations of melting
("energy.basal_melt.use_grounded_cell_fraction", "false"),
("calving.methods", "ocean_kill,thickness_calving,eigen_calving"), # FIXME 
("calving.eigen_calving.K", 1e16), # FIXME 1e17
("calving.thickness_calving.threshold", 50), # 200
#("calving.ocean_kill.file", os.path.join(input_data_dir,
#                      "bedmap2_albmap_racmo_wessem_tillphi_pism_"+grid_id+".nc")),
#("calving.float_kill.calve_near_grounding_line", "false"), # FIXME remove? Keep one shelf cell
# the following four options are equivalent to command line option -pik
# if all set to true
("stress_balance.calving_front_stress_bc", "true"),
("geometry.part_grid.enabled", "true"),
("geometry.remove_icebergs", "true"),
("geometry.grounded_cell_fraction", "true"),
("ocean.pico.exclude_ice_rises", "yes"),
#("hydrology.set_tillwat_ocean", "yes"), # use Mattias tillwat fix
## Include limit for the nomass runs! FIXME nomass only! And for Bedmachine because of convergence errors
("stress_balance.ssa.fd.max_speed", 10e3),
("stress_balance.sia.limit_diffusivity", "yes"),
("stress_balance.sia.max_diffusivity", 10),
#("hydrology.use_const_bmelt", "yes"),
])





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




