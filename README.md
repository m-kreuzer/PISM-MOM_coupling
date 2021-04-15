# PISM-MOM_coupling
code for coupling the Parallel Ice Sheet Model PISM with the Modular Ocean Model MOM

How to run a coupled MOM-PISM setup?
-------------------------------------
This is a guide how to run a coupled setup of the Modular Ocean Model MOM and the
ice-sheet model PISM on the PIK cluster HLRS2015.


## 0. Prerequisites


Most of the coupling specific scripts are written in Python and require a 
suitable environment:

Python 3 with packages 

-  numpy
-  netcdf4
-  xarray

To create such an environment with the name `py3_netcdf_xarray` locally, use 
Anaconda:

	conda create -n py3_netcdf_xarray python=3 numpy scipy netcdf4 tqdm nc-time-axis cftime=1.2.1 xarray

Anaconda should be available by default, if not type: 

	module load anaconda

(The coupling main script `run_coupled.sh` will activate the environment 
`py3_netcdf_xarray`  when it is invoked. Thus, make sure you adapt this in the 
script in case you name your environment differently or already have an 
environment fulfilling the requirements that you want to use.)


## 1. Adjust settings

Cluster specific settings like paths to MOM and PISM template directories, 
working directory for coupled runs, SLURM commands and specifications are
set in `create_experiments/pikcluster_settings.py`. 
Settings for a specific coupled run, like the coupling time step and coupling
iterations or flags for ocean anomaly approach, ocean time series forcing, ...
are set in `create_experiments/settings.py`.

## 2. create coupled run experiment from template

After adapting the setting files, a coupled experiment is created via

    python create_experiment/create_run.py

A directory with the given experiment name in `settings.py` is created in the
path specified by `working_dir` in `pikcluster_settings.py`.

## 3. run coupled setup

The created experiment includes the file `run_coupled.sh` and can be started 
via SLURM with `sbatch run_coupled.sh`

Output of the coupling routine is written in `sbatch.$SLURM_JOBID.out` which 
also contains runtime statics of the coupled setup.
Model output files are in `POEM/history/` and `PISM/results` for both models 
respectively.

