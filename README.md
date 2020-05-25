# PISM-MOM_coupling
code for coupling the Parallel Ice Sheet Model PISM with the Modular Ocean Model MOM

-------------------------------------
How to run a coupled POEM-PISM setup?
-------------------------------------
      
This is a guide how to run a coupled setup of PIK's climate model POEM and the
ice-sheet model PISM on the PIK cluster HLRS2015.


## 0. Prerequisites


Most of the coupling specific scripts are written in Python and require a 
suitable environment:

Python 3 with packages 

-  numpy  
-  netcdf4

To create such an environment with the name `py3_netcdf` locally, use 
Anaconda:

	conda create -n py3_netcdf python=3 numpy netcdf4

Anaconda should be available by default, if not type: 

	module load anaconda

(The coupling main script `run_coupled.sh` will activate the environment 
`py3_netcdf` through `source activate py3_netcdf` when it is invoked. Thus, 
make sure you adapt this in the script in case you name your environment 
differently or already have an environment fulfilling the requirements that 
you want to use.)


## 1. Setting up POEM

The coupling setup has been developed for coupling the climate model POEM 
(MOM5 ocean incl. SIS sea ice + Aeolus atmosphere + Atlantes land components) 
with the ice-sheet model PISM via it's submodel, the sub-shelf cavity model 
PICO. Using any other setup combination based on GFDL's MOM5 ocean model 
including the FMS coupler (minor changes have been made to FMS) instead of POEM
should also work (e.g. different atmosphere model).

A detailed description of how to run POEM has been written by Stefan Petri 
(PIK). It can by found in PIK's subversion system:

    http://subversion.pik-potsdam.de/xsvn/climber3/docs/HOWTO_run_POEM.txt

or in PIK's redmine system (in the POEM project):

    https://redmine.pik-potsdam.de/dmsf/files/1559/view

As the coupling infrastructure of POEM and PISM uses the directory structure of 
POEM, this has to be set up first. So follow the instructions in 
`HOWTO_run_POEM.txt` and 

- checkout the code from PIK's subversion server to your preferred location
    (the location will be refered to as `POEM_root/` in the following)
- compile a model configuration (MOM-coupled-with-Aeolus has been tested) 
    (minor changes have been made to `src/coupler/flux_exchange.F90` to allow 
     additional fluxes from the ice to the ocean model)
- setup an experiment in `work/` subdirectory
- make a test run of your setup (and be sure the model does not crash)
    

According to the How-to, you have a repository now which looks like this:

	 POEM_root
	 +-- bin/
	 +-- data/ 
	 +-- exec/
	 +-- exp/
	 |   +-- mom4_Aeolus_lad-init-21-mar-clim/
	 |
	 +-- src/
	 +-- work/
	     +-- run-001/
	         +-- data_table
	         +-- diag_table
	         +-- field_table
	         +-- fms_MOM_LAD_AEOLUS.x-par
	         +-- history/
	         +-- INPUT/
	         +-- RESTART
	         +-- submit-par.slurm


To continue with the coupled setup, the ice-sheet model PISM is required.


## 2. Setting up PISM

The ice-sheet model PISM can be downloaded from github.
Create a new directory for that, possibly next to `POEM_root/`    
Accordingly this location will be called PISM_root in the following.

	mkdir PISM_root

and clone the github directory
	
	git clone https://github.com/pism/pism.git PISM_root

See also the documentation for more information:

    http://pism-docs.org/wiki/doku.php?id=stable_version

At PIK cluster you can use a predefined cluster environment module which loads 
all required prerequisites for compiling like compilers, netCDF library, PETSC, 
etc.

	module load pism/stable1.2

For compiling copy the following compile script to your folder:

	cp /p/projects/pism/kreuzer/software/pism1.2.1/BuildPismOnCluster.sh
	cd PISM_root/
	./BuildPismOnCluster.sh > BuildPismOnCluster.log 2>&1 &


## 3. Coupled Setup



### 3.1 POEM in Coupled Setup


In your `POEM_root` directory create a new directory in `work/` for the coupled 
setup and copy the template from `exp/`:

    mkdir -p POEM_root/work/POEM_PISM_coupled/run_001
    cd POEM_root/work/POEM_PISM_coupled/run_001

    cp -rp ../../../exp/POEM_PISM_coupled/* .

Now you have the directory structure for the coupled setup including the 
coupling scripts in your newly created folder, where you can run experiments. 
Your directory should now look like this:
 
	 POEM_root/work/POEM_PISM_coupled/run_001
	 +-- inter-model-processing/
	 |   +-- PISM-to-MOM_processing.py
	 |   +-- regriddedMOM-to-PISM_processing.py
	 |
	 +-- PISM/
	 |   + initdata/
	 |   + prerun/
	 |   + results/
	 |
	 +-- POEM/
	 |   + data_table -> data_table-dummy
	 |   + data_table-dummy
	 |   + diag_table-prerun
	 |   + input.nml -> input.nml-regular
	 |   + input.nml-regular
	 |
	 +-- pre-processing/
	 |   +-- find_edge.py
	 |   +-- PISMbasin-to-MOMcell_mapping.py
	 |
	 +-- run_coupled.sh
	 +-- README-POEM_PISM_coupled.txt
	 +-- x_MOM-to-PISM/
	 +-- x_PISM-to-MOM/
	
In directories POEM/ and PISM/ input and output files are kept, for POEM also 
configuration files. `POEM/` corresponds to a regular folder for running 
experiments like
	
    POEM_root/work/mom4_Aeolus_lad-init-21-mar-clim/run-001

if you followed 3.1.1 in `HOWTO_run_POEM.txt` 

The bash script `run_coupled.sh` is the main coupling file and *.py files are 
used for pre- or inter-modelling processing. 

Like in the previous uncoupled POEM case (see 1.), copy your compiled POEM 
executable to the POEM/ folder in POEM_PISM_coupled/run_001:

    cd POEM
    cp -p ../../../../exec/pik-hlrs2015-ifort/MOM_LAD_AEOLUS/fms_MOM_LAD_AEOLUS.x .

 as well as the configuration files from the master experiment template
 
    cp -av ../../../../exp/mom4_Aeolus_lad-init-21-mar-clim/diag_table* .
    cp -av ../../../../exp/mom4_Aeolus_lad-init-21-mar-clim/field_table .

 Create a symbolic link to the diag_table which specifies the variables and 
 frequency of POEM diagnostic file output:
 
    ln -s diag_table-monthly diag_table

Make sure all configuration files are set:

	data_table -> data_table-dummy
	diag_table -> diag_table-monthly
	field_table
	input.nml -> input.nml-regular
	
In the coupled setup special files are needed for `input.nml` and `data_table`. 
These are provided in `exp/POEM_PISM_coupled/POEM/` and are referenced via 
symbolic links.
    
*input.nml-regular*
> The coupling main script will copy this file and modify it by 
> overwriting the model integration time. So you can also use your own 
> input namelist or modify it of course, but it is important that there is 
> the symbolic link
> 
>     input.nml -> input.nml-regular
 
*data_table-dummy*
> This data table specifies paths with `FILE_NAME_DUMMY` inside, which 
> are replaced by the main coupling script `run_coupled.sh`. Also a 
> modified version is used for POEM preruns. Likewise to 
> input.nml-regular it is important that there points a symbolic link 
> named
>
>     data_table -> data_table-dummy
>
> (Make sure this was not overwritten when copying files from the master 
> experiment template.)
        
*diag_table-prerun*
> This flavor of the diag_table is part of the POEM_PISM_coupled template 
> (see directory structure above). It specifies a minimal set of 
> variables to be written out needed for pre-processing. The main 
> coupling script run_coupled.sh will put a symbolic link of diag_table to 
> it when needed.



Like in `HOWTO_run_POEM.txt`, create sub-directories for input and output files in 
`POEM/`:
 
    mkdir INPUT RESTART history

 and make symbolic links for the input files:
 
    cd INPUT
    ln -s ../../../../../exp/mom4_Aeolus_lad-init-21-mar-clim/INPUT/* .
    cd ../..

The POEM model uses now the restart conditions linked in the experiment master 
folder.


### 3.2 PISM in Coupled Setup

Also PISM needs input files to run. An example set can be downloaded from Torsten 
Albrecht's webserver:

    cd PISM/initdata
    wget -r -np -A *.nc,*.cdl -nH --cut-dirs=3 http://www.pik-potsdam.de/~albrecht/pism_equilibrium/init_data/ 
    cd ../..


### 3.3 Coupling Settings

Now a few adjustments need to be done in the main coupling script. Open 
`run_coupled.sh` in your favorite editor and 

 - set the project paths to your parent folders of `POEM_root/` and `PISM_root/` (see 1. and 2.) 
	 
	    POEM_PROJ_DIR=/path/to/POEM
	    PISM_PROJ_DIR=/path/to/PISM

    (Use absolute path names.)


 - make sure your POEM executable (or a link to it) in  

        POEM_root/work/POEM_PISM_coupled/run_001/POEM 

    is named 

        fms_MOM_LAD_AEOLUS.x-par

    as this is called in `run_coupled.sh` to execute POEM.



 - set PISM grid parameters and input file links

        # PISM geometric parameters
        P_Mx=381
        P_My=381
        P_Lz=6000
        P_Lbz=2000
        P_Mz=81
        P_Mbz=21

    You can also change the initial restart file for PISM and initial ocean 
    input file for PICO through variables:
  
        PISM_RESTART_FILE
        PISM_OCEAN_START_FILE
  
    By default they are pointing to the files in the `PISM/initdata/` directory 
    which have been downloaded in Section 2.

    Further PISM options are set directly as command line arguments after `srun` 
    in `pism_run()` and `pism_prerun()`.
    After PISM argument `-extra_vars`, variables are mentioned to write to the 
    file specified with `-extra_file` in frequency of `-extra_times`. This is used 
    to ouput the `basal_mass_flux_floating` which is used in the inter-model 
    processing scripts to calculate energy and mass flux between POEM and PISM. 

 - set parameters for the coupling time step (in years) and the maximum number 
   of coupling iterations
   
        CPL_TIMESTEP=10
        MAX_CPL_ITERATION=10

 - set number of cores to use. Make sure the arbitrary function allocates the 
   same amount of CPUs on as many nodes as required (without exceeding the 
   maximum number of CPUs per node)
  
        #SBATCH --tasks 48
        local arb="`arbitrary 16 16 16"


As in the POEM standalone run, parameters of the ocean, atmosphere, sea ice or 
land model, as well as for the POEM internal coupler FMS, can be made in POEM's 
namelist in

        POEM_root/work/POEM_PISM_coupled/run_001/POEM/input.nml-regular

The number of cores specified there must not exceed the maximum requested CPUs 
with `#SBATCH --tasks`


### 3.4 Run the coupled setup


To run the coupled setup on the cluster, the `run_coupled.sh` script must be 
submitted via SLURM:

    sbatch run_coupled.sh

It creates an output file of coupled_main.sh in 

    POEM_root/work/POEM_PISM_coupled/run_001
    
named `sbatch.SLURMJOBID.out` which also contains runtime statics of the coupled
setup.

The output files are in `POEM/history/` and `PISM/results` for both models 
respectively.


## 4. Postprocessing

When analysing the results of both models, it might be useful to synchronise 
output files to the same time axis. For this the command following NCO command
can be used:

    ncap2 -s time-=100.0 -s time_bounds-=100.0 infile.nc -O -o outfile.nc

which shifts the variables time and time_bounds of `infile.nc` by 100 units. 
For big shifts (in unit seconds) put a floating point even when shifting an 
integer (prevents integer overflow). The NCO module is required, so do
`module load nco` if not yet loaded.


