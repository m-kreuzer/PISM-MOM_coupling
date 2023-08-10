### machine-dependent settings ###

# ----------------------------------- paths ------------------------------------
# base directory where experiments will be located (submit/input/output files)
#working_dir = "/p/projects/pism/kreuzer/coupled_PISM_MOM/experiments"
working_dir = "/p/tmp/kreuzer/coupled_PISM_MOM/experiments"

# directory with coupling infrastructure template 
#  -> clone of https://github.com/m-kreuzer/PISM-MOM_coupling.git
coupl_template_dir = "/p/projects/pism/kreuzer/software/PISM-MOM_coupling"

### PISM paths
#pism_code_dir = "/p/projects/pism/kreuzer/software/pism1.2.1"
#pism_code_dir = "/p/projects/climber3/kreuzer/PISM/pism1.1"
#pism_code_dir = "/p/projects/pism/kreuzer/software/pism1.1.4_extra-double"
pism_code_dir = "/home/kreuzer/projects/phd/pism_1.0_precip_scale_hash/"
#pism_code_dir = "/home/reese/pism_code/pism2019/pism"
#pism_input_root_dir = "/p/projects/pism/kreuzer/pism_input_2020"
pism_input_root_dir = "/p/tmp/kreuzer/coupled_PISM_MOM/experiments/pism1.0_precipscale_hash_equi_16km_run10/initdata"


### POEM paths
# POEM project directory with tools like mppnccombine, time_stamp.csh, ...
#poem_proj_dir = "/p/projects/climber3/kreuzer/POEM/trunk"
poem_proj_dir = "/p/projects/climber3/kreuzer/POEM/branches/kreuzer_mom5.1"
# POEM template directory with all namelist & input files
poem_template_dir = "/p/projects/climber3/kreuzer/POEM/POEM_PISM_coupling_templates/CM2M_coarse-EM3-no_antarctic_river_runoff-basal"



# ---------------------------------- commands ---------------------------------- 
pism_exec = "pismr"
pism_mpi_do = "srun -n"
submit_command="sbatch run_coupled.sh"


# ---------------------------- resource management ----------------------------- 
slurm_nodes = 4
#slurm_tasks = 32
#slurm_partition="broadwell"
slurm_exclusive = True

#slurm_qos = 'short'
#slurm_time = "23:00:00"
slurm_qos = 'medium'
slurm_time = "4-23:00:00"
# Acceptable time formats include "minutes", "minutes:seconds", "hours:minutes:seconds",
#                   "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"

#slurm_add_directives = "#SBATCH --constraint=haswell,tasksmax"
#slurm_add_directives = "#SBATCH --constraint=tasksmax"
slurm_add_directives = ""
