#!/bin/bash

{% if settings.is_pikcluster: -%}
#SBATCH --qos=short
#SBATCH --time=0-23:50:00
#SBATCH --job-name={{settings.experiment}}
#SBATCH --account=ice
#SBATCH --output=./slurm_out.out
#SBATCH --error=./slurm_error.err
#SBATCH --ntasks=160
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=reese@pik-potsdam.de
{% if settings.pik_partition=="broadwell" -%}
#SBATCH --tasks-per-node=32
#SBATCH --partition=broadwell
{% else %}
#SBATCH --tasks-per-node=16
{% endif %}

module purge
module load pism/stable1.0
runname=`echo $PWD | awk -F/ '{print $NF}'`
outdir={{settings.working_dir}}/$runname

mkdir -p $outdir/log/

# make the PISM execution script aware that it is on compute nodes.
export PISM_ON_CLUSTER=1
./pism_run.sh $SLURM_NTASKS > $outdir/log/pism.out

{% else %}

#@ job_name = smuc_$(cluster)_$(stepid)
#@ class = micro
#@ group = pn69ru
#@ notify_user = reese@pik-potsdam.de
#@ job_type = MPICH
#@ output = ./loadl.out
#@ error  = ./loadl.err
#@ wall_clock_limit = 47:59:00
#@ notification=always
#@ network.MPI = sn_all,not_shared,us
#@ node = 7
#@ tasks_per_node = 28
#@ island_count = 1
#@ energy_policy_tag = albrecht_pism_2015
#@ minimize_time_to_solution = yes
#@ queue

source /etc/profile.d/modules.sh
module unload mpi.ibm intel mkl netcdf
module load mkl/2018 intel/18.0 mpi.intel/2018
module load hdf5/mpi/1.10.1
module load gcc/8
module load gsl/2.4
module load fftw/mpi/v3_mkl2018

#NETCDF_ROOT=/home/hpc/pn69ru/di36lav2/software/netcdf-4.6.1-intel-2018
NETCDF_ROOT=/home/hpc/pn69ru/di36lav2/software/netcdf-c-4.6.2-intel-2018
export PATH=$NETCDF_ROOT/bin:$PATH

export PETSC_DIR=/home/hpc/pn69ru/di36lav2/petsc/petsc-3.9.2-intel2018
export PETSC_ARCH=arch-linux2-c-opt

# get user and computer specific variables like paths
runname=`echo $PWD | awk -F/ '{print $NF}'`
outdir={{settings.working_dir}}/$runname

echo $LOADL_STEP_ID
echo $HOME
echo $LOADL_PROCESSOR_LIST
echo $LOADL_TOTAL_TASKS
echo $outdir

number_of_cores=`echo $LOADL_PROCESSOR_LIST | wc -w`
echo $number_of_cores

mkdir -p $outdir/log/

export PISM_ON_CLUSTER=1
#run_smoothing_nomass.sh $LOADL_TOTAL_TASKS > $outdir/log/smoothing_nomass.out
pism_run.sh $LOADL_TOTAL_TASKS >> $outdir/log/pism.out
{% endif %}

