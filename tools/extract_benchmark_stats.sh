#!/bin/sh

#ROOT_WORK_DIR=$(readlink -f $PWD/..)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
ROOT_WORK_DIR=$(readlink -f $SCRIPT_DIR/..)
echo ROOT_WORK_DIR: $ROOT_WORK_DIR

SLURM_JOBID=$1
OUTFILE=$ROOT_WORK_DIR/runtime-stat.$SLURM_JOBID.out

echo SLURM_JOBID $SLURM_JOBID > $OUTFILE
echo $PWD >> $OUTFILE

tail -n 26 $ROOT_WORK_DIR/sbatch.$SLURM_JOBID.out | head -n 20 >> $OUTFILE




