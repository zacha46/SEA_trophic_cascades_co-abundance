#!/bin/bash -l

# SPECIFY YOUR GROUP NAME
#SBATCH --account=a_luskin_ecl

## SPECIFY YOUR COMPUTATIONAL REQUIREMENTS 

# Select 1 node per job
#SBATCH --nodes=1

# Select 1 task per CPU (b/c R is not MPI)
#SBATCH --ntasks=1

# Select 4 CPUS (aka threads) per node (one per MCMC chain)
#SBATCH --cpus-per-task=4

# Select 50,000 MB (50 GB) of memory per node 
#SBATCH --mem=50000

# Ensure we are in the general queue, not AI, debug, or GPU
#SBATCH --partition=general

# Select 8 hours (h:m:s format) of walltime for small mods
#SBATCH --time=8:00:00

# SPECIFY THE JOB ARRAY-
#SBATCH --array=1-64

# SPECIFY THE JOB NAME
#SBATCH --job-name=CF3

# SPECIFY .err AND .out FILE LOCATIONS
#SBATCH --output=OE/output/slurm-%A_%a.out
#SBATCH --error=OE/error/slurm-%A_%a.err

# LOAD THE RELEVANT MODULE
module load rjags/4-10-foss-2021a-r-4.1.0 
if [ $? -ne 0 ]; then 
  echo "Module load failed"
  exit 3
fi

# SET THE 'setting', 'pref' & 'gb' VARIABLES THAT WILL BE LOADED IN R
export SETTING="SHORT" 

# SET THE 'counter' VARIABLE THAT WILL BE LOADED IN R
export COUNTER="counterfactual3_isolate_altitude"

# SPECIFY THE PBS WORKING DIRECTORY AND PRINT TO VERIFY
cd $SLURM_SUBMIT_DIR
pwd

# LOAD THE R SCRIPT, SUBMIT JOB FROM /scratch/user/uqzamir/, AND SPECIFY THE ARRAY INDEX 
srun Rscript code/HPC_co-abundance_model_counterfactuals.R $SLURM_ARRAY_TASK_ID