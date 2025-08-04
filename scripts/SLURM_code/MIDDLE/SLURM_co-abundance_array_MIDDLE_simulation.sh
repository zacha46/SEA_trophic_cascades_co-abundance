#!/bin/bash
# SPECIFY YOUR GROUP NAME
#SBATCH --account=a_luskin_ecl

## SPECIFY YOUR COMPUTATIONAL REQUIREMENTS 

# Select 1 node per job
#SBATCH --nodes=1

# Select 1 task per node (b/c R is not MPI)
#SBATCH --ntasks=1

# Select 4 CPUS per node (one per MCMC chain)
#SBATCH --cpus-per-task=3

# Select 50,000 MB (50 GB) of memory per node 
#SBATCH --mem=50000

# Ensure we are in the general queue, not AI, debug, or GPU
#SBATCH --partition=general

# Select 2 hours (h:m:s format) of walltime for mid-setting models 
#SBATCH --time=2:00:00

# SPECIFY THE JOB ARRAY-
#SBATCH --array=1-30

# SPECIFY THE JOB NAME
#SBATCH --job-name=CoA-SIM

# SPECIFY .err AND .out FILE LOCATIONS
#SBATCH --output=OE/output/slurm-%A_%a.out
#SBATCH --error=OE/error/slurm-%A_%a.err

# LOAD THE RELEVANT MODULE
module load rjags/4-10-foss-2021a-r-4.1.0 
if [ $? -ne 0 ]; then 
  echo "Module load failed"
  exit 3
fi

# SET THE 'setting' VARIABLE THAT WILL BE LOADED IN R
export SETTING="MIDDLE" 

# SPECIFY THE PBS WORKING DIRECTORY AND PRINT TO VERIFY
cd $SLURM_SUBMIT_DIR
pwd

# LOAD THE R SCRIPT, SUBMIT JOB FROM /scratch/user/uqzamir/, AND SPECIFY THE ARRAY INDEX 
Rscript code/HPC_co-abundance_model_simulations.R $SLURM_ARRAY_TASK_ID
