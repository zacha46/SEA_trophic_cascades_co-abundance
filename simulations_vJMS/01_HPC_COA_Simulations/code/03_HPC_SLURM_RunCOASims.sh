#!/bin/bash
# SPECIFY YOUR GROUP NAME
#SBATCH --account=a_luskin_ecl

## SPECIFY YOUR COMPUTATIONAL REQUIREMENTS 

# Select 1 node per job
#SBATCH --nodes=1

# Select 1 task per node (b/c R is not MPI)
#SBATCH --ntasks=1

# Select 4 CPUS per node (one per MCMC chain)
#SBATCH --cpus-per-task=4

# Select 50000 MB (50 GB) of memory per node 
#SBATCH --mem=50000

# Ensure we are in the general queue, not AI, debug, or GPU
#SBATCH --partition=general

# Select 4 hours (h:m:s format) of walltime for small mods; we should not need that long
#SBATCH --time=20:00:00

# SPECIFY THE JOB ARRAY-
#SBATCH --array=1-924

# SPECIFY THE JOB NAME
#SBATCH --job-name=COASims

# SPECIFY .err AND .out FILE LOCATIONS
#SBATCH --output=/scratch/user/s4633921/01_HPC_COA_Simulations/results/OUT/slurm-%A_%a.out
#SBATCH --error=/scratch/user/s4633921/01_HPC_COA_Simulations/results/ERRORS/slurm-%A_%a.err

# SET THE 'setting' VARIABLE THAT WILL BE LOADED IN R
export SETTING="SIMULATION" 

# LOAD THE RELEVANT MODULE
module load rjags/4-10-foss-2021a-r-4.1.0 
if [ $? -ne 0 ]; then 
  echo "Module load failed"
  exit 3
fi

# SPECIFY THE PBS WORKING DIRECTORY AND PRINT TO VERIFY
cd $SLURM_SUBMIT_DIR
pwd

# Compute iteration range for this array task
start=$(( ($SLURM_ARRAY_TASK_ID - 1) * 25 + 1 ))
end=$(( $SLURM_ARRAY_TASK_ID * 25 ))

echo "Running iterations $start to $end in job ID $SLURM_JOB_ID"

# LOAD THE R SCRIPT, SUBMIT JOB AND SPECIFY THE ARRAY INDEX
for i in $(seq $start $end); do
   Rscript /scratch/user/s4633921/01_HPC_COA_Simulations/code/02_Simulations_FitModelsHPC.R $i 
done




