#!/bin/bash -l

# This needs to be updated to the account of your cluster
#SBATCH --account sverakv

# This needs to be updated based on your cluster
#SBATCH --partition agsmall
#SBATCH --nodes 1
#SBATCH --ntasks 16
#SBATCH --cpus-per-task 8

# This likely does not need to be updated
#SBATCH --mail-type ALL
#SBATCH --job-name branch_points_verification
#SBATCH --output HPC/logs/%x.o
#SBATCH --error HPC/logs/%x.e

if [ -z "${CGL_SLURM_MEM_PER_NODE}" ]; then
    # This is the amount of memory to use per node in GB. It needs to
    # be tuned to the cluster.
    export CGL_SLURM_MEM_PER_NODE=440
fi

time julia --project=. HPC/scripts/branch_points_verification.jl "$@"
