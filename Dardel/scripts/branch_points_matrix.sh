#!/bin/bash -l

# This needs to be updated to the account of your cluster
#SBATCH --account naiss2024-22-1038

# This needs to be updated based on your cluster
#SBATCH --partition main
#SBATCH --nodes 1
#SBATCH --ntasks 16
#SBATCH --cpus-per-task 16

# This likely does not need to be updated
#SBATCH --mail-type ALL
#SBATCH --job-name branch_points_matrix
#SBATCH --output Dardel/logs/%x.o
#SBATCH --error Dardel/logs/%x.e

if [ -z "${CGL_SLURM_MEM_PER_NODE}" ]; then
    # This is the amount of memory to use per node in GB. It needs to
    # be tuned to the cluster.
    export CGL_SLURM_MEM_PER_NODE=220
fi

time julia --project=. Dardel/scripts/branch_points_matrix.jl "$@"
