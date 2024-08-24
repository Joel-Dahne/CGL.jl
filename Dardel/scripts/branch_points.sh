#!/bin/bash -l
# See https://www.kth.se/blogs/pdc/2018/08/getting-started-with-slurm/
#SBATCH --job-name branch_points
#SBATCH --account naiss2024-22-1038
#SBATCH --mail-type=ALL
#SBATCH --time 01:00:00

#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=16

#SBATCH -o Dardel/logs/branch_points.o
#SBATCH -e Dardel/logs/branch_points.e

if [ -z "${CGL_SLURM_MEM_PER_NODE}" ]; then
    # This is the amount of memory to use per node in GB. It needs to
    # be tuned to the cluster.
    export CGL_SLURM_MEM_PER_NODE=220
fi
time julia --project=. Dardel/scripts/branch_points.jl "$@"
