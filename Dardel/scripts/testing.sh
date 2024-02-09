#!/bin/bash -l
# See https://www.kth.se/blogs/pdc/2018/08/getting-started-with-slurm/
#SBATCH --job-name run_proof
#SBATCH --account naiss2024-22-13
#SBATCH --mail-type=ALL
#SBATCH --time 00:10:00

#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=16

#SBATCH -o Dardel/logs/testing.o
#SBATCH -e Dardel/logs/testing.e

if [ -z "${CGL_SLURM_MEM_PER_NODE}" ]; then
    # This is the amount of memory to use per node in GB. It needs to
    # be tuned to the cluster.
    export CGL_SLURM_MEM_PER_NODE=220
fi
time julia --project=. Dardel/scripts/testing.jl "$@"
