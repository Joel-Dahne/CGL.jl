#!/bin/bash -l
# See https://www.kth.se/blogs/pdc/2018/08/getting-started-with-slurm/
#SBATCH --job-name run_proof
#SBATCH --account naiss2024-22-13
#SBATCH --mail-type=ALL
#SBATCH --time 04:00:00

#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=4

#SBATCH -o Dardel/logs/testing.o
#SBATCH -e Dardel/logs/testing.e

if [ -z "${CGL_SLURM_MEM_PER_NODE}" ]; then
    # This is the amount of memory per node in GB. It needs to be
    # tuned to the cluster.
    export CGL_SLURM_MEM_PER_NODE=256
fi
time julia --project=. Dardel/scripts/testing.jl "$@"
