#!/bin/bash -l
# See https://www.kth.se/blogs/pdc/2018/08/getting-started-with-slurm/
#SBATCH --job-name branch
#SBATCH --account naiss2024-22-13
#SBATCH --mail-type=ALL
#SBATCH --time 01:00:00

#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=16

#SBATCH -o Dardel/logs/branch.o
#SBATCH -e Dardel/logs/branch.e

time julia --project=. Dardel/scripts/branch_existence.jl "$@"
