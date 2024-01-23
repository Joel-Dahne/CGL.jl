#!/bin/bash -l
# See https://www.kth.se/blogs/pdc/2018/08/getting-started-with-slurm/
#SBATCH --job-name branch_points
#SBATCH --account naiss2024-22-13
#SBATCH --mail-type=ALL
#SBATCH --time 01:00:00

#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=8

#SBATCH -o Dardel/logs/branch_points.o
#SBATCH -e Dardel/logs/branch_points.e

time julia --project=. Dardel/scripts/branch_points.jl "$@"
