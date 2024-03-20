# Computations on Dardel
Some of the computations for the computer assisted part of the proof
were done on the
[Dardel<](https://www.pdc.kth.se/hpc-services/computing-systems/about-the-dardel-hpc-system-1.1053338)
computer cluster, This directory contains the scripts used for those
computations.

The computations were enabled by resources provided by the National
Academic Infrastructure for Supercomputing in Sweden (NAISS) at PDC
partially funded by the Swedish Research Council through grant
agreement no. 2022-06725.

The Dardel cluster makes use of Slurm for job scheduling. The scripts
for running the computations are written in terms of Slurm scripts.
From a fresh clone of this repository an initial setup is done with
the following commands, executed from the root of this repository.

``` shell
# Install all required packages for Julia - this takes some time
julia --project=. --eval 'using Pkg; Pkg.instantiate()'
```

The computations are run using several scripts found in
`Dardel/scripts/`.

## Verification of points

``` shell
# The case d = 1
# Fix epsilon (the default)
sbatch --time=0:15:00 -p main --job-name=branch_points_1 -o Dardel/logs/branch_points_1.o -e Dardel/logs/branch_points_1.e Dardel/scripts/branch_points.sh 1

# Fix kappa
sbatch --time=0:15:00 -p main --job-name=branch_points_1_fix_kappa -o Dardel/logs/branch_points_1_fix_kappa.o -e Dardel/logs/branch_points_1_fix_kappa.e Dardel/scripts/branch_points.sh 1 1

# The case d = 3
# Fix epsilon (the default)
sbatch --time=0:30:00 -p main --job-name=branch_points_3 -o Dardel/logs/branch_points_3.o -e Dardel/logs/branch_points_3.e Dardel/scripts/branch_points.sh 3
# Fix kappa
sbatch --time=1:00:00 -p main --job-name=branch_points_3_fix_kappa -o Dardel/logs/branch_points_3_fix_kappa.o -e Dardel/logs/branch_points_3_fix_kappa.e Dardel/scripts/branch_points.sh 3 1
```

## The branch d = 1, j = 3

### Existence of branch

``` shell
# The case d = 1, j = 3
sbatch --time=0:15:00 -p main --job-name=branch_existence_3_top -o Dardel/logs/branch_existence_3_top.o -e Dardel/logs/branch_existence_3_top.e Dardel/scripts/branch_existence.sh 3 1 top

sbatch --time=2:00:00 -p main --job-name=branch_existence_3_turn -o Dardel/logs/branch_existence_3_turn.o -e Dardel/logs/branch_existence_3_turn.e Dardel/scripts/branch_existence.sh 3 1 turn

sbatch --time=4:00:00 -p main --job-name=branch_existence_3_bottom -o Dardel/logs/branch_existence_3_bottom.o -e Dardel/logs/branch_existence_3_bottom.e Dardel/scripts/branch_existence.sh 3 1 bottom
```

### Continuation of branch

``` shell
# The case d = 1, j = 3
sbatch --time=0:15:00 -p main --job-name=branch_continuation_3_top -o Dardel/logs/branch_continuation_3_top.o -e Dardel/logs/branch_continuation_3_top.e Dardel/scripts/branch_continuation.sh 3 1 top

sbatch --time=1:00:00 -p main --job-name=branch_continuation_3_turn -o Dardel/logs/branch_continuation_3_turn.o -e Dardel/logs/branch_continuation_3_turn.e Dardel/scripts/branch_continuation.sh 3 1 turn

sbatch --time=2:00:00 -p main --job-name=branch_continuation_3_bottom -o Dardel/logs/branch_continuation_3_bottom.o -e Dardel/logs/branch_continuation_3_bottom.e Dardel/scripts/branch_continuation.sh 3 1 bottom
```

To then construct the proof witness running

``` julia
using Arblib, CGL

j, d = 3, 1
μ, γ, κ, ϵ, ξ₁, λ = CGL.sverak_params(Arb, j, d)

base_dirname = "Dardel/output/branch_continuation/"

part_filename_top = "branch_continuation_j=$(j)_d=$(d)_part=top.csv.gz"
part_filename_turn = "branch_continuation_j=$(j)_d=$(d)_part=turn.csv.gz"
part_filename_bottom = "branch_continuation_j=$(j)_d=$(d)_part=bottom.csv.gz"

dirnames = sort(readdir(base_dirname))

i_top = findlast(dirnames) do dirname
    in(part_filename_top, readdir(joinpath(base_dirname, dirname)))
end
i_turn = findlast(dirnames) do dirname
    in(part_filename_turn, readdir(joinpath(base_dirname, dirname)))
end
i_bottom = findlast(dirnames) do dirname
    in(part_filename_bottom, readdir(joinpath(base_dirname, dirname)))
end

filename_top = joinpath(base_dirname, dirnames[i_top], part_filename_top)
filename_turn = joinpath(base_dirname, dirnames[i_turn], part_filename_turn)
filename_bottom = joinpath(base_dirname, dirnames[i_bottom], part_filename_bottom)

parameters, data_top, data_turn, data_bottom =
    CGL.construct_proof_witness(filename_top, filename_turn, filename_bottom, ξ₁, λ);

CGL.check_proof_witness(parameters, data_top, data_turn, data_bottom)

directory = "proof/data/branch_j=$(j)_d=$(d)"

CGL.write_proof_witness(directory, parameters, data_top, data_turn, data_bottom)

# It can then be read with
#parameters, data_top, data_turn, data_bottom = CGL.read_proof_witness(directory)
```

## Other branches

### Existence of branches

``` shell
sbatch --time=0:45:00 -p main --job-name=branch_existence_1_top -o Dardel/logs/branch_existence_1_top.o -e Dardel/logs/branch_existence_1_top.e Dardel/scripts/branch_existence.sh 1 1 top

sbatch --time=1:00:00 -p main --job-name=branch_existence_2_top -o Dardel/logs/branch_existence_2_top.o -e Dardel/logs/branch_existence_2_top.e Dardel/scripts/branch_existence.sh 2 1 top

sbatch --time=1:00:00 -p main --job-name=branch_existence_4_top -o Dardel/logs/branch_existence_4_top.o -e Dardel/logs/branch_existence_4_top.e Dardel/scripts/branch_existence.sh 4 1 top

sbatch --time=1:00:00 -p main --job-name=branch_existence_5_top -o Dardel/logs/branch_existence_5_top.o -e Dardel/logs/branch_existence_5_top.e Dardel/scripts/branch_existence.sh 5 1 top

sbatch --time=0:45:00 -p main --job-name=branch_existence_6_top -o Dardel/logs/branch_existence_6_top.o -e Dardel/logs/branch_existence_6_top.e Dardel/scripts/branch_existence.sh 6 1 top

sbatch --time=0:45:00 -p main --job-name=branch_existence_7_top -o Dardel/logs/branch_existence_7_top.o -e Dardel/logs/branch_existence_7_top.e Dardel/scripts/branch_existence.sh 7 1 top

sbatch --time=3:00:00 -p main --job-name=branch_existence_8_top -o Dardel/logs/branch_existence_8_top.o -e Dardel/logs/branch_existence_8_top.e Dardel/scripts/branch_existence.sh 8 1 top
```

### Continuation of branches

``` shell
sbatch --time=0:30:00 -p main --job-name=branch_continuation_1_top -o Dardel/logs/branch_continuation_1_top.o -e Dardel/logs/branch_continuation_1_top.e Dardel/scripts/branch_continuation.sh 1 1 top

sbatch --time=0:30:00 -p main --job-name=branch_continuation_2_top -o Dardel/logs/branch_continuation_2_top.o -e Dardel/logs/branch_continuation_2_top.e Dardel/scripts/branch_continuation.sh 2 1 top

sbatch --time=0:30:00 -p main --job-name=branch_continuation_4_top -o Dardel/logs/branch_continuation_4_top.o -e Dardel/logs/branch_continuation_4_top.e Dardel/scripts/branch_continuation.sh 4 1 top

sbatch --time=0:30:00 -p main --job-name=branch_continuation_5_top -o Dardel/logs/branch_continuation_5_top.o -e Dardel/logs/branch_continuation_5_top.e Dardel/scripts/branch_continuation.sh 5 1 top

sbatch --time=0:30:00 -p main --job-name=branch_continuation_6_top -o Dardel/logs/branch_continuation_6_top.o -e Dardel/logs/branch_continuation_6_top.e Dardel/scripts/branch_continuation.sh 6 1 top

sbatch --time=0:30:00 -p main --job-name=branch_continuation_7_top -o Dardel/logs/branch_continuation_7_top.o -e Dardel/logs/branch_continuation_7_top.e Dardel/scripts/branch_continuation.sh 7 1 top

sbatch --time=0:30:00 -p main --job-name=branch_continuation_8_top -o Dardel/logs/branch_continuation_8_top.o -e Dardel/logs/branch_continuation_8_top.e Dardel/scripts/branch_continuation.sh 8 1 top
```

## Running on non-SLURM systems
It is also possible to run the computations on a non-SLURM system
(running Linux). In this case the number of workers and threads needs
to be set explicitly using the environmental variables `CGL_WORKERS`
and `CGL_THREADS` respectively. For example the following code could
be used on a system with 256 threads.

``` shell
# Points d = 1
CGL_WORKERS=32 CGL_THREADS=4 sh Dardel/scripts/branch_points.sh 1

# Points d = 3
CGL_WORKERS=32 CGL_THREADS=4 sh Dardel/scripts/branch_points.sh 3
```
