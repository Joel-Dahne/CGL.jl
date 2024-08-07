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

## Pointwise verification of branches
This verifies the existence of a solution for all points that make up
the numerical approximations of the branches. The `branch_points.sh`
script takes 4 arguments:
1. `d = 1`: The $d$ value for which to compute branches. For $d = 1$ it
   computes branches 1 to 8, and for $d = 3$ it computes branches 1 to
   5.
2. `fix_kappa = 0`: Should be either `0` or `1`. With the value `0`
   the verification is done by fixing $\epsilon$, with the value `1`
   the verification is done by fixing $\kappa$.
3. `N = 0`: Number of points to verify, the value `0` indicates that
   all points should be verified. This is mainly used in development
   to allow for running a shorter run.
4. `scaling = 1.0`: Optionally scale the values using the equations
   scaling symmetry. This maps `ω` to `ω * scaling^2`, with the
   corresponding scaling for the other values. This is mainly used as
   a sanity test to see that it returns compatible results for
   different scalings.

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

setprecision(Arb, 128)

j, d = 3, 1

parameters, data_top, data_turn, data_bottom, data_connection_points =
    CGL.construct_proof_witness(j, d);

CGL.check_proof_witness(parameters, data_top, data_turn, data_bottom, data_connection_points)

directory = joinpath(dirname(pathof(CGL)), "../proof/data/branch_j=$(j)_d=$(d)")

CGL.write_proof_witness(directory, parameters, data_top, data_turn, data_bottom, data_connection_points)

# It can then be read with
#parameters, data_top, data_turn, data_bottom, data_connection_points =
#    CGL.read_proof_witness(directory)
```

## Other branches

### Existence of branches

Top
``` shell
# 10 min
sbatch --time=0:15:00 -p main --job-name=branch_existence_1_top -o Dardel/logs/branch_existence_1_top.o -e Dardel/logs/branch_existence_1_top.e Dardel/scripts/branch_existence.sh 1 1 top

# 10 min
sbatch --time=0:15:00 -p main --job-name=branch_existence_2_top -o Dardel/logs/branch_existence_2_top.o -e Dardel/logs/branch_existence_2_top.e Dardel/scripts/branch_existence.sh 2 1 top

# 10 min
#sbatch --time=0:15:00 -p main --job-name=branch_existence_3_top -o Dardel/logs/branch_existence_3_top.o -e Dardel/logs/branch_existence_3_top.e Dardel/scripts/branch_existence.sh 3 1 top # Done above

# 15 min
sbatch --time=0:20:00 -p main --job-name=branch_existence_4_top -o Dardel/logs/branch_existence_4_top.o -e Dardel/logs/branch_existence_4_top.e Dardel/scripts/branch_existence.sh 4 1 top

# 20 min
sbatch --time=0:30:00 -p main --job-name=branch_existence_5_top -o Dardel/logs/branch_existence_5_top.o -e Dardel/logs/branch_existence_5_top.e Dardel/scripts/branch_existence.sh 5 1 top

# 45 min
sbatch --time=1:00:00 -p main --job-name=branch_existence_6_top -o Dardel/logs/branch_existence_6_top.o -e Dardel/logs/branch_existence_6_top.e Dardel/scripts/branch_existence.sh 6 1 top

# 45 min
sbatch --time=1:00:00 -p main --job-name=branch_existence_7_top -o Dardel/logs/branch_existence_7_top.o -e Dardel/logs/branch_existence_7_top.e Dardel/scripts/branch_existence.sh 7 1 top

# 2 h 25 min
sbatch --time=3:00:00 -p main --job-name=branch_existence_8_top -o Dardel/logs/branch_existence_8_top.o -e Dardel/logs/branch_existence_8_top.e Dardel/scripts/branch_existence.sh 8 1 top
```

Turn
``` shell
# 80 min
sbatch --time=1:30:00 -p main --job-name=branch_existence_1_turn -o Dardel/logs/branch_existence_1_turn.o -e Dardel/logs/branch_existence_1_turn.e Dardel/scripts/branch_existence.sh 1 1 turn

# 12 min
sbatch --time=0:15:00 -p main --job-name=branch_existence_2_turn -o Dardel/logs/branch_existence_2_turn.o -e Dardel/logs/branch_existence_2_turn.e Dardel/scripts/branch_existence.sh 2 1 turn

# 20 min
#sbatch --time=0:25:00 -p main --job-name=branch_existence_3_turn -o Dardel/logs/branch_existence_3_turn.o -e Dardel/logs/branch_existence_3_turn.e Dardel/scripts/branch_existence.sh 3 1 turn # Done above

# 35 min
sbatch --time=00:45:00 -p main --job-name=branch_existence_4_turn -o Dardel/logs/branch_existence_4_turn.o -e Dardel/logs/branch_existence_4_turn.e Dardel/scripts/branch_existence.sh 4 1 turn

# 90 min
sbatch --time=1:40:00 -p main --job-name=branch_existence_5_turn -o Dardel/logs/branch_existence_5_turn.o -e Dardel/logs/branch_existence_5_turn.e Dardel/scripts/branch_existence.sh 5 1 turn

 # TODO
sbatch --time=3:00:00 -p main --job-name=branch_existence_6_turn -o Dardel/logs/branch_existence_6_turn.o -e Dardel/logs/branch_existence_6_turn.e Dardel/scripts/branch_existence.sh 6 1 turn

 # TODO
sbatch --time=3:00:00 -p main --job-name=branch_existence_7_turn -o Dardel/logs/branch_existence_7_turn.o -e Dardel/logs/branch_existence_7_turn.e Dardel/scripts/branch_existence.sh 7 1 turn

 # TODO: Long time
sbatch --time=0:30:00 -p main --job-name=branch_existence_8_turn -o Dardel/logs/branch_existence_8_turn.o -e Dardel/logs/branch_existence_8_turn.e Dardel/scripts/branch_existence.sh 8 1 turn
```

Bottom
``` shell
# 77% in 2 hours
# 83% in 3:40
# 85% in 4 hours
# 91% in 8 hours
# 22 min for 80 segments
sbatch --time=0:30:00 -p main --job-name=branch_existence_1_bottom -o Dardel/logs/branch_existence_1_bottom.o -e Dardel/logs/branch_existence_1_bottom.e Dardel/scripts/branch_existence.sh 1 1 bottom 80

# 285 min
sbatch --time=5:00:00 -p main --job-name=branch_existence_2_bottom -o Dardel/logs/branch_existence_2_bottom.o -e Dardel/logs/branch_existence_2_bottom.e Dardel/scripts/branch_existence.sh 2 1 bottom

# 160 min
sbatch --time=3:00:00 -p main --job-name=branch_existence_3_bottom -o Dardel/logs/branch_existence_3_bottom.o -e Dardel/logs/branch_existence_3_bottom.e Dardel/scripts/branch_existence.sh 3 1 bottom

# 200 min
sbatch --time=4:00:00 -p main --job-name=branch_existence_4_bottom -o Dardel/logs/branch_existence_4_bottom.o -e Dardel/logs/branch_existence_4_bottom.e Dardel/scripts/branch_existence.sh 4 1 bottom

# TODO: Full estimate > 4 hours
# > 120 min for 14 segments
# 45 min for 8 segments
sbatch --time=2:00:00 -p main --job-name=branch_existence_5_bottom -o Dardel/logs/branch_existence_5_bottom.o -e Dardel/logs/branch_existence_5_bottom.e Dardel/scripts/branch_existence.sh 5 1 bottom 8

# TODO: Estimate > 8 hours
sbatch --time=2:00:00 -p main --job-name=branch_existence_6_bottom -o Dardel/logs/branch_existence_6_bottom.o -e Dardel/logs/branch_existence_6_bottom.e Dardel/scripts/branch_existence.sh 6 1 bottom

# TODO: Estimate > 24 hours
sbatch --time=2:00:00 -p main --job-name=branch_existence_7_bottom -o Dardel/logs/branch_existence_7_bottom.o -e Dardel/logs/branch_existence_7_bottom.e Dardel/scripts/branch_existence.sh 7 1 bottom

# TODO: Estimate > 24 hours
sbatch --time=2:00:00 -p main --job-name=branch_existence_8_bottom -o Dardel/logs/branch_existence_8_bottom.o -e Dardel/logs/branch_existence_8_bottom.e Dardel/scripts/branch_existence.sh 8 1 bottom
```

### Continuation of branches

Top
``` shell
# 5 min
sbatch --time=0:10:00 -p main --job-name=branch_continuation_1_top -o Dardel/logs/branch_continuation_1_top.o -e Dardel/logs/branch_continuation_1_top.e Dardel/scripts/branch_continuation.sh 1 1 top

# 4 min
sbatch --time=0:10:00 -p main --job-name=branch_continuation_2_top -o Dardel/logs/branch_continuation_2_top.o -e Dardel/logs/branch_continuation_2_top.e Dardel/scripts/branch_continuation.sh 2 1 top

# 2 min
sbatch --time=0:10:00 -p main --job-name=branch_continuation_3_top -o Dardel/logs/branch_continuation_3_top.o -e Dardel/logs/branch_continuation_3_top.e Dardel/scripts/branch_continuation.sh 3 1 top

# 4 min
sbatch --time=0:10:00 -p main --job-name=branch_continuation_4_top -o Dardel/logs/branch_continuation_4_top.o -e Dardel/logs/branch_continuation_4_top.e Dardel/scripts/branch_continuation.sh 4 1 top

# 5 min
sbatch --time=0:10:00 -p main --job-name=branch_continuation_5_top -o Dardel/logs/branch_continuation_5_top.o -e Dardel/logs/branch_continuation_5_top.e Dardel/scripts/branch_continuation.sh 5 1 top

# 4 min
sbatch --time=0:10:00 -p main --job-name=branch_continuation_6_top -o Dardel/logs/branch_continuation_6_top.o -e Dardel/logs/branch_continuation_6_top.e Dardel/scripts/branch_continuation.sh 6 1 top

# 6 min
sbatch --time=0:10:00 -p main --job-name=branch_continuation_7_top -o Dardel/logs/branch_continuation_7_top.o -e Dardel/logs/branch_continuation_7_top.e Dardel/scripts/branch_continuation.sh 7 1 top

# 7 min
sbatch --time=0:10:00 -p main --job-name=branch_continuation_8_top -o Dardel/logs/branch_continuation_8_top.o -e Dardel/logs/branch_continuation_8_top.e Dardel/scripts/branch_continuation.sh 8 1 top
```

Turn
``` shell
# 30 min
sbatch --time=0:40:00 -p main --job-name=branch_continuation_1_turn -o Dardel/logs/branch_continuation_1_turn.o -e Dardel/logs/branch_continuation_1_turn.e Dardel/scripts/branch_continuation.sh 1 1 turn

# 5 min
sbatch --time=0:10:00 -p main --job-name=branch_continuation_2_turn -o Dardel/logs/branch_continuation_2_turn.o -e Dardel/logs/branch_continuation_2_turn.e Dardel/scripts/branch_continuation.sh 2 1 turn

# 3 min
sbatch --time=0:10:00 -p main --job-name=branch_continuation_3_turn -o Dardel/logs/branch_continuation_3_turn.o -e Dardel/logs/branch_continuation_3_turn.e Dardel/scripts/branch_continuation.sh 3 1 turn

# 6 min
sbatch --time=0:10:00 -p main --job-name=branch_continuation_4_turn -o Dardel/logs/branch_continuation_4_turn.o -e Dardel/logs/branch_continuation_4_turn.e Dardel/scripts/branch_continuation.sh 4 1 turn

# 20 min
sbatch --time=0:30:00 -p main --job-name=branch_continuation_5_turn -o Dardel/logs/branch_continuation_5_turn.o -e Dardel/logs/branch_continuation_5_turn.e Dardel/scripts/branch_continuation.sh 5 1 turn

# TODO
sbatch --time=0:10:00 -p main --job-name=branch_continuation_6_turn -o Dardel/logs/branch_continuation_6_turn.o -e Dardel/logs/branch_continuation_6_turn.e Dardel/scripts/branch_continuation.sh 6 1 turn

# TODO
sbatch --time=0:10:00 -p main --job-name=branch_continuation_7_turn -o Dardel/logs/branch_continuation_7_turn.o -e Dardel/logs/branch_continuation_7_turn.e Dardel/scripts/branch_continuation.sh 7 1 turn

# TODO
sbatch --time=0:10:00 -p main --job-name=branch_continuation_8_turn -o Dardel/logs/branch_continuation_8_turn.o -e Dardel/logs/branch_continuation_8_turn.e Dardel/scripts/branch_continuation.sh 8 1 turn
```

Bottom
``` shell
# 5 min for 80 segments
sbatch --time=2:00:00 -p main --job-name=branch_continuation_1_bottom -o Dardel/logs/branch_continuation_1_bottom.o -e Dardel/logs/branch_continuation_1_bottom.e Dardel/scripts/branch_continuation.sh 1 1 bottom

# 60 min
sbatch --time=1:20:00 -p main --job-name=branch_continuation_2_bottom -o Dardel/logs/branch_continuation_2_bottom.o -e Dardel/logs/branch_continuation_2_bottom.e Dardel/scripts/branch_continuation.sh 2 1 bottom

# 20 min
sbatch --time=0:30:00 -p main --job-name=branch_continuation_3_bottom -o Dardel/logs/branch_continuation_3_bottom.o -e Dardel/logs/branch_continuation_3_bottom.e Dardel/scripts/branch_continuation.sh 3 1 bottom

# 35 min
sbatch --time=0:40:00 -p main --job-name=branch_continuation_4_bottom -o Dardel/logs/branch_continuation_4_bottom.o -e Dardel/logs/branch_continuation_4_bottom.e Dardel/scripts/branch_continuation.sh 4 1 bottom

# 7 min for 8 segments
sbatch --time=2:00:00 -p main --job-name=branch_continuation_5_bottom -o Dardel/logs/branch_continuation_5_bottom.o -e Dardel/logs/branch_continuation_5_bottom.e Dardel/scripts/branch_continuation.sh 5 1 bottom

# TODO
sbatch --time=0:10:00 -p main --job-name=branch_continuation_6_bottom -o Dardel/logs/branch_continuation_6_bottom.o -e Dardel/logs/branch_continuation_6_bottom.e Dardel/scripts/branch_continuation.sh 6 1 bottom

# TODO
sbatch --time=0:10:00 -p main --job-name=branch_continuation_7_bottom -o Dardel/logs/branch_continuation_7_bottom.o -e Dardel/logs/branch_continuation_7_bottom.e Dardel/scripts/branch_continuation.sh 7 1 bottom

# TODO
sbatch --time=0:10:00 -p main --job-name=branch_continuation_8_bottom -o Dardel/logs/branch_continuation_8_bottom.o -e Dardel/logs/branch_continuation_8_bottom.e Dardel/scripts/branch_continuation.sh 8 1 bottom
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
