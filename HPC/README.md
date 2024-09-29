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
`Dardel/scripts/`. If you want to run these scripts yourself on a
Slurm system you will need to update a few parameters. The top part of
these scripts are given by (something similar to)

``` shell
# This needs to be updated to the account of your cluster
#SBATCH --account naiss2024-22-1038

# This needs to be updated based on your cluster
#SBATCH --partition main
#SBATCH --nodes 1
#SBATCH --ntasks 16
#SBATCH --cpus-per-task 16
```

You need to update the `--account` value to your Slurm account. In
addition you need to updated the `--partition`, `--nodes`, `--ntasks`
and `--cpus-per-task` based on the configuration of your system.

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
# d = 1, fix epsilon
sbatch -t 0:10:00 -J branch_points_1_fix_epsilon HPC/scripts/branch_points.sh 1 0

# d = 1, fix kappa
sbatch -t 0:10:00 -J branch_points_1_fix_kappa HPC/scripts/branch_points.sh 1 1

# d = 3, fix epsilon
sbatch -t 0:10:00 -J branch_points_3_fix_epsilon HPC/scripts/branch_points.sh 3 0

# d = 3, fix kappa
sbatch -t 0:10:00 -J branch_points_3_fix_kappa HPC/scripts/branch_points.sh 3 1
```

### Counting critical points

The case `d = 1`
``` shell
# 4 min
sbatch -t 0:10:00 -J branch_critical_points_1_1_points HPC/scripts/branch_critical_points.sh 1 1 points

# 4 min
sbatch -t 0:10:00 -J branch_critical_points_2_1_points HPC/scripts/branch_critical_points.sh 2 1 points

# 4 min
sbatch -t 0:10:00 -J branch_critical_points_3_1_points HPC/scripts/branch_critical_points.sh 3 1 points

# 4 min
sbatch -t 0:10:00 -J branch_critical_points_4_1_points HPC/scripts/branch_critical_points.sh 4 1 points

# 4 min
sbatch -t 0:10:00 -J branch_critical_points_5_1_points HPC/scripts/branch_critical_points.sh 5 1 points

# 4 min
sbatch -t 0:10:00 -J branch_critical_points_6_1_points HPC/scripts/branch_critical_points.sh 6 1 points

# 4 min
sbatch -t 0:10:00 -J branch_critical_points_7_1_points HPC/scripts/branch_critical_points.sh 7 1 points

# 4 min
sbatch -t 0:10:00 -J branch_critical_points_8_1_points HPC/scripts/branch_critical_points.sh 8 1 points
```

The case `d = 3`
``` shell
# 4 min
sbatch -t 0:10:00 -J branch_critical_points_1_3_points HPC/scripts/branch_critical_points.sh 1 3 points

# 4 min
sbatch -t 0:10:00 -J branch_critical_points_2_3_points HPC/scripts/branch_critical_points.sh 2 3 points

# 4 min
sbatch -t 0:10:00 -J branch_critical_points_3_3_points HPC/scripts/branch_critical_points.sh 3 3 points

# 4 min
sbatch -t 0:10:00 -J branch_critical_points_4_3_points HPC/scripts/branch_critical_points.sh 4 3 points

# 4 min
sbatch -t 0:10:00 -J branch_critical_points_5_3_points HPC/scripts/branch_critical_points.sh 5 3 points
```

## Branches for $d = 1$

### Existence of branches

Top
``` shell
# 5 min
sbatch -t 0:10:00 -J branch_existence_1_1_top HPC/scripts/branch_existence.sh 1 1 top

# 5 min
sbatch -t 0:10:00 -J branch_existence_2_1_top HPC/scripts/branch_existence.sh 2 1 top

# 5 min
sbatch -t 0:10:00 -J branch_existence_3_1_top HPC/scripts/branch_existence.sh 3 1 top

# 10 min
sbatch -t 0:15:00 -J branch_existence_4_1_top HPC/scripts/branch_existence.sh 4 1 top

# 10 min
sbatch -t 0:15:00 -J branch_existence_5_1_top HPC/scripts/branch_existence.sh 5 1 top

# 15 min
sbatch -t 0:20:00 -J branch_existence_6_1_top HPC/scripts/branch_existence.sh 6 1 top

# 20 min
sbatch -t 0:30:00 -J branch_existence_7_1_top HPC/scripts/branch_existence.sh 7 1 top

# 60 min
sbatch -t 1:15:00 -J branch_existence_8_1_top HPC/scripts/branch_existence.sh 8 1 top
```

Turn
``` shell
# 15 min
sbatch -t 0:20:00 -J branch_existence_1_1_turn HPC/scripts/branch_existence.sh 1 1 turn

# 10 min
sbatch -t 0:15:00 -J branch_existence_2_1_turn HPC/scripts/branch_existence.sh 2 1 turn

# 15 min
sbatch -t 0:20:00 -J branch_existence_3_1_turn HPC/scripts/branch_existence.sh 3 1 turn

# 30 min
sbatch -t 00:40:00 -J branch_existence_4_1_turn HPC/scripts/branch_existence.sh 4 1 turn

# 60 min
sbatch -t 1:15:00 -J branch_existence_5_1_turn HPC/scripts/branch_existence.sh 5 1 turn

 # 2 h 10 min
sbatch -t 2:30:00 -J branch_existence_6_1_turn HPC/scripts/branch_existence.sh 6 1 turn

 # 32% in 3 h
sbatch -t 3:00:00 -J branch_existence_7_1_turn HPC/scripts/branch_existence.sh 7 1 turn

# 2% in 3 h
sbatch -t 3:00:00 -J branch_existence_8_1_turn HPC/scripts/branch_existence.sh 8 1 turn
```

Bottom
``` shell
# 77% in 2 hours
# 83% in 3:40
# 85% in 4 hours
# 91% in 8 hours
# 25 min for 80 segments
sbatch -t 0:30:00 -J branch_existence_1_1_bottom HPC/scripts/branch_existence.sh 1 1 bottom 80

# 300 min for full
# 20 min for 80 segments
sbatch -t 5:30:00 -J branch_existence_2_1_bottom HPC/scripts/branch_existence.sh 2 1 bottom 80

# 195 min
sbatch -t 3:30:00 -J branch_existence_3_1_bottom HPC/scripts/branch_existence.sh 3 1 bottom

# 145 min
sbatch -t 3:00:00 -J branch_existence_4_1_bottom HPC/scripts/branch_existence.sh 4 1 bottom

# 20 min for 8 segments
sbatch -t 2:00:00 -J branch_existence_5_1_bottom HPC/scripts/branch_existence.sh 5 1 bottom 8

# 55% in 2 hours
# 35 min for 16 segments
sbatch -t 0:45:00 -J branch_existence_6_1_bottom HPC/scripts/branch_existence.sh 6 1 bottom 16

# 0% in 1 hour
sbatch -t 1:00:00 -J branch_existence_7_1_bottom HPC/scripts/branch_existence.sh 7 1 bottom

# 0% in 1 hour
sbatch -t 1:00:00 -J branch_existence_8_1_bottom HPC/scripts/branch_existence.sh 8 1 bottom
```

### Continuation of branches

Top
``` shell
# 2 min
sbatch -t 0:10:00 -J branch_continuation_1_1_top HPC/scripts/branch_continuation.sh 1 1 top

# 2 min
sbatch -t 0:10:00 -J branch_continuation_2_1_top HPC/scripts/branch_continuation.sh 2 1 top

# 3 min
sbatch -t 0:10:00 -J branch_continuation_3_1_top HPC/scripts/branch_continuation.sh 3 1 top

# 2 min
sbatch -t 0:10:00 -J branch_continuation_4_1_top HPC/scripts/branch_continuation.sh 4 1 top

# 2 min
sbatch -t 0:10:00 -J branch_continuation_5_1_top HPC/scripts/branch_continuation.sh 5 1 top

# 2 min
sbatch -t 0:10:00 -J branch_continuation_6_1_top HPC/scripts/branch_continuation.sh 6 1 top

# 2 min
sbatch -t 0:10:00 -J branch_continuation_7_1_top HPC/scripts/branch_continuation.sh 7 1 top

# 2 min
sbatch -t 0:10:00 -J branch_continuation_8_1_top HPC/scripts/branch_continuation.sh 8 1 top
```

Turn
``` shell
# 6 min
sbatch -t 0:10:00 -J branch_continuation_1_1_turn HPC/scripts/branch_continuation.sh 1 1 turn

# 2 min
sbatch -t 0:10:00 -J branch_continuation_2_1_turn HPC/scripts/branch_continuation.sh 2 1 turn

# 3 min
sbatch -t 0:10:00 -J branch_continuation_3_1_turn HPC/scripts/branch_continuation.sh 3 1 turn

# 2 min
sbatch -t 0:10:00 -J branch_continuation_4_1_turn HPC/scripts/branch_continuation.sh 4 1 turn

# 2 min
sbatch -t 0:10:00 -J branch_continuation_5_1_turn HPC/scripts/branch_continuation.sh 5 1 turn

# 2 min
sbatch -t 0:10:00 -J branch_continuation_6_1_turn HPC/scripts/branch_continuation.sh 6 1 turn

# TODO
#sbatch -t 0:10:00 -J branch_continuation_7_1_turn HPC/scripts/branch_continuation.sh 7 1 turn

# TODO
#sbatch -t 0:10:00 -J branch_continuation_8_1_turn HPC/scripts/branch_continuation.sh 8 1 turn
```

Bottom
``` shell
# 15 min for 80 segments
sbatch -t 0:25:00 -J branch_continuation_1_1_bottom HPC/scripts/branch_continuation.sh 1 1 bottom

# 5 min for 80 segments
sbatch -t 0:10:00 -J branch_continuation_2_1_bottom HPC/scripts/branch_continuation.sh 2 1 bottom

# 5 min
#sbatch -t 0:10:00 -J branch_continuation_3_1_bottom HPC/scripts/branch_continuation.sh 3 1 bottom

# 5 min
sbatch -t 0:10:00 -J branch_continuation_4_1_bottom HPC/scripts/branch_continuation.sh 4 1 bottom

# 5 min for 8 segments
sbatch -t 0:10:00 -J branch_continuation_5_1_bottom HPC/scripts/branch_continuation.sh 5 1 bottom

# 5 min for 16 segments
sbatch -t 0:10:00 -J branch_continuation_6_1_bottom HPC/scripts/branch_continuation.sh 6 1 bottom

# TODO
#sbatch -t 0:10:00 -J branch_continuation_7_1_bottom HPC/scripts/branch_continuation.sh 7 1 bottom

# TODO
#sbatch -t 0:10:00 -J branch_continuation_8_1_bottom HPC/scripts/branch_continuation.sh 8 1 bottom
```

### Counting critical points

Top
``` shell
# 2 min
sbatch -t 1:00:00 -J branch_critical_points_1_1_top HPC/scripts/branch_critical_points.sh 1 1 top

# 3 min
sbatch -t 0:10:00 -J branch_critical_points_2_1_top HPC/scripts/branch_critical_points.sh 2 1 top

# 3 min
sbatch -t 0:20:00 -J branch_critical_points_3_1_top HPC/scripts/branch_critical_points.sh 3 1 top

# 2 min
sbatch -t 0:10:00 -J branch_critical_points_4_1_top HPC/scripts/branch_critical_points.sh 4 1 top

# 3 min
sbatch -t 0:10:00 -J branch_critical_points_5_1_top HPC/scripts/branch_critical_points.sh 5 1 top

# 5 min
sbatch -t 0:10:00 -J branch_critical_points_6_1_top HPC/scripts/branch_critical_points.sh 6 1 top

# 10 min
sbatch -t 0:20:00 -J branch_critical_points_7_1_top HPC/scripts/branch_critical_points.sh 7 1 top

# 163 min
sbatch -t 3:00:00 -J branch_critical_points_8_1_top HPC/scripts/branch_critical_points.sh 8 1 top
```

Turn
``` shell
# 40% in 60 min
sbatch -t 4:00:00 -J branch_critical_points_1_1_turn HPC/scripts/branch_critical_points.sh 1 1 turn

# 45 min
sbatch -t 1:00:00 -J branch_critical_points_2_1_turn HPC/scripts/branch_critical_points.sh 2 1 turn

# 18 min
sbatch -t 1:00:00 -J branch_critical_points_3_1_turn HPC/scripts/branch_critical_points.sh 3 1 turn

# 72% in 60 min
sbatch -t 1:00:00 -J branch_critical_points_4_1_turn HPC/scripts/branch_critical_points.sh 4 1 turn

# 27% in 60 min
sbatch -t 1:00:00 -J branch_critical_points_5_1_turn HPC/scripts/branch_critical_points.sh 5 1 turn

# 14% in 60 min
sbatch -t 1:00:00 -J branch_critical_points_6_1_turn HPC/scripts/branch_critical_points.sh 6 1 turn

# TODO
#sbatch -t 1:00:00 -J branch_critical_points_7_1_turn HPC/scripts/branch_critical_points.sh 7 1 turn

# TODO
#sbatch -t 1:00:00 -J branch_critical_points_8_1_turn HPC/scripts/branch_critical_points.sh 8 1 turn
```

Bottom
``` shell
# 39% in 60 min
sbatch -t 4:00:00 -J branch_critical_points_1_1_bottom HPC/scripts/branch_critical_points.sh 1 1 bottom

# 45 min
sbatch -t 1:00:00 -J branch_critical_points_2_1_bottom HPC/scripts/branch_critical_points.sh 2 1 bottom

# TODO
#sbatch -t 1:00:00 -J branch_critical_points_3_1_bottom HPC/scripts/branch_critical_points.sh 3 1 bottom

# 72% in 60 min
sbatch -t 1:00:00 -J branch_critical_points_4_1_bottom HPC/scripts/branch_critical_points.sh 4 1 bottom

# 27% in 60 min
sbatch -t 1:00:00 -J branch_critical_points_5_1_bottom HPC/scripts/branch_critical_points.sh 5 1 bottom

# 14% in 60 min
sbatch -t 1:00:00 -J branch_critical_points_6_1_bottom HPC/scripts/branch_critical_points.sh 6 1 bottom

# TODO
#sbatch -t 1:00:00 -J branch_critical_points_7_1_bottom HPC/scripts/branch_critical_points.sh 7 1 bottom

# TODO
#sbatch -t 1:00:00 -J branch_critical_points_8_1_bottom HPC/scripts/branch_critical_points.sh 8 1 bottom
```

### Proof witness
Proof witnesses can then be constructed with

``` julia
using Arblib, CGL

setprecision(Arb, 128)

d = 1

for j = 1:8
    @info "j = $j"
    parameters, data_top, data_turn, data_bottom, data_connection_points =
        CGL.construct_proof_witness(j, d)

    CGL.check_proof_witness(
        parameters,
        data_top,
        data_turn,
        data_bottom,
        data_connection_points,
    )

    directory = joinpath(dirname(pathof(CGL)), "../proof/data/branch_j=$(j)_d=$(d)")

    CGL.write_proof_witness(
        directory,
        parameters,
        data_top,
        data_turn,
        data_bottom,
        data_connection_points,
    )
end
```

## Branches for $d = 3$

Top
``` shell
#sbatch -t 0:15:00 -J branch_existence_3_1_top HPC/scripts/branch_existence.sh 1 3 top

#sbatch -t 0:15:00 -J branch_existence_3_2_top HPC/scripts/branch_existence.sh 2 3 top

#sbatch -t 0:15:00 -J branch_existence_3_3_top HPC/scripts/branch_existence.sh 3 3 top

#sbatch -t 0:15:00 -J branch_existence_3_4_top HPC/scripts/branch_existence.sh 4 3 top

#sbatch -t 0:15:00 -J branch_existence_3_5_top HPC/scripts/branch_existence.sh 5 3 top
```

Turn
``` shell
# 22% in 60 min
sbatch -t 1:00:00 -J branch_existence_3_1_turn HPC/scripts/branch_existence.sh 1 3 turn 800:1301

# 20 min
sbatch -t 1:00:00 -J branch_existence_3_2_turn HPC/scripts/branch_existence.sh 2 3 turn 900:908

# 44% in 60 min
sbatch -t 1:00:00 -J branch_existence_3_3_turn HPC/scripts/branch_existence.sh 3 3 turn 1500:2001

# 62% in 60 min
sbatch -t 1:00:00 -J branch_existence_3_4_turn HPC/scripts/branch_existence.sh 4 3 turn 1100:1108

#
#sbatch -t 1:00:00 -J branch_existence_3_5_turn HPC/scripts/branch_existence.sh 5 3 turn 1300:1308
```

Bottom
``` shell
# 40 min from 1300
sbatch -t 1:00:00 -J branch_existence_3_1_bottom HPC/scripts/branch_existence.sh 1 3 bottom 1300:2433

#sbatch -t 1:00:00 -J branch_existence_3_2_bottom HPC/scripts/branch_existence.sh 2 3 bottom

# 9% in 60 min for 2000:11043
# 33 min for 2000:2500
sbatch -t 1:00:00 -J branch_existence_3_3_bottom HPC/scripts/branch_existence.sh 3 3 bottom 2000:2500

#sbatch -t 1:00:00 -J branch_existence_3_4_bottom HPC/scripts/branch_existence.sh 4 3 bottom

# 0% in 60 min
#sbatch -t 1:00:00 -J branch_existence_3_5_bottom HPC/scripts/branch_existence.sh 5 3 bottom 1271:3000
```

## Running on Agate
- `cpus-per-task` should be half of Dardel
- `CGL_SLURM_MEM_PER_NODE` should be set to a larger value

## Running on non-SLURM systems
It is also possible to run the computations on a non-SLURM system
(running Linux). In this case the number of workers and threads needs
to be set explicitly using the environmental variables `CGL_WORKERS`
and `CGL_THREADS` respectively. For example the following code could
be used on a system with 256 threads.

``` shell
# Points d = 1
CGL_WORKERS=32 CGL_THREADS=4 sh HPC/scripts/branch_points.sh 1

# Points d = 3
CGL_WORKERS=32 CGL_THREADS=4 sh HPC/scripts/branch_points.sh 3
```
