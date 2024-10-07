# HPC computations
The computations that are done to prove the existence of branches of
solutions to the CGL equation require a significant amount of
computational time, and have for that reason been done on an HPC
cluster. Two different clusters have been used,
[Dardel](https://www.pdc.kth.se/hpc-services/computing-systems/about-the-dardel-hpc-system-1.1053338)
and
[Agate](https://msi.umn.edu/about-msi-services/high-performance-computing/agate).
This directory contains the description of how these computations were
run.

For Dardel, the computations were enabled by resources provided by the
National Academic Infrastructure for Supercomputing in Sweden (NAISS)
and the Swedish National Infrastructure for Computing (SNIC) at the
PDC Center for High Performance Computing, KTH Royal Institute of
Technology, partially funded by the Swedish Research Council through
grant agreements no. 2022-06725 and no. 2018-05973. For Agate, the
computations were enabled by resources provided by Minnesota
Supercomputing Institute (MSI) at the University of Minnesota.

Both clusters makes use of Slurm for job scheduling. The scripts for
running the computations are written in terms of Slurm scripts. From a
fresh clone of this repository an initial setup is done with the
following commands, executed from the root of this repository.

``` shell
# Install all required packages for Julia - this takes some time
julia --project=. --eval 'using Pkg; Pkg.instantiate()'
```

The computations are run using several scripts found in
`HPC/scripts/`. If you want to run these scripts yourself on a Slurm
system you will need to update a few parameters. The top part of these
scripts are given by (something similar to)

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

## Overview of scripts
The [`script/`](/HPC/script/) directory contains four scripts, each
with one Julia script and one Slurm script, that are responsible for
different parts of the computations,
- `branch_points` - Handles pointwise verification along the branches.
- `branch_existence` - Handles local existence along the branches, but
  dynamically splitting them into small segments.
- `branch_continuation` - Takes the output of `branch_existence` and
  verifies that the segments can be joined together into one
  continuous curve, bisecting segments further if required.
- `branch_critical_points` - Takes the output of
   `branch_continuation.jl` and computes the number of critical points
   for each segment.

There is also a script `branch_points_verification`, that computes
pointwise enclosures using several different parameters and scalings
and verifies that all results agree.

The scripts are mostly thin wrappers around the functions found in the
[`src/orchestration/`](/src/orchestration) directory.

## Running scripts
The computations are split into three families
1. Pointwise verification of all branches
2. Full verification of branches for $d = 1$ (called **Case I** in the
   paper)
3. Full verification of branches for $d = 3$ (called **Case II** in
   the paper)

The data is written to subdirectories of `HPC/output/` and the logs
are found in `HPC/logs/`.

The instructions given below are intended for running the code on a
Slurm cluster. To run the computations on another system an adjusted
approach is required. In this case the number of workers and threads
needs to be set explicitly using the environmental variables
`CGL_WORKERS` and `CGL_THREADS` respectively. For example, code given
below for the pointwise verification could be run on a non-Slurm
cluster with 256 threads as

``` shell
# d = 1, fix epsilon
CGL_WORKERS=32 CGL_THREADS=4 sh HPC/scripts/branch_points.sh 1 0

# d = 3, fix epsilon
CGL_WORKERS=32 CGL_THREADS=4 sh HPC/scripts/branch_points.sh 3 0
```

### Pointwise verification of branches
This uses the [`branch_points.sh`](/HPC/scripts/branch_points.sh)
script and verifies the existence of a solution for all points that
make up the numerical approximations of the branches. The script takes
4 arguments:
1. `d = 1`: The $d$ value for which to compute branches. For $d = 1$ it
   computes branches 1 to 8, and for $d = 3$ it computes branches 1 to
   5.
2. `fix_kappa = 0`: Should be either `0` or `1`. With the value `0`
   the verification is done by fixing $\epsilon$, with the value `1`
   the verification is done by fixing $\kappa$.
3. `N = 0`: Number of points to verify, the value `0` indicates that
   all points should be verified. This option is mainly used in
   development to allow for testing shorter runs.
4. `scaling = 1.0`: Optionally scale the values using the equations
   scaling symmetry. This maps `ω` to `ω * scaling^2`, with the
   corresponding scaling for the other values. This can be used as a
   sanity test to see that it returns compatible results for different
   scalings.

We do the computations for $d = 1$ and $d = 3$, using both $fix_kappa
= 0$ and $fix_kappa = 1$. On a Slurm cluster the computations are
started using the following commands, ran from the root directory of
the repository.

``` shell
# d = 1, fix epsilon
sbatch -t 0:10:00 -J branch_points_1_fix_epsilon HPC/scripts/branch_points.sh 1 0

# d = 1, fix kappa
sbatch -t 0:10:00 -J branch_points_1_fix_kappa HPC/scripts/branch_points.sh 1 1

# d = 3, fix epsilon
sbatch -t 0:15:00 -J branch_points_3_fix_epsilon HPC/scripts/branch_points.sh 3 0

# d = 3, fix kappa
sbatch -t 0:15:00 -J branch_points_3_fix_kappa HPC/scripts/branch_points.sh 3 1
```

### Full verification of branches for $d = 1$
The verification consists of three parts:
1. Running `branch_existence` to get initial enclosures of the curves
2. Running `branch_continuation` to get enclosures that can be joined
   into on continuous curve.
3. Running `branch_critical_points` to count the number of critical
   points.

Each branch is also split into three parts, a top part, a turning part
and a bottom part. For each part the three above computations have to
be run. As a final step the top, turn and bottom parts are joined
together using the function `proof_witness`.

#### Existence of branches
For step 1 the following code is run.

For top parts
``` shell
# 5 min
sbatch -t 0:10:00 -J branch_existence_1_1_top HPC/scripts/branch_existence.sh 1 1 top

# 5 min
sbatch -t 0:10:00 -J branch_existence_1_2_top HPC/scripts/branch_existence.sh 2 1 top

# 5 min
sbatch -t 0:10:00 -J branch_existence_1_3_top HPC/scripts/branch_existence.sh 3 1 top

# 10 min
sbatch -t 0:15:00 -J branch_existence_1_4_top HPC/scripts/branch_existence.sh 4 1 top

# 10 min
sbatch -t 0:15:00 -J branch_existence_1_5_top HPC/scripts/branch_existence.sh 5 1 top

# 15 min
sbatch -t 0:20:00 -J branch_existence_1_6_top HPC/scripts/branch_existence.sh 6 1 top

# 20 min
sbatch -t 0:30:00 -J branch_existence_1_7_top HPC/scripts/branch_existence.sh 7 1 top

# 60 min
sbatch -t 1:15:00 -J branch_existence_1_8_top HPC/scripts/branch_existence.sh 8 1 top
```

For turning parts (note that $j = 7$ and $j = 8$ are excluded)
``` shell
# 15 min
sbatch -t 0:20:00 -J branch_existence_1_1_turn HPC/scripts/branch_existence.sh 1 1 turn

# 10 min
sbatch -t 0:15:00 -J branch_existence_1_2_turn HPC/scripts/branch_existence.sh 2 1 turn

# 15 min
sbatch -t 0:20:00 -J branch_existence_1_3_turn HPC/scripts/branch_existence.sh 3 1 turn

# 25 min
sbatch -t 00:40:00 -J branch_existence_1_4_turn HPC/scripts/branch_existence.sh 4 1 turn

# 60 min
sbatch -t 1:20:00 -J branch_existence_1_5_turn HPC/scripts/branch_existence.sh 5 1 turn

 # 30 min
sbatch -t 1:00:00 -J branch_existence_1_6_turn HPC/scripts/branch_existence.sh 6 1 turn 100
```

For bottom parts (note that $j = 6$, $j = 7$ and $j = 8$ are excluded)
``` shell
# 20 min for 80 segments
sbatch -t 0:35:00 -J branch_existence_1_1_bottom HPC/scripts/branch_existence.sh 1 1 bottom 80

# 20 min for 80 segments
sbatch -t 0:30:00 -J branch_existence_1_2_bottom HPC/scripts/branch_existence.sh 2 1 bottom 80

# 50 min for 40 segments
sbatch -t 3:30:00 -J branch_existence_1_3_bottom HPC/scripts/branch_existence.sh 3 1 bottom 40

# 20 min for 40 segments
sbatch -t 0:45:00 -J branch_existence_1_4_bottom HPC/scripts/branch_existence.sh 4 1 bottom 40

# 20 min for 8 segments
sbatch -t 2:00:00 -J branch_existence_1_5_bottom HPC/scripts/branch_existence.sh 5 1 bottom 8
```

#### Continuation of branches
For step 2 the following code is run.

For top part
``` shell
# 2 min
sbatch -t 0:10:00 -J branch_continuation_1_1_top HPC/scripts/branch_continuation.sh 1 1 top

# 2 min
sbatch -t 0:10:00 -J branch_continuation_1_2_top HPC/scripts/branch_continuation.sh 2 1 top

# 3 min
sbatch -t 0:10:00 -J branch_continuation_1_3_top HPC/scripts/branch_continuation.sh 3 1 top

# 2 min
sbatch -t 0:10:00 -J branch_continuation_1_4_top HPC/scripts/branch_continuation.sh 4 1 top

# 2 min
sbatch -t 0:10:00 -J branch_continuation_1_5_top HPC/scripts/branch_continuation.sh 5 1 top

# 2 min
sbatch -t 0:10:00 -J branch_continuation_1_6_top HPC/scripts/branch_continuation.sh 6 1 top

# 2 min
sbatch -t 0:10:00 -J branch_continuation_1_7_top HPC/scripts/branch_continuation.sh 7 1 top

# 2 min
sbatch -t 0:10:00 -J branch_continuation_1_8_top HPC/scripts/branch_continuation.sh 8 1 top
```

For turning parts (note that $j = 7$ and $j = 8$ are excluded)
``` shell
# 6 min
sbatch -t 0:10:00 -J branch_continuation_1_1_turn HPC/scripts/branch_continuation.sh 1 1 turn

# 2 min
sbatch -t 0:10:00 -J branch_continuation_1_2_turn HPC/scripts/branch_continuation.sh 2 1 turn

# 3 min
sbatch -t 0:10:00 -J branch_continuation_1_3_turn HPC/scripts/branch_continuation.sh 3 1 turn

# 2 min
sbatch -t 0:10:00 -J branch_continuation_1_4_turn HPC/scripts/branch_continuation.sh 4 1 turn

# 2 min
sbatch -t 0:10:00 -J branch_continuation_1_5_turn HPC/scripts/branch_continuation.sh 5 1 turn

# 2 min
sbatch -t 0:10:00 -J branch_continuation_1_6_turn HPC/scripts/branch_continuation.sh 6 1 turn
```

For bottom parts (note that $j = 6$, $j = 7$ and $j = 8$ are excluded)
``` shell
# 15 min for 80 segments
sbatch -t 0:25:00 -J branch_continuation_1_1_bottom HPC/scripts/branch_continuation.sh 1 1 bottom

# 5 min for 80 segments
sbatch -t 0:10:00 -J branch_continuation_1_2_bottom HPC/scripts/branch_continuation.sh 2 1 bottom

# 5 min
sbatch -t 0:10:00 -J branch_continuation_1_3_bottom HPC/scripts/branch_continuation.sh 3 1 bottom

# 2 min
sbatch -t 0:10:00 -J branch_continuation_1_4_bottom HPC/scripts/branch_continuation.sh 4 1 bottom

# 5 min for 8 segments
sbatch -t 0:10:00 -J branch_continuation_1_5_bottom HPC/scripts/branch_continuation.sh 5 1 bottom
```

#### Counting critical points
For step 3 the following code is run.

For top part
``` shell
# 2 min
sbatch -t 1:00:00 -J branch_critical_points_1_1_top HPC/scripts/branch_critical_points.sh 1 1 top

# 3 min
sbatch -t 0:10:00 -J branch_critical_points_1_2_top HPC/scripts/branch_critical_points.sh 2 1 top

# 3 min
sbatch -t 0:20:00 -J branch_critical_points_1_3_top HPC/scripts/branch_critical_points.sh 3 1 top

# 2 min
sbatch -t 0:10:00 -J branch_critical_points_1_4_top HPC/scripts/branch_critical_points.sh 4 1 top

# 3 min
sbatch -t 0:10:00 -J branch_critical_points_1_5_top HPC/scripts/branch_critical_points.sh 5 1 top

# 3 min
sbatch -t 0:10:00 -J branch_critical_points_1_6_top HPC/scripts/branch_critical_points.sh 6 1 top

# 5 min
sbatch -t 0:20:00 -J branch_critical_points_1_7_top HPC/scripts/branch_critical_points.sh 7 1 top

# 160 min
sbatch -t 3:00:00 -J branch_critical_points_1_8_top HPC/scripts/branch_critical_points.sh 8 1 top
```

For turning parts (note that $j = 7$ and $j = 8$ are excluded)
``` shell
# 398 min
sbatch -t 8:00:00 -J branch_critical_points_1_1_turn HPC/scripts/branch_critical_points.sh 1 1 turn

# 45 min
sbatch -t 1:00:00 -J branch_critical_points_1_2_turn HPC/scripts/branch_critical_points.sh 2 1 turn

# 18 min
sbatch -t 0:30:00 -J branch_critical_points_1_3_turn HPC/scripts/branch_critical_points.sh 3 1 turn

# 104 min
sbatch -t 2:00:00 -J branch_critical_points_1_4_turn HPC/scripts/branch_critical_points.sh 4 1 turn

# 548 min
sbatch -t 10:00:00 -J branch_critical_points_1_5_turn HPC/scripts/branch_critical_points.sh 5 1 turn

# 522 min
sbatch -t 12:00:00 -J branch_critical_points_1_6_turn HPC/scripts/branch_critical_points.sh 6 1 turn
```

For bottom parts (note that $j = 6$, $j = 7$ and $j = 8$ are excluded)
``` shell
# 82 min
sbatch -t 1:40:00 -J branch_critical_points_1_1_bottom HPC/scripts/branch_critical_points.sh 1 1 bottom

# 45 min
sbatch -t 1:00:00 -J branch_critical_points_1_2_bottom HPC/scripts/branch_critical_points.sh 2 1 bottom

# 35 min
sbatch -t 4:00:00 -J branch_critical_points_1_3_bottom HPC/scripts/branch_critical_points.sh 3 1 bottom

# 54 min
sbatch -t 1:10:00 -J branch_critical_points_1_4_bottom HPC/scripts/branch_critical_points.sh 4 1 bottom

# 85 min
sbatch -t 2:00:00 -J branch_critical_points_1_5_bottom HPC/scripts/branch_critical_points.sh 5 1 bottom
```

#### Proof witness
With all of the above computations done, the data for the top, turn
and bottom parts of the branches can be joined together by running the
following Julia code.

``` julia
using Arblib, CGL

setprecision(Arb, 128)

for j = 1:8
    @info "j = $j"
    parameters, data_top, data_turn, data_bottom, data_connection_points =
        CGL.construct_proof_witness(j, 1)

    CGL.check_proof_witness(
        parameters,
        data_top,
        data_turn,
        data_bottom,
        data_connection_points,
    )

    directory = joinpath(dirname(pathof(CGL)), "../proof/data/branch_d=1_j=$(j)")

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

### Full verification of branches for $d = 3$
The verification for $d = 3$ is similar to the above case, except that
the number of critical points is not computed. We also only verify a
very small part of the curves. The top part is not verified for any of
the curves, only the turn and bottom parts are considered.

#### Existence

For the turn we compute parts for $j = 1, 2, 3, 4$.
``` shell
# 205 min
sbatch -t 12:00:00 -J branch_existence_3_1_turn HPC/scripts/branch_existence.sh 1 3 turn 800:1301

# 20 min
sbatch -t 0:40:00 -J branch_existence_3_2_turn HPC/scripts/branch_existence.sh 2 3 turn 900:908

# 117 min
sbatch -t 6:00:00 -J branch_existence_3_3_turn HPC/scripts/branch_existence.sh 3 3 turn 1500:2001

# 62 min
sbatch -t 4:00:00 -J branch_existence_3_4_turn HPC/scripts/branch_existence.sh 4 3 turn 1100:1108
```

For the bottom we only compute parts for $j = 1$ and $j = 3.
``` shell
# 40 min
sbatch -t 1:00:00 -J branch_existence_3_1_bottom HPC/scripts/branch_existence.sh 1 3 bottom 1300:2433

# 288 min
sbatch -t 8:00:00 -J branch_existence_3_3_bottom HPC/scripts/branch_existence.sh 3 3 bottom 2000:3500
```

#### Continuation

Turn
``` shell
# 135 min
sbatch -t 3:00:00 -J branch_continuation_3_1_turn HPC/scripts/branch_continuation.sh 1 3 turn

# 8 min
sbatch -t 0:20:00 -J branch_continuation_3_2_turn HPC/scripts/branch_continuation.sh 2 3 turn

# 2 min
sbatch -t 0:10:00 -J branch_continuation_3_3_turn HPC/scripts/branch_continuation.sh 3 3 turn

# 30 min
sbatch -t 1:00:00 -J branch_continuation_3_4_turn HPC/scripts/branch_continuation.sh 4 3 turn
```

Bottom
``` shell
# 34 min
sbatch -t 0:50:00 -J branch_continuation_3_1_bottom HPC/scripts/branch_continuation.sh 1 3 bottom

# 20 min
sbatch -t 0:30:00 -J branch_continuation_3_3_bottom HPC/scripts/branch_continuation.sh 3 3 bottom
```

#### Proof witness
Proof witnesses can then be constructed with

``` julia
using Arblib, CGL

setprecision(Arb, 128)

for j = 1:4
    @info "j = $j"
    parameters, data_top, data_turn, data_bottom, data_connection_points =
        CGL.construct_proof_witness(j, 3)

    CGL.check_proof_witness(
        parameters,
        data_top,
        data_turn,
        data_bottom,
        data_connection_points,
    )

    directory = joinpath(dirname(pathof(CGL)), "../proof/data/branch_d=3_j=$(j)")

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
