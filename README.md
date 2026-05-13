# Self-similar singular solutions to the nonlinear Schrödinger and the complex Ginzburg-Landau equation

This repository contains the code for the computer assisted parts of
the proofs for the paper [Self-similar singular solutions to the
nonlinear Schrödinger and the complex Ginzburg-Landau equation](TODO).

The results of the paper are presented in 5 different
[Pluto.jl](https://plutojl.org/) notebooks, found in the
[`proof`](proof) directory. These notebooks are responsible for
generating all the numbers and figures that appear in the paper. It is
possible to view the results of notebooks without running any code by
opening the corresponding html-files found in the [`proof`](proof)
directory. They can be opened in any browser such as Firefox. The five
notebooks are:
- `NLS-example.jl` - This is the notebook version of Section 5.1 in
  the paper, and contains a detailed walk through of the proof the
  existence of a self-similar singular solution to the 3D cubic NLS
  equation. The notebook follows the same layout as Section 5.1, but
  includes a bit more details that are related to the evaluation of
  the code. All numbers and figures in Section 5.1 are from this
  notebook.
- `NLS.jl` - This notebook contains the full results for the NLS
  equation. Theorems 3.1, 3.2 and 3.3 as well as Table 1 and 2 in
  Section 3 are based on this notebook.
- `CGL-example.jl` - The CGL version of `NLS-example.jl`, it
  corresponds to Section 6.1 in the paper. All numbers and figures in
  Section 6.1 are from this notebook.
- `CGL.jl` - The CGL version of `NLS.jl`, it presents the full results
  for the CGL equation. The results of theorems 4.1 and 4.4, as well
  as all the figures in Sections 4 and 6.2 are from this notebook.
- `CGL-pointwise.jl` - This notebook gives pointwise results for the
  branches of the CGL equation in Case II. It corresponds to Section 9
  in the paper.

The example notebooks, `NLS-example.jl` and `CGL-example.jl`, contain
a fairly detailed description of the computations performed, generally
following the same structure as the corresponding sections in the
paper. The other notebooks are primarily concerned with generating the
results and corresponding figures, and contain much less prose.
Explanations of the data and figures they generate are found in the
paper.

## Reproducing the proofs

The proofs were generated with Julia version 1.12.6. This repository
contains the same `Manifest.toml` file as was used when running the
proofs, this allows installing exactly the same versions of the Julia
packages. Part of the computations are done using the C++ library
[CAPD](https://github.com/CAPDGroup/CAPD). There is no Julia package
wrapping the CAPD library and programs therefore have to be compiled
outside of Julia, which complicates the setup process. Information
about how to install CAPD and compile the required code is found in
[`CAPD/README.md`](CAPD/README.md).

Once CAPD has been installed and the required programs compiled, the
Julia part of the code can be setup by starting Julia from this
directory and running

``` julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

This will likely take some time the first time you run it since it has
to download and compile all the packages. You can see if it seems to
work by running `Pkg.test()`. If the tests run successfully then you
should be good to go!

### Reproducing proofs for the NLS equation

The proofs for the NLS equation require a relatively small amount of
computational time, and all of the computations are performed inside
the notebooks `NLS-example.jl` and `NLS.jl` (found in the `proof`
directory). To run the notebooks you first need to start Pluto,
starting Julia from this directory you can run

``` julia
using Pkg
Pkg.activate(".")
using Pluto
Pluto.run()
```

which should open a Pluto tab in your browser. Now you can open the
notebooks inside the `proof` directory through this and it should run
the proof. In the `NLS-example` notebook a lot of the intermediate
steps in the computations are shown and explained, for the `NLS`
notebook only final results are presented.

### Reproducing the proofs for the CGL equation

The proofs for the CGL equation require significantly more
computational time, and it would not be feasible to do all of the
computations inside of the notebooks. The computations are instead run
on a separate HPC system, the output of these computations is then
saved and stored in the `proof/data` directory. The notebooks read
this precomputed data and generate the figures and numbers that
appear in the paper. In the case of the `CGL-example.jl` notebook a
small part of the computations are also performed in the notebook, to
show what it looks like.

The procedure for generating the precomputed data is found in
[`HPC/README.md`](HPC/README.md).

## Notes about the implementation

The code in this repository is spread out over four directories:

1. `src/` - This directory contains the vast majority of the code
   (around 10000 lines of code, though a large part of that is
   documentation). More information about it is given below.
2. `CAPD/src/` - This directory contains the CAPD code, more
   information about it is found in
   [`CAPD/README.md`](CAPD/README.md).
3. `HPC/scripts/` - This directory contains the scripts that are used
   for running computations on an HPC system. It is primarily glue
   code and does not handle any actual computations.
4. `proof/` - As mentioned, this directory contains the notebooks. The
   code implemented here is primarily for presenting the generated
   results.

Some good information to have:
- Some of the parameters are fixed throughout the computations, these
  are $d$, $\omega$, $\sigma$ and $\delta$. Their values depend on if
  we are considering Case I or Case II in the paper (though the code
  supports other values than that). To simplify the code we store all
  of these in a struct, `CGLParams`, which is passed around under the
  name `Λ`.
- Many of the implemented functions contain documentation that
  explains what they do. These occasionally refer to equations or
  lemmas in the paper. For these references we use the notation "REF",
  e.g. "Lemma REF(lemma:U)". We use labels, e.g. `lemma:U`, rather
  than numbers to make it easier to keep the code and the paper in
  sync. The labels are not directly seen in the PDF version of the
  paper, but can be accessed by downloading the LaTeX code for the
  paper from ArXiv.
- Many of the functions have associated tests in [`test`](test) that
  serve to increase the confidence in the implementation. All of these
  tests can be run with `Pkg.test()` as discussed above.

The parts of the code responsible for generating the proofs are found
in `src/` and `CAPD/src/`. The CAPD code is, as mentioned, discussed
in [`CAPD/README.md`](CAPD/README.md). We here go through the code in
the [`src/`](src) directory. The directory consists of 5
subdirectories

1. [`CGLBranch/`](src/CGLBranch)
2. [`Q_zero/`](src/Q_zero)
3. [`Q_infinity/`](src/Q_infinity)
4. [`branch/`](src/branch)
5. [`orchestration/`](src/orchestration)

as well as a number of separate files found directly in [`src/`](src).
Let us go through the most important parts.

1. [`CGLBranch/`](src/CGLBranch) - This directory contains the
   [BifurcationKit.jl](https://github.com/bifurcationkit/BifurcationKit.jl)
   code used for computing numerical approximations of the branches
   and is completely self contained (it is implemented as a sub
   module). It supports computing the branches using both $\epsilon$
   and $\kappa$ as continuation parameters, as well as two different
   approaches for approximating $Q_\infty$. The paper only makes use
   of the continuation in $\epsilon$ and one of the approaches for
   approximating $Q_\infty$, the other versions were used in the
   development phase.
2. [`Q_zero/`](src/Q_zero) - This directory is related to Section 8 in
   the paper and contains all the Julia code used for evaluating
   $Q_0$. It consists of 5 files
   - `Q_zero/equation.jl` - Implements the equation for the ODE as
     well as the recursion for computing the coefficients of the
     Taylor expansion.
   - `Q_zero/Q_float.jl` - Contains the implementation of the
     non-rigorous evaluation of $Q_0$.
   - `Q_zero/Q_taylor.jl` - Contains the code for enclosing $Q_0$
     using the Taylor expansion at zero.
   - `Q_zero/Q_capd.jl` - Contains the Julia code that is responsible
     for calling the CAPD code, which then computes enclosures of
     $Q_0$.
   - `Q_zero/Q.jl` - Defines wrapper methods that call either the
     `Q_float.jl` methods or `Q_capd.jl` methods depending on the type
     of the input. This is the interface that is primarily used by the
     other parts of the code. It defines the functions `Q_zero`,
     `Q_zero_jacobian_kappa` and `Q_zero_jacobian_epsilon`.
3. [`Q_infinity/`](src/Q_infinity) - The Julia code used for
   evaluating $Q_\infty$ is found in this directory. It consists of 9
   files, and is the most delicate part of the code.
   - `Q_infinity/parameters.jl` - Contains code for computing $a$,
     $b$, $c$, $B_W$, $B_{W,\kappa}$ and $B_{W,\epsilon}$.
   - `Q_infinity/functions.jl` - Contains code for evaluating
     $P(\xi)$, $E(\xi)$ and related functions. The functions follow
     the same naming scheme as in the paper, with the derivatives
     denoted by a postfix, for example function `P_dξ_dξ_dκ` is used
     to compute $P_\kappa''$.
   - `Q_infinity/function_bounds.jl` - Responsible for computing the
     bounds that occur in Lemma 7.3 in the paper. The bounds are
     collected into one struct `C::FunctionBounds`, from which the
     bounds can be accessed. For example, `C.E_dξ_dξ_dξ` corresponds
     to $C_{E'''}$ in the paper.
   - `Q_infinity/I_bounds.jl` - The bounds from Sections 7.3.2 and
     7.3.3 in the paper are implemented here. The bounds are collected
     into one struct `CI::IBounds`, from which the bounds can be
     accessed. For example, `CI.I_P_2_4` corresponds to $C_{I_P,2,4}$
     in the paper.
   - `Q_infinity/norm_bounds.jl` and
     `Q_infinity/norm_bounds_constants.jl` - Contains the code for
     computing initial bounds of the norms of $Q_\infty$ and its
     derivatives. The bounds are collected into one struct
     `norms::NormBounds`, from which the bounds can be accessed. For
     example, `norms.Q_dϵ_dξ` is a bound for
     $\|Q_\epsilon'\|_\mathrm{v}$. For $\|Q\|_\mathrm{v}$, the bounds
     are based on the fixed-point approach in Section 7.2. The other
     bounds are based on the results in Section 7.3.3. The constants
     are given in `Q_infinity/norm_bounds_constants.jl` and the bounds
     are computed in `Q_infinity/norm_bounds.jl`.
   - `Q_infinity/I_P.jl` - This code computes enclosures, not only
     bounds, for $I_P$, $I_{P,\gamma}$, $I_{P,\kappa}$, and
     $I_{P,\epsilon}$. It is based on the approach discussed in
     Section 7.3.4.
   - `Q_infinity/Q.jl` - Puts all of the above code together to
     compute enclosures for $Q$ and its Jacobian. It defines the
     functions `Q_infinity`, `Q_infinity_jacobian_kappa` and
     `Q_infinity_jacobian_epsilon`.
   - `Q_infinity/verify_monotonicity.jl` - Implements the check for
     monotonicity on the interval $(\xi_2, \infty)$ that is described
     in Section 5.1. In particular it implements the bounds found in
     Section 7.4.
4. [`branch/`](src/branch) - This directory contains the main logic
   for extending the methods applied to the NLS equation to handle the
   branches for the CGL equation. A lot of the code is tasked with
   handling dynamic bisection and distributing tasks to workers. It is
   structured to allow for efficient computations on an HPC cluster.
   - `branch/branch_points.jl` - Handles pointwise verification along
     the branches. The results in Section 9 are produced with this.
   - `branch/branch_existence.jl` and
     `branch/branch_segment_existence.jl` - Dynamically splits the
     branches into small segments and proves local existence for each
     segment.
   - `branch/branch_continuation.jl` - Takes the output of
     `branch/branch_existence.jl` and checks that the segments can be
     joined together into one continuous curve, bisects segments
     further if required.
   - `branch/branch_critical_points.jl` - Takes the output of
     `branch/branch_continuation.jl` and computes the number of
     critical points for each segment.
   - `branch/proof_witness.jl` - Responsible for taking the output of
     all above parts and joining it into one dataset that contains all
     the data that is required for the proofs.
   - `branch/data_handling.jl` - Contains functions for creating,
     writing and reading dataframes containing the produced data.
5. [`orchestration/`](src/orchestration) - Contains a significant amount
   of gluing code to make calling the functions in `branch/` simpler.
   In particular these functions handle loading of precomputed data as
   well as formatting and storing the output. The are primarily
   intended to be used through the HPC interface, see
   [`HPC/README.md`](HPC/README.md).
6. The files in the root of the [`src/`](src) directory handle a wide
   variety of tasks.
   - `arb.jl`, `helper.jl`, `special-functions.jl` - These files
     implement a variety of convenience functions, such as formatting
     of output and progress logging, as well as some basic
     mathematical functions such as the rising factorial.
   - `verify_and_refine_root.jl` - Implements the functions related to
     the interval Newton method. See the individual functions
     documentation for more details.
   - `U.jl`, `U_expansion.jl` - Implements evaluation and expansion of
     the confluent hypergeometric function $U$. In particular this
     makes use of the bounds from Lemmas 7.1 and 7.2.
   - `CGLParams.jl` - Defines the type `CGLParams` that is used for
     storing the parameters that are held fixed throughout the
     computations, that is $d$, $\omega$, $\sigma$ and $\delta$. It
     also contains code to generate the initial guesses for the NLS
     equation.
   - `refine_approximation.jl` - Contains the code for non-rigorous
     refinement of initial approximations.
   - `G.jl` - Contains the code for enclosing the function $G$.
   - `G_solve.jl` - Wrapper functions for computing roots of $G$ with
     the help of functions from `verify_and_refine_root.jl`.
   - `count_critical_points.jl` - The procedure that is detailed in
     Section 5.1 for counting the number of critical points of $|Q|$
     is implemented here.

### Diagrams

Here are some incomplete, but potentially helpful, diagrams showing
some of the dependencies in the code.

Dependencies of `G_solve.jl`:

``` mermaid
graph TD;
  Q_zero/Q.jl --> G.jl
  Q_infinity/Q.jl --> G.jl
  G.jl --> G_solve.jl
  verify_and_refine_root.jl --> G_solve.jl
```
