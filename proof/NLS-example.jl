### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try
            Base.loaded_modules[Base.PkgId(
                Base.UUID("6e696c72-6542-2067-7265-42206c756150"),
                "AbstractPlutoDingetjes",
            )].Bonds.initial_value
        catch
            b -> missing
        end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ d18c2f76-73b5-11ef-27a6-6d984691fd6c
begin
    using Pkg, Revise
    Pkg.activate("..", io = devnull)
    using Arblib
    using CGL
    using IntervalArithmetic: interval
    using LaTeXStrings
    using Plots
    using PlutoUI

    # Set the bits of precision used for the computations with Arblib
    setprecision(Arb, 128)

    ENV["LD_LIBRARY_PATH"] = "/home/joeldahne/Programs/capd/lib/"

    nothing
end

# ╔═╡ 8d4fdc56-6f0b-4a20-b07b-f5ac07c235b2
md"""
# A detail example for the NLS equation
In this notebook we go through a detailed example with the computations that are required to prove the existence of a self-similar singular solution to the non-linear Schrödinger equation, as well as to determine the number of critical points of the corresponding profile.

The notebook corresponds to Section 5.1 in the paper with which this repository is associated, it follows the same structure, and all numbers in that section of the paper are coming from this notebook. The notebook serves to makes it explicit what code is run to get these results. Since the output of the notebook is used in the paper, part of the notebook is concerned with formatting the output in a way that is suitable for inclusion in LaTeX.
"""

# ╔═╡ c3053bb7-4415-40e5-8774-3e72ae6e527f
TableOfContents()

# ╔═╡ 8ebc0fbe-8309-4eb0-ad70-d09657316253
md"""
The section in the paper goes through the details for the first self-simila singular solution of the 3D cubic non-linear Schrödinger equation. This corresponds to $j = 1$ and $d = 3$.
"""

# ╔═╡ 2c73a770-ec32-4c0b-9c74-31676f85eb32
j, d = 1, 3

# ╔═╡ 259a5c21-1497-4343-a67c-4b8d930d8f44
md"""
Check this box to set the code to save the figures.
- Save figures $(@bind save CheckBox(default = false))
"""

# ╔═╡ c88671ba-f5a2-4293-a3b9-9cd8e5ad7cf0
guidefontsize, tickfontsize = if save
    pgfplotsx()
    24, 24
else
    gr()
    11, 11
end

# ╔═╡ 08229fb6-98af-4581-af0a-2587e8861be2
md"""
## Mathematical setup
Recall that the proof of existence is based on proving the existence of a zero of the function $G(\mu, \gamma, \kappa)$ defined in Section 2 of the paper (we fix $\epsilon = 0$ since we are considering the NLS case). To prove the existence of a root we make use of the so called **interval Newton method**. By splitting $\gamma$ into real and imaginary parts, we can treat $G$ as a function from $\mathbb{R}^4$ to $\mathbb{R}^4$. To apply the interval Newton method we need to find a set $X \subseteq \mathbb{R}^4$, such that

$$\operatorname{mid}(X) - J_G^{-1}(X) G(\operatorname{mid}(X)) \subsetneq X,$$

where $\operatorname{mid}(X)$ denotes the midpoint of $X$ and $J_G^{-1}$ denotes the inverse of the jacobian of $G$. If we find $X$ satisfying this condition, then the function $G$ has a unique zero in $X$, and this zero is contained in the set given by the left hand side of the above expression.

The proof only requires us to verify this condition, but a large part of the code is involved in first finding this set $X$.
"""

# ╔═╡ 4583b371-c009-42cd-8956-6be3b1de3b7b
md"""
## Getting an approximate solution
The first step in finding the set $X$ is to find a good numerical approximation, this is handled by the function `sverak_params`. This function starts with a rough, hard coded, approximation of the zero, which is internally refined using a few Newton iterations, through the `refine_approximation_fix_epsilon` function.
"""

# ╔═╡ d9078fa5-0e06-49a2-b947-54a8370e0952
μ₀, γ₀, κ₀, ϵ, ξ₁, λ = CGL.sverak_params(Arb, j, d)

# ╔═╡ d102bce4-0e13-4644-89f9-f0044bdff092
md"""
The variables `µ₀`, `γ₀` and `κ₀` holds our approximate zero
"""

# ╔═╡ 8487b9c4-bd77-4272-9d6b-3ab043d9e305
μ₀, γ₀, κ₀

# ╔═╡ 2fc0db5b-0d71-4f93-8e7b-0312103eaeeb
md"""
Note that even though these approximations are printed as intervals above, they are in fact exact (floating point) numbers. To see the exact numbers we need to manually specify a larger number of digits to print.
"""

# ╔═╡ 42d3826b-c8f6-47f1-9e4b-19d6c2e55887
string(μ₀, digits = 64)

# ╔═╡ 9abf6a82-301e-4065-9913-54e6da57596f
string(γ₀, digits = 64)

# ╔═╡ f6a4ac88-a959-4e5e-9d3d-9f9b81033c1f
string(κ₀, digits = 64)

# ╔═╡ aec77170-7845-4390-91f8-c8cba67c0115
md"""
The value for $ϵ$ is zero, since we are in the NLS case.
"""

# ╔═╡ c9ae053a-ca79-4188-9c01-f9d8953e0cb9
iszero(ϵ) # We are in the NLS case, so ϵ = 0

# ╔═╡ ae76061d-4435-4376-9456-0c79dfc6cd6e
md"""
The value for `ξ₁` determines at which point the solution from zero and the solution at infinity are matched.
"""

# ╔═╡ 037453b2-caa8-47a2-aa48-83f08a7c948f
md"""
The variable `λ` holds the parameters $d$, $\omega$, $\sigma$ and $\delta$ that, together with $\kappa$ and $\epsilon$, define our ODE.
"""

# ╔═╡ 0d425114-1f17-4297-ab51-1aef6df765b8
λ

# ╔═╡ 9795fa90-9821-485a-9583-04634d0f7f0f
md"""
We can evaluate the function $G$ on our approximation, giving us
"""

# ╔═╡ 467a5c66-9a07-4fbc-b738-39dce6d2be1a
G_approximation = CGL.G(μ₀, real(γ₀), imag(γ₀), κ₀, ϵ, ξ₁, λ)

# ╔═╡ 5cb7cb1f-66e6-466a-bfbc-2e1a1d47e0e0
md"""
As seen the approximation is good enough so that the error bounds, coming from the evaluation of $Q_0$ and $Q_\infty$, make the resulting enclosure straddle zero. Note that by default only significant digits of the results are printed, we can get a more precise view of the result by forcing it to print `more` digits.
"""

# ╔═╡ 17eee8c0-15c3-4986-8b5c-a7422ab16ece
string.(G_approximation, more = true)

# ╔═╡ 09578453-f174-4957-a0a6-238c88d56a57
md"""
Similarly, we can also compute the jacobian of $G$ at our approximation.
"""

# ╔═╡ 7fbc5c46-2164-445a-99e0-f2d3750a11e6
J_G_approximation = CGL.G_jacobian_kappa(μ₀, real(γ₀), imag(γ₀), κ₀, ϵ, ξ₁, λ)

# ╔═╡ f2e891bd-7669-43e6-9ab1-286dee91d221
md"""
### Formatting output for paper
To include the results in the paper we want to format it as latex code.
"""

# ╔═╡ 06b734d2-7362-4f78-b737-68fcc94e603c
µ₀_latex = string(µ₀, no_radius = true, digits = 16)

# ╔═╡ 8b61869b-5c95-4392-a07d-5de8dda85bd7
γ₀_latex = replace(string(γ₀, no_radius = true, digits = 16), "+ -" => "-", "im" => "i")

# ╔═╡ d9137251-fa69-4942-bd70-983f75adad30
κ₀_latex = string(κ₀, no_radius = true, digits = 16)

# ╔═╡ b3f55a5b-9a79-4446-8f56-adf774634c86
approximaton_latex = """
\\begin{align*}
  \\mu_0 &= $(μ₀_latex),\\\\
  \\gamma_0 &= $(γ₀_latex),\\\\
  \\kappa_0 &= $(κ₀_latex).
\\end{align*}
"""

# ╔═╡ 5a756d4f-297e-4f4c-bfc8-eda7063a51b9
md"""
The initial approximation is, rounded to 16 digits, given by
"""

# ╔═╡ 85b544a0-3480-4bb1-9565-ef889b686ed6
latexstring(approximaton_latex)

# ╔═╡ ef43800a-59dd-4724-af00-df370792e58b
print(approximaton_latex)

# ╔═╡ a8c8b1e0-cf1b-462c-8ae1-6b898c35f8ab
G_approximation_latex = let
    str_1 =
        "\\big(" * join(CGL.format_interval_precise.(G_approximation[1:2]), ", ") * ",\\\\"
    str_2 = join(CGL.format_interval_precise.(G_approximation[3:4]), ", ") * "\\big)"

    """
    \\begin{multline*}
      G(\\mu_0, \\gamma_0, \\kappa_0) \\in 
      $str_1
      $str_2.
    \\end{multline*}
    """
end

# ╔═╡ 54af636d-7b4b-4cc8-9847-1186fd66df86
md"""
We have
"""

# ╔═╡ 36ebd7ae-ec94-4b43-83da-bb0fd836ce40
latexstring(G_approximation_latex)

# ╔═╡ fcc09ed6-29bd-45af-b731-77e9e6ccf8b8
print(G_approximation_latex)

# ╔═╡ 3b0b26c9-f0a0-4564-b635-636f40d8279e
J_G_approximation_latex = let
    J_G_approximation_rows_latex = map(
        row -> join(row, " & "),
        eachrow(CGL.format_interval_precise.(J_G_approximation, min_digits = 4)),
    )

    """
    \\begin{equation*}
      J_G(\\mu_0, \\gamma_0, \\kappa_0) \\in
      \\begin{pmatrix}
        $(J_G_approximation_rows_latex[1]) \\\\
        $(J_G_approximation_rows_latex[2]) \\\\
        $(J_G_approximation_rows_latex[3]) \\\\
        $(J_G_approximation_rows_latex[4]) \\\\
      \\end{pmatrix}.
    \\end{equation*}
    """
end

# ╔═╡ 4c14198f-531f-4bf6-82d4-94e315034d15
md"""
Similarly, we have
"""

# ╔═╡ a3568944-2444-456c-b262-de32dedf7c06
latexstring(J_G_approximation_latex)

# ╔═╡ 939c60ab-b10a-41d1-8f67-823840f5c728
print(J_G_approximation_latex)

# ╔═╡ 520c152b-3196-4d7c-9950-2e25c6898354
md"""
## Enclosing the root
With our approximation of the root, we next want to compute a rigorous enclosure of the root. This is handled by the function `G_solve_fix_epsilon`, it takes as argument our numerical approximation, our parameters and some flags telling it exactly what to compute.

What the function does is that it takes numerical approximation and then dynamically tries to find a set $X \subseteq \mathbb{R}^4$ around this approximation, such that the condition of the interval Newton method is satisfied. One has to be careful when choosing $X$, both too small and too large sets will make the required inclusion not be satisfied. Exactly how it determines what set to use is outside the scope of this notebook, see Section 10 in the paper. However, from the values returned by the function we can a posteriori verify that the required condition is satisfied (this is automatically done by the function, we only do the a posteriori check here to show the procedure).
"""

# ╔═╡ 82d39edb-ecbf-496c-9d09-476d369a7b24
root, root_uniqueness = CGL.G_solve_fix_epsilon(
    μ₀,
    real(γ₀),
    imag(γ₀),
    κ₀,
    ϵ,
    ξ₁,
    λ,
    return_uniqueness = Val(true),
    try_expand_uniqueness = false,
    verbose = true,
)

# ╔═╡ 62a35f0e-dab9-468f-b52a-42b2c2de8b34
all(isfinite, root) # Verify that it succeeded enclosing the root

# ╔═╡ 9b5f5289-9256-4d99-94b7-53aad0df45dd
md"""
Here the variable `root_uniqueness` corresponds to the set $X$ and `root` corresponds to the left hand side of the required interval Newton inclusion. In terms of the variables $\mu$, $\gamma$ and $\kappa$ we have the following enclosures for the root.
"""

# ╔═╡ 8f7e5899-94da-4f5b-b2c3-b6395be86711
μ, γ, κ = root[1], Acb(root[2], root[3]), root[4]

# ╔═╡ 2194a21f-09bb-4cff-883f-3b26f289ce45
md"""
And the root is unique in the following sets.
"""

# ╔═╡ 9fde3fa4-0c85-4ece-a5b0-e09f44ccb370
μ_uniqueness, γ_uniqueness, κ_uniqueness =
    root_uniqueness[1], Acb(root_uniqueness[2], root_uniqueness[3]), root_uniqueness[4]

# ╔═╡ 2aeddc6f-05b8-4545-931c-08c55c5b767e
md"""
Our goal is now to verify the condition for the interval Newton method. We then have the set $X$ given by
"""

# ╔═╡ 71c8a3aa-374b-4eb1-bfd7-41d989d3ec65
X = root_uniqueness

# ╔═╡ 7190d92b-7f62-4805-bae1-d23e869c736b
md"""
For $\operatorname{mid}(X)$ we get
"""

# ╔═╡ b3a3e9d0-3ea0-4523-8141-8ee0cb68741a
X_mid = midpoint.(Arb, root_uniqueness)

# ╔═╡ 4a12122c-c7dd-4e14-aff7-621c3badc910
md"""
We can compute $G(\operatorname{mid}(X))$ as
"""

# ╔═╡ b35ec786-9c68-4d09-88c8-a3a316bbba31
G_X_mid = CGL.G(X_mid..., ϵ, ξ₁, λ)

# ╔═╡ 2867d01e-4839-4200-9fa9-938856d6e521
md"""
And $J_G(X)$ as
"""

# ╔═╡ cf7a3fab-ea03-4eab-a3d3-291611d74672
J_G_X = CGL.G_jacobian_kappa(X..., ϵ, ξ₁, λ)

# ╔═╡ eddb69c4-13a0-4dde-8eb1-1bc5249cde09
md"""
To rigorously compute the inverse of $J_G(X)$ we need to use some of the lower level methods of Arblib.jl, this requires us to first convert the matrix to an `ArbMatrix`, apply the `Arblib.inv!` method and then convert back.
"""

# ╔═╡ 2a603618-19ef-4645-8977-2ec4b2d9b630
J_G_X_inv = oftype(J_G_X, inv(ArbMatrix(J_G_X)))

# ╔═╡ dbfa101a-dec9-4b27-a2b8-1c70cd2e128c
md"""
However, to for computing $J_G^{-1}(X) G(\operatorname{mid}(X))$ it is, as always with numerical methods, better to not directly compute the inverse of $J_G$ and left multiply with it, but rather to solve the corresponding linear system directly. In Julia this is generally done using the backslash operator `\`, however, to get rigorous enclosures we again need to reach for the lower level methods of Arblib.jl.
"""

# ╔═╡ 41ef9478-0f84-4eba-bc54-9a383f15565d
J_G_X_ldiv_G_X_mid = let
    J_G_X_ldiv_G_X_mid = ArbMatrix(4, 1)
    success =
        !iszero(Arblib.solve!(J_G_X_ldiv_G_X_mid, ArbMatrix(J_G_X), ArbMatrix(G_X_mid)))
    @assert success # Verify that the operation succeeded
    oftype(G_X_mid, J_G_X_ldiv_G_X_mid)
end

# ╔═╡ 8b43f3a6-7315-433d-97ea-2e21ccadd899
md"""
With this we get that the left hand side of the required inclusion is given by
"""

# ╔═╡ 2b1fc470-59a2-4da5-9f9d-d72ba7b96216
X_newton = X_mid - J_G_X_ldiv_G_X_mid

# ╔═╡ 1181d465-27c9-41fa-89c6-dc70bca6c286
md"""
And we can verify that this indeed is strictly included in $X$.
"""

# ╔═╡ d31e6dc5-79dc-4101-958c-0c44919ac865
all(Arblib.contains_interior.(X, X_newton))

# ╔═╡ d1e6fabf-6f4f-4608-884d-7fd025d62345
md"""
Note that the right hand side is exactly the `root` result returned by `G_solve_fix_epsilon`, they are computed in the same way.
"""

# ╔═╡ 95d3cba6-1cdb-47f8-bbdf-db0620ac5ea6
isequal(root, X_newton)

# ╔═╡ f7bce6c5-8b02-467a-ae18-7115540301d7
md"""
### Formatting output for paper
To include the results in the paper we want to format it as latex code.
"""

# ╔═╡ e5067db1-e7db-4703-a162-ff5eafbc153b
# TODO: Should we add ... notation for upper and lower bounds?

# ╔═╡ 06473469-0f0c-483f-af67-e1b4c17652f9
# TODO: These should be rounded inwards

# ╔═╡ 83e16fe2-6d33-4dcf-a8e1-abffc043d869
μ_uniqueness_latex = CGL.format_interval_precise(μ_uniqueness, min_digits = 4)

# ╔═╡ 11d9ac67-075b-47f5-bdc0-efe990dc40d2
γ_real_uniqueness_latex = CGL.format_interval_precise(real(γ_uniqueness), min_digits = 4)

# ╔═╡ b2f87105-4083-469c-ba8d-556d980804ad
γ_imag_uniqueness_latex = CGL.format_interval_precise(imag(γ_uniqueness), min_digits = 4)

# ╔═╡ 38dbeb52-743b-4020-9aba-36649bbd9c48
κ_uniqueness_latex = CGL.format_interval_precise(κ_uniqueness, min_digits = 4)

# ╔═╡ 1a0ef6c0-293a-4f73-a013-d60adbefbfa2
root_uniqueness_latex = """
\\begin{equation*}
  X = $(μ_uniqueness_latex) \\times 
  $(γ_real_uniqueness_latex) \\times
  $(γ_imag_uniqueness_latex) \\times
  $(κ_uniqueness_latex).
\\end{equation*}
"""

# ╔═╡ 1f13b882-5e96-4fee-84aa-c6480e6f30c8
md"""
Moreover, this root is unique for $(\mu, \operatorname{Re}(\gamma), \operatorname{Im}(\gamma), \kappa)$ in the set
"""

# ╔═╡ ea662336-745c-438c-bdc2-acbfb151ed34
latexstring(root_uniqueness_latex)

# ╔═╡ 861405c9-92ea-46ee-b61a-e937ebaf7a45
print(root_uniqueness_latex)

# ╔═╡ ec85ad10-2ed8-4695-b2cd-e550e4994bbc
# TODO: Add support for intervals containing zero to format_interval_precise

# ╔═╡ 2de7e6fd-a952-43aa-833e-f66e3749b462
G_X_mid_latex = let
    str = "\\left(" * join(CGL.format_interval_precise.(G_X_mid), ", ") * "\\right)"

    """
    \\begin{equation*}
      G(\\operatorname{mid}(X)) \\in $str,
    \\end{equation*}
    """
end

# ╔═╡ bc60cae3-001a-453c-be3b-e5eaaea38efe
J_G_X_latex = let
    J_G_X_rows_latex = map(
        row -> join(row, " & "),
        eachrow(CGL.format_interval_precise.(J_G_X, min_digits = 4)),
    )

    """
    \\begin{equation*}
      J_G(X) \\subseteq
      \\begin{pmatrix}
        $(J_G_X_rows_latex[1]) \\\\
        $(J_G_X_rows_latex[2]) \\\\
        $(J_G_X_rows_latex[3]) \\\\
        $(J_G_X_rows_latex[4]) \\\\
      \\end{pmatrix}.
    \\end{equation*}
    """
end

# ╔═╡ 6dd19686-d49e-4acb-b301-59c74969b307
md"""
We have
"""

# ╔═╡ b1ff87ee-2e7b-450c-aaea-dd76b4e11718
latexstring(G_X_mid_latex)

# ╔═╡ c1bfa2ec-a316-4cd6-82b5-cf59b1ae2a6a
md"""
and
"""

# ╔═╡ 0f012491-884f-4930-b95c-0a7c255113bf
latexstring(J_G_X_latex)

# ╔═╡ 919b631a-c4a7-4504-b346-1a0040a490ad
print(G_X_mid_latex)

# ╔═╡ 478381cc-0f8a-4f4c-9919-fac3744a5a2c
print(J_G_X_latex)

# ╔═╡ 2ddd1d45-38f9-4ecd-baef-885e8988b618
J_G_X_inv_latex = let
    J_G_X_inv_rows_latex = map(
        row -> join(row, " & "),
        eachrow(CGL.format_interval_precise.(J_G_X_inv, min_digits = 4)),
    )

    """
    \\begin{equation*}
      J_G^{-1}(X) \\subseteq
      \\begin{pmatrix}
        $(J_G_X_inv_rows_latex[1]) \\\\
        $(J_G_X_inv_rows_latex[2]) \\\\
        $(J_G_X_inv_rows_latex[3]) \\\\
        $(J_G_X_inv_rows_latex[4]) \\\\
      \\end{pmatrix}.
    \\end{equation*}
    """
end

# ╔═╡ 188adc35-bbb6-41a1-b189-9a2955281abc
md"""
We have
"""

# ╔═╡ 66c23812-1859-405b-816d-08f36f958080
latexstring(J_G_X_inv_latex)

# ╔═╡ af0f44fd-17df-4f0c-b7ca-84644fe08276
print(J_G_X_inv_latex)

# ╔═╡ 7900e0e5-9817-403b-bbe2-9866138efb1e
J_G_X_ldiv_G_X_mid_latex = let
    str =
        "\\left(" *
        join(CGL.format_interval_precise.(J_G_X_ldiv_G_X_mid), ", ") *
        "\\right)"

    """
    \\begin{equation*}
      J_G^{-1}(X) G(\\operatorname{mid}(X))\\subseteq
      $str.
    \\end{equation*}
    """
end

# ╔═╡ 9897916c-ea93-46d4-9ffd-a74c160f1301
latexstring(J_G_X_ldiv_G_X_mid_latex)

# ╔═╡ df6e6400-b9c4-4494-adc2-98028ee86a6d
print(J_G_X_ldiv_G_X_mid_latex)

# ╔═╡ 5e796977-0bff-49d7-aeec-322714d9471d
X_newton_latex = let
    str =
        "\\left(" *
        join(CGL.format_interval_precise.(X_newton, min_digits = 4), ", ") *
        "\\right)"

    """
    \\begin{equation*}
      \\operatorname{mid}(X) - J_G^{-1}(X) G(\\operatorname{mid}(X))\\subseteq
      $str.
    \\end{equation*}
    """
end

# ╔═╡ 3bf6078f-9f62-4d24-bdcc-769842c78f67
md"""
We have
"""

# ╔═╡ 48f98796-75c4-4851-817d-558caba65bbf
latexstring(X_newton_latex)

# ╔═╡ d376330e-1793-4652-815d-0a091feb9024
print(X_newton_latex)

# ╔═╡ b2b46e6c-f3a7-4713-b911-27a2462eadb5
μ_latex = CGL.format_interval_precise(μ, min_digits = 4)

# ╔═╡ fe7e0856-a161-4431-996e-6ffb7438a4b4
γ_latex = CGL.format_interval_precise(γ, min_digits = 4)

# ╔═╡ 314e4e16-6e8f-46a3-8d99-2df74d73560a
κ_latex = CGL.format_interval_precise(κ, min_digits = 4)

# ╔═╡ 588ffcb8-2865-4c46-af76-dae8e8dc1e2a
root_latex = """
\\begin{equation*}
  \\mu \\in $(μ_latex),\\quad
  \\gamma \\in $(γ_latex),\\quad
  \\kappa \\in $(κ_latex).
\\end{equation*}
"""

# ╔═╡ 2e881d99-aa92-4ef6-afac-65990aff2487
md"""
There exists a root $(\mu, \gamma, \kappa)$ of $G$, with
"""

# ╔═╡ 27375192-be08-48c0-a872-dbf630e03ca4
latexstring(root_latex)

# ╔═╡ 4b9a7d4c-8828-41e3-931e-b5ac76fbf2cb
print(root_latex)

# ╔═╡ 88a12187-3dcd-4dbf-81e5-d80232d63dbf
md"""
## Plotting the profile
"""

# ╔═╡ dbea5cae-4150-4155-910e-e52b6952019d
md"""
With the computed encloses for $\mu$, $\gamma$ and $\kappa$ we can, with the help of the rigorous numerical integrator implemented by CAPD, compute enclosures of $Q$ on the whole interval $[0, \xi_1]$. Below we plot first the real and imaginary parts of $Q$, and then the absolute value of $Q$.
"""

# ╔═╡ f09c6a45-7116-44ac-ac7f-696b084935a1
ξ₁s_plot, Q₀s_plot = CGL.Q_zero_capd_curve(μ, κ, ϵ, ξ₁, λ)[1:2]

# ╔═╡ db3ad4c3-7c9c-4dad-a097-eb63f22ade38
let pl = plot(xlabel = L"\xi"; guidefontsize, tickfontsize)
    plot!(
        pl,
        vcat.(interval.(ξ₁s_plot), interval.(getindex.(Q₀s_plot, 1))),
        label = L"Re(Q)", #L"\operatorname{Re}(Q)", # operatorname doesn't work
    )
    plot!(
        pl,
        vcat.(interval.(ξ₁s_plot), interval.(getindex.(Q₀s_plot, 2))),
        label = L"Im(Q)", #L"\operatorname{Im}(Q)", # operatorname doesn't work
    )

    save && savefig(pl, "figures/NLS-example-real-imaginary.pdf")

    pl
end

# ╔═╡ 1f3a56af-070e-4204-a54f-d7aa7e9c47fc
let pl = plot(xlabel = L"\xi"; guidefontsize, tickfontsize)
    abs_Q₀s_plot = abs.(Acb.(getindex.(Q₀s_plot, 1), getindex.(Q₀s_plot, 2)))
    plot!(pl, vcat.(interval.(ξ₁s_plot), interval.(abs_Q₀s_plot)), label = L"|Q|")

    save && savefig(pl, "figures/NLS-example-abs.pdf")

    pl
end

# ╔═╡ 2e0ff47e-ea61-424b-bac5-b8022e1de9ce
md"""
## Counting the number of critical points
The above plot of $|Q|$ indicates that it has
"""

# ╔═╡ 4f299a4d-a7a8-4353-a0b1-c771e6d77dfc
L"j - 1 = %$(j - 1)"

# ╔═╡ 94a05218-0ae9-44d1-a5f2-439c5b1cbd4c
md"""
critical points on the interval $(0, \xi_1]$. We will prove that this is the case, and it holds true also on the full interval $(0, \infty)$. To prove this we want to count the number of zeros of $|Q|'$, for which we split the interval $(0, \infty)$ into three parts
- $(0, \xi_0).$
- $[\xi_0, \xi_2],$
- $(\xi_2, \infty).$
On the first and last part we prove that $|Q|'$ is non-zero, and on the middle part we prove that it has exactly $j - 1$ zeros.

In practice we work with $\frac{d}{d\xi}|Q|^2$ instead of $|Q|$' directly, since that gives slightly simpler formulas.
"""

# ╔═╡ 833eb4e0-7c7c-44f6-94e7-02d4859b7c4c
md"""
Let us start with the interval $(\xi_2, \infty)$. From Lemma 7.2 in the paper, we have that for $\xi > \xi_1$,

$$\frac{d}{d\xi}|Q|^{2}
    = p_{X}(\lambda, \xi)\xi^{-\frac{2}{\sigma} - 1}
    + R_{X}(\lambda, \xi)\xi^{(2\sigma + 1)\mathrm{v} - \frac{2}{\sigma} - 3}.$$

with

$$|p_{X}(\lambda, \xi)| \geq C_{p_{X}}(\lambda)
  \text{ and }
  |R_{X}(\lambda, \xi)| \leq C_{R_{X}}(\lambda).$$

Taking $\xi_{2} \geq \xi_{1}$ such that

$$\xi_{2} \geq
  \left(\frac{C_{p_{X}}}{C_{R_{X}}}\right)^{\frac{1}{(2\sigma + 1)\mathrm{v} - 2}},$$

is then enough to assert that $\frac{d}{d\xi}|Q|^{2}$ is non-zero on
the interval $(\xi_{2}, \infty)$.
"""

# ╔═╡ 001294b1-866c-43d7-ac9f-e6c602363d25
md"""
The function responsible for checking this condition is `verify_monotonicity_infinity`.
"""

# ╔═╡ a74ec90f-4673-4720-a34a-696d19015321
ξ₂, C_p_X, C_R_X, ξ₂_lower_bound =
    CGL.verify_monotonicity_infinity(γ, κ, ϵ, ξ₁, λ, return_coefficients = Val(true))

# ╔═╡ 3f730b42-2456-4b49-a50e-27b3acbd4573
md"""
It returns
- `ξ₂`: Value for $\xi_2$ such that the above inequality is guaranteed to be satisfied for all $\xi > \xi_2$.
- `C_p_X` and `C_R_X`: Enclosures of $C_{p_X}$ and $C_{R_X}$.
- `ξ₂_lower_bound`: Enclosure of $\left(\frac{C_{p_{X}}}{C_{R_{X}}}\right)^{\frac{1}{(2\sigma + 1)\mathrm{v} - 2}}$. 
Taking the maximum of `ξ₂_lower_bound` and $\xi_1$ gives us $\xi_2$. In the case of NLS it will generally return $\xi_2 = \xi_1$, though along the branches for CGL it is sometimes necessary to use a larger value for $\xi_2$.
"""

# ╔═╡ 745c36db-56c6-46af-8b21-ea8ea8aa7a20
isequal(ξ₁, ξ₂)

# ╔═╡ d8c514e4-cf96-4431-9041-1ca4a497a399
md"""
With the monotonicity checked on $(\xi_2, \infty)$, what remains is checking $(0, \xi_0)$ and $[\xi_0, \xi_2]$. This is done by using the rigorous numerical integrator to enclose the function on the interval $[0, \xi_2]$, similar to how the above plots of $Q$ were produced.
"""

# ╔═╡ e0217811-4851-47a6-ac5a-e4c163a919ef
ξ₁s, _, _, abs2_Q_derivatives, abs2_Q_derivative2s = CGL.Q_zero_capd_curve(μ, κ, ϵ, ξ₂, λ)

# ╔═╡ 94dbdfe3-a5ce-47c0-bf6b-c582603c10cc
md"""
TODO: Write about return values
"""

# ╔═╡ 68508799-c37b-4c4f-932f-0750d915f892
md"""
To determine $\xi_0$ we look for the first enclosure of $\frac{d}{d\xi}|Q|^{2}$ which is non-zero. The $\xi_0$ value is then given by the upper bound of $\xi$ for the preceeding enclosure.
"""

# ╔═╡ a356d6d9-a004-431c-8c85-d5ba6f79ba66
i = findfirst(!Arblib.contains_zero, abs2_Q_derivatives)

# ╔═╡ 2b043bd4-aacf-43be-8141-578bb0928f86
ξ₀ = ubound(Arb, ξ₁s[i-1])

# ╔═╡ 4f9e593b-d884-4601-88c2-0d806c4e8998
md"""
We can now verify the condition on the interval $[ξ₀, ξ₂]$. In the case $j = 1$ it is enough to verify that the enclosure of $\frac{d}{d\xi}|Q|^{2}$ is non-zero for all indices above `i`. For $j > 1$ we need to count the number of zeros, this is straight forward but slightly tedious to implement. The implementation of that is given in `count_critical_points`, which is discussed further down. In this notebook we only treat the case $j = 1$.
"""

# ╔═╡ 7b2938c0-6831-49ae-bb68-82ce688147b6
if j == 1
    all(!Arblib.contains_zero, abs2_Q_derivatives[i:end])
end

# ╔═╡ 89a29d15-8711-4168-9d9a-cf8adc86b83a
md"""
What remains is to handle the interval $(0, \xi_0)$. Due to the boundary conditions at $\xi = 0$ we have that $\frac{d}{d\xi}|Q|^{2}$ is zero at $\xi = 0$, to verify that it is non-zero on $(0, \xi_0)$ it is therefore enough to verify that $\frac{d^2}{d\xi^2}|Q|^{2}$ is non-zero on $[0, \xi_0]$. The `Q_zero_capd_curve` function used for enclosing the curve on $[0, \xi_2]$ also returns an enclosure of $\frac{d^2}{d\xi^2}|Q|^{2}$, however, it is generally (at least in the case $d = 3$) not good enough to verify that it is non-zero.
"""

# ╔═╡ d47f8a46-e517-4251-b30b-c8c91a51e483
Arblib.contains_zero(reduce(Arblib.union, abs2_Q_derivative2s[1:i-1]))

# ╔═╡ 97614c71-26c6-4dc6-974d-240fd281743c
md"""
Instead we enclose the $\frac{d^2}{d\xi^2}|Q|^{2}$ on the interval $[0, \xi_0]$ by Taylor expanding at zero.
"""

# ╔═╡ 392baa05-6032-4b6f-a1d3-15548aa59ecb
abs2_Q_derivative2_ξ₀ = if d == 3
    # Compute enclosures of a = real(Q) and b = imag(Q) as well as their
    # first and second order derivatives (da, db, d2a and d2b).
    (a_ξ₀, b_ξ₀, da_ξ₀, db_ξ₀), (d2a_ξ₀, d2b_ξ₀) =
        CGL.Q_zero_taylor(μ, κ, ϵ, ξ₀, λ, enclose_curve = Val{true}())

    2(d2a_ξ₀ * a_ξ₀ + da_ξ₀^2 + d2b_ξ₀ * b_ξ₀ + db_ξ₀^2)
else
    reduce(Arblib.union, abs2_Q_derivative2s[1:i-1])
end

# ╔═╡ 41e0c536-08b2-4f18-86f6-a46222a097bc
md"""
We can now verify that the enclosure indeed does not contain zero.
"""

# ╔═╡ af1e9e82-14e9-479e-81d3-1f7a2a4b8a62
Arblib.contains_zero(abs2_Q_derivative2_ξ₀)

# ╔═╡ 57b6013f-6fe1-4798-a65d-80fc246f0eee
md"""
This whole procedure for counting the number of critical points is implemented in the function `count_critical_points`.
"""

# ╔═╡ 83bb7eab-921a-44f6-9f81-736624982e44
success, critical_points, verified =
    CGL.count_critical_points(μ, γ, κ, ϵ, ξ₁, λ, verbose = true)

# ╔═╡ 0ee00ce0-c395-4de9-aaaa-021991bf4885
md"""
The `success` variables indicates if the proof succeeded or not.
"""

# ╔═╡ 352e4f92-58e1-4d5b-a44b-52051aaf5589
success # Check that the computation succeeded

# ╔═╡ 8a1909f8-a48f-4479-b0cf-3a0398020864
md"""
The `critical_points` variable contains a list of enclosures of all the $\xi_1$ values for the critical points, the length of the list gives the number of critical points. If there are no critical points the list would then be empty. The `verified` variable is used to indicate which of the critical points are fully verified, if `success` is true then it should always be a list of just true values.
"""

# ╔═╡ b0d4bb71-9c87-4002-bee9-e5db160a29d6
num_critical_points = length(critical_points)

# ╔═╡ fd7e44f4-554b-4ee1-a3e9-4a965a84cb41
md"""
### Formatting output for paper
To include the results in the paper we want to format it as latex code.
"""

# ╔═╡ 624bdfbd-9547-418c-8678-86d19f5c32b4
C_p_X_latex = CGL.format_interval_precise(C_p_X)

# ╔═╡ 716f515d-fcb6-41a5-b3ab-f0aedd489836
C_R_X_latex = CGL.format_interval_precise(C_R_X)

# ╔═╡ 288e71ec-adb2-4b21-8d29-dd6ca0fb3e20
monotonicity_infinity_latex = """
\\begin{equation*}
  C_{p_X} \\in $(C_p_X_latex) \\text{ and }
  C_{R_X} \\in $(C_R_X_latex),
\\end{equation*}
"""

# ╔═╡ 91b6cc42-6e5d-4dec-ab34-bec5fe73be48
ξ₂_lower_bound_latex = CGL.format_interval_precise(ξ₂_lower_bound)

# ╔═╡ ee456f31-40c3-40d4-930d-8f6c9930cfd7
ξ₂_lower_bound_equation_latex = """
\\begin{equation*}
  \\left(\\frac{C_{p_{X}}}{C_{R_{X}}}\\right)^{\\frac{1}{(2\\sigma + 1)\\mathrm{v} - 2}} \\in $(ξ₂_lower_bound_latex).
\\end{equation*}
"""

# ╔═╡ 4bffedfe-f887-4540-b0c2-208509682b03
md"""
We have
"""

# ╔═╡ 5ff56d71-170a-494a-a890-2a2f3e5f6755
latexstring(monotonicity_infinity_latex)

# ╔═╡ 8b5c0fbb-38da-450e-b521-e0307af5722b
md"""
giving us
"""

# ╔═╡ 1b45e60c-7cb7-4109-85da-57eef2b27a0b
latexstring(ξ₂_lower_bound_equation_latex)

# ╔═╡ bac5b3d7-a5e2-4323-90a9-85ea2bdfd416
print(monotonicity_infinity_latex)

# ╔═╡ 139ba5bc-f06d-485c-83c8-be141adbc80f
print(ξ₂_lower_bound_equation_latex)

# ╔═╡ b2c1830d-f694-4df6-a416-826e8d1d279a
let pl = plot(xlabel = L"\xi"; guidefontsize, tickfontsize)
    plot!(
        pl,
        vcat.(interval.(ξ₁s), interval.(abs2_Q_derivatives)),
        label = L"\frac{d}{d\xi}|Q|^2",
    )

    save && savefig(pl, "figures/NLS-example-derivative.pdf")

    pl
end

# ╔═╡ a214d93b-caee-42f9-b2af-088f9ad9d253
let pl = plot(xlabel = L"\xi", xlims = (NaN, 2); guidefontsize, tickfontsize)
    plot!(
        pl,
        vcat.(interval.(ξ₁s), interval.(abs2_Q_derivatives)),
        label = L"\frac{d}{d\xi}|Q|^2",
    )

    save && savefig(pl, "figures/NLS-example-derivative-zoom.pdf")

    pl
end

# ╔═╡ f7b2cfa4-a59d-43d2-83aa-44abcd72bdb5
ξ₀_latex = CGL.format_interval_precise(ξ₀)

# ╔═╡ e36acd5f-3a6f-4415-a805-9f3cfa615333
ξ₀_rounded_latex = string(ξ₀, digits = 10, no_radius = true)

# ╔═╡ 0a62baa4-04bc-4df8-9232-1fbc76b75e21
ξ₀_equation_latex = """
\\begin{equation*}
  \\xi_0 = $(ξ₀_latex)
\\end{equation*}
"""

# ╔═╡ 3049ee1c-cc91-4eb9-8eee-bb1df4e6713d
ξ₀_rounded_equation_latex = """
\\begin{equation*}
  \\xi_0 = $(ξ₀_rounded_latex)
\\end{equation*}
"""

# ╔═╡ bc1745d7-2026-42af-b18a-bf9669db543c
md"""
We have
"""

# ╔═╡ 3a951310-9437-4c3c-ad72-cfa9620cbb20
latexstring(ξ₀_equation_latex)

# ╔═╡ 4c7227d2-3322-4dad-98d0-80b7b118ae5d
md"""
which rounded to 10 digits is
"""

# ╔═╡ 540ab5f3-226b-45ab-956a-2b371ec975af
latexstring(ξ₀_rounded_equation_latex)

# ╔═╡ d9a797ed-54f5-4a10-8ab9-8c7ad6aff38a
print(ξ₀_equation_latex)

# ╔═╡ 3bbffd81-268a-408d-a150-c9a65f074e5a
print(ξ₀_rounded_equation_latex)

# ╔═╡ 2cf059d7-a966-4dd8-b5a2-bb3b01752e08
abs2_Q_derivative2_ξ₀_latex = CGL.format_interval_precise(abs2_Q_derivative2_ξ₀)

# ╔═╡ 89c62802-6242-4f7a-ab6e-cba6851d5ff5
abs2_Q_derivative2_ξ₀_equation_latex = """
\\begin{equation*}
  $(abs2_Q_derivative2_ξ₀_latex).
\\end{equation*}
"""

# ╔═╡ 7c898b59-e6af-4197-ba06-9ed474373eb9
md"""
On the interval $[0, \xi_0]$ we have that $\frac{d^2}{d\xi^2}|Q|^{2}$ is contained in the interval
"""

# ╔═╡ 589ebba1-0acf-4f62-a1c2-0a8ab081016b
latexstring(abs2_Q_derivative2_ξ₀_equation_latex)

# ╔═╡ 28b0b3b6-e5c5-42ec-b38b-50354bc724b4
print(abs2_Q_derivative2_ξ₀_equation_latex)

# ╔═╡ Cell order:
# ╟─8d4fdc56-6f0b-4a20-b07b-f5ac07c235b2
# ╠═d18c2f76-73b5-11ef-27a6-6d984691fd6c
# ╠═c3053bb7-4415-40e5-8774-3e72ae6e527f
# ╟─8ebc0fbe-8309-4eb0-ad70-d09657316253
# ╠═2c73a770-ec32-4c0b-9c74-31676f85eb32
# ╟─259a5c21-1497-4343-a67c-4b8d930d8f44
# ╟─c88671ba-f5a2-4293-a3b9-9cd8e5ad7cf0
# ╟─08229fb6-98af-4581-af0a-2587e8861be2
# ╟─4583b371-c009-42cd-8956-6be3b1de3b7b
# ╠═d9078fa5-0e06-49a2-b947-54a8370e0952
# ╟─d102bce4-0e13-4644-89f9-f0044bdff092
# ╠═8487b9c4-bd77-4272-9d6b-3ab043d9e305
# ╟─2fc0db5b-0d71-4f93-8e7b-0312103eaeeb
# ╠═42d3826b-c8f6-47f1-9e4b-19d6c2e55887
# ╠═9abf6a82-301e-4065-9913-54e6da57596f
# ╠═f6a4ac88-a959-4e5e-9d3d-9f9b81033c1f
# ╟─aec77170-7845-4390-91f8-c8cba67c0115
# ╠═c9ae053a-ca79-4188-9c01-f9d8953e0cb9
# ╟─ae76061d-4435-4376-9456-0c79dfc6cd6e
# ╟─037453b2-caa8-47a2-aa48-83f08a7c948f
# ╠═0d425114-1f17-4297-ab51-1aef6df765b8
# ╟─9795fa90-9821-485a-9583-04634d0f7f0f
# ╠═467a5c66-9a07-4fbc-b738-39dce6d2be1a
# ╟─5cb7cb1f-66e6-466a-bfbc-2e1a1d47e0e0
# ╠═17eee8c0-15c3-4986-8b5c-a7422ab16ece
# ╟─09578453-f174-4957-a0a6-238c88d56a57
# ╠═7fbc5c46-2164-445a-99e0-f2d3750a11e6
# ╟─f2e891bd-7669-43e6-9ab1-286dee91d221
# ╠═06b734d2-7362-4f78-b737-68fcc94e603c
# ╠═8b61869b-5c95-4392-a07d-5de8dda85bd7
# ╠═d9137251-fa69-4942-bd70-983f75adad30
# ╟─b3f55a5b-9a79-4446-8f56-adf774634c86
# ╟─5a756d4f-297e-4f4c-bfc8-eda7063a51b9
# ╠═85b544a0-3480-4bb1-9565-ef889b686ed6
# ╠═ef43800a-59dd-4724-af00-df370792e58b
# ╟─a8c8b1e0-cf1b-462c-8ae1-6b898c35f8ab
# ╟─54af636d-7b4b-4cc8-9847-1186fd66df86
# ╠═36ebd7ae-ec94-4b43-83da-bb0fd836ce40
# ╟─fcc09ed6-29bd-45af-b731-77e9e6ccf8b8
# ╟─3b0b26c9-f0a0-4564-b635-636f40d8279e
# ╟─4c14198f-531f-4bf6-82d4-94e315034d15
# ╟─a3568944-2444-456c-b262-de32dedf7c06
# ╟─939c60ab-b10a-41d1-8f67-823840f5c728
# ╟─520c152b-3196-4d7c-9950-2e25c6898354
# ╠═82d39edb-ecbf-496c-9d09-476d369a7b24
# ╠═62a35f0e-dab9-468f-b52a-42b2c2de8b34
# ╟─9b5f5289-9256-4d99-94b7-53aad0df45dd
# ╠═8f7e5899-94da-4f5b-b2c3-b6395be86711
# ╟─2194a21f-09bb-4cff-883f-3b26f289ce45
# ╠═9fde3fa4-0c85-4ece-a5b0-e09f44ccb370
# ╟─2aeddc6f-05b8-4545-931c-08c55c5b767e
# ╠═71c8a3aa-374b-4eb1-bfd7-41d989d3ec65
# ╟─7190d92b-7f62-4805-bae1-d23e869c736b
# ╠═b3a3e9d0-3ea0-4523-8141-8ee0cb68741a
# ╟─4a12122c-c7dd-4e14-aff7-621c3badc910
# ╠═b35ec786-9c68-4d09-88c8-a3a316bbba31
# ╟─2867d01e-4839-4200-9fa9-938856d6e521
# ╠═cf7a3fab-ea03-4eab-a3d3-291611d74672
# ╟─eddb69c4-13a0-4dde-8eb1-1bc5249cde09
# ╠═2a603618-19ef-4645-8977-2ec4b2d9b630
# ╟─dbfa101a-dec9-4b27-a2b8-1c70cd2e128c
# ╠═41ef9478-0f84-4eba-bc54-9a383f15565d
# ╟─8b43f3a6-7315-433d-97ea-2e21ccadd899
# ╠═2b1fc470-59a2-4da5-9f9d-d72ba7b96216
# ╟─1181d465-27c9-41fa-89c6-dc70bca6c286
# ╠═d31e6dc5-79dc-4101-958c-0c44919ac865
# ╟─d1e6fabf-6f4f-4608-884d-7fd025d62345
# ╠═95d3cba6-1cdb-47f8-bbdf-db0620ac5ea6
# ╟─f7bce6c5-8b02-467a-ae18-7115540301d7
# ╠═e5067db1-e7db-4703-a162-ff5eafbc153b
# ╠═06473469-0f0c-483f-af67-e1b4c17652f9
# ╟─83e16fe2-6d33-4dcf-a8e1-abffc043d869
# ╟─11d9ac67-075b-47f5-bdc0-efe990dc40d2
# ╟─b2f87105-4083-469c-ba8d-556d980804ad
# ╟─38dbeb52-743b-4020-9aba-36649bbd9c48
# ╟─1a0ef6c0-293a-4f73-a013-d60adbefbfa2
# ╟─1f13b882-5e96-4fee-84aa-c6480e6f30c8
# ╟─ea662336-745c-438c-bdc2-acbfb151ed34
# ╠═861405c9-92ea-46ee-b61a-e937ebaf7a45
# ╠═ec85ad10-2ed8-4695-b2cd-e550e4994bbc
# ╟─2de7e6fd-a952-43aa-833e-f66e3749b462
# ╟─bc60cae3-001a-453c-be3b-e5eaaea38efe
# ╟─6dd19686-d49e-4acb-b301-59c74969b307
# ╟─b1ff87ee-2e7b-450c-aaea-dd76b4e11718
# ╟─c1bfa2ec-a316-4cd6-82b5-cf59b1ae2a6a
# ╟─0f012491-884f-4930-b95c-0a7c255113bf
# ╠═919b631a-c4a7-4504-b346-1a0040a490ad
# ╠═478381cc-0f8a-4f4c-9919-fac3744a5a2c
# ╟─2ddd1d45-38f9-4ecd-baef-885e8988b618
# ╟─188adc35-bbb6-41a1-b189-9a2955281abc
# ╟─66c23812-1859-405b-816d-08f36f958080
# ╠═af0f44fd-17df-4f0c-b7ca-84644fe08276
# ╟─7900e0e5-9817-403b-bbe2-9866138efb1e
# ╟─9897916c-ea93-46d4-9ffd-a74c160f1301
# ╠═df6e6400-b9c4-4494-adc2-98028ee86a6d
# ╟─5e796977-0bff-49d7-aeec-322714d9471d
# ╟─3bf6078f-9f62-4d24-bdcc-769842c78f67
# ╟─48f98796-75c4-4851-817d-558caba65bbf
# ╠═d376330e-1793-4652-815d-0a091feb9024
# ╟─b2b46e6c-f3a7-4713-b911-27a2462eadb5
# ╟─fe7e0856-a161-4431-996e-6ffb7438a4b4
# ╟─314e4e16-6e8f-46a3-8d99-2df74d73560a
# ╟─588ffcb8-2865-4c46-af76-dae8e8dc1e2a
# ╟─2e881d99-aa92-4ef6-afac-65990aff2487
# ╟─27375192-be08-48c0-a872-dbf630e03ca4
# ╠═4b9a7d4c-8828-41e3-931e-b5ac76fbf2cb
# ╟─88a12187-3dcd-4dbf-81e5-d80232d63dbf
# ╟─dbea5cae-4150-4155-910e-e52b6952019d
# ╠═f09c6a45-7116-44ac-ac7f-696b084935a1
# ╟─db3ad4c3-7c9c-4dad-a097-eb63f22ade38
# ╟─1f3a56af-070e-4204-a54f-d7aa7e9c47fc
# ╟─2e0ff47e-ea61-424b-bac5-b8022e1de9ce
# ╟─4f299a4d-a7a8-4353-a0b1-c771e6d77dfc
# ╟─94a05218-0ae9-44d1-a5f2-439c5b1cbd4c
# ╟─833eb4e0-7c7c-44f6-94e7-02d4859b7c4c
# ╟─001294b1-866c-43d7-ac9f-e6c602363d25
# ╠═a74ec90f-4673-4720-a34a-696d19015321
# ╟─3f730b42-2456-4b49-a50e-27b3acbd4573
# ╠═745c36db-56c6-46af-8b21-ea8ea8aa7a20
# ╟─d8c514e4-cf96-4431-9041-1ca4a497a399
# ╠═e0217811-4851-47a6-ac5a-e4c163a919ef
# ╟─94dbdfe3-a5ce-47c0-bf6b-c582603c10cc
# ╟─68508799-c37b-4c4f-932f-0750d915f892
# ╠═a356d6d9-a004-431c-8c85-d5ba6f79ba66
# ╠═2b043bd4-aacf-43be-8141-578bb0928f86
# ╟─4f9e593b-d884-4601-88c2-0d806c4e8998
# ╠═7b2938c0-6831-49ae-bb68-82ce688147b6
# ╟─89a29d15-8711-4168-9d9a-cf8adc86b83a
# ╠═d47f8a46-e517-4251-b30b-c8c91a51e483
# ╟─97614c71-26c6-4dc6-974d-240fd281743c
# ╠═392baa05-6032-4b6f-a1d3-15548aa59ecb
# ╟─41e0c536-08b2-4f18-86f6-a46222a097bc
# ╠═af1e9e82-14e9-479e-81d3-1f7a2a4b8a62
# ╟─57b6013f-6fe1-4798-a65d-80fc246f0eee
# ╠═83bb7eab-921a-44f6-9f81-736624982e44
# ╟─0ee00ce0-c395-4de9-aaaa-021991bf4885
# ╠═352e4f92-58e1-4d5b-a44b-52051aaf5589
# ╟─8a1909f8-a48f-4479-b0cf-3a0398020864
# ╠═b0d4bb71-9c87-4002-bee9-e5db160a29d6
# ╟─fd7e44f4-554b-4ee1-a3e9-4a965a84cb41
# ╟─624bdfbd-9547-418c-8678-86d19f5c32b4
# ╟─716f515d-fcb6-41a5-b3ab-f0aedd489836
# ╟─288e71ec-adb2-4b21-8d29-dd6ca0fb3e20
# ╠═91b6cc42-6e5d-4dec-ab34-bec5fe73be48
# ╟─ee456f31-40c3-40d4-930d-8f6c9930cfd7
# ╟─4bffedfe-f887-4540-b0c2-208509682b03
# ╟─5ff56d71-170a-494a-a890-2a2f3e5f6755
# ╟─8b5c0fbb-38da-450e-b521-e0307af5722b
# ╟─1b45e60c-7cb7-4109-85da-57eef2b27a0b
# ╠═bac5b3d7-a5e2-4323-90a9-85ea2bdfd416
# ╠═139ba5bc-f06d-485c-83c8-be141adbc80f
# ╟─b2c1830d-f694-4df6-a416-826e8d1d279a
# ╟─a214d93b-caee-42f9-b2af-088f9ad9d253
# ╠═f7b2cfa4-a59d-43d2-83aa-44abcd72bdb5
# ╠═e36acd5f-3a6f-4415-a805-9f3cfa615333
# ╟─0a62baa4-04bc-4df8-9232-1fbc76b75e21
# ╟─3049ee1c-cc91-4eb9-8eee-bb1df4e6713d
# ╟─bc1745d7-2026-42af-b18a-bf9669db543c
# ╟─3a951310-9437-4c3c-ad72-cfa9620cbb20
# ╟─4c7227d2-3322-4dad-98d0-80b7b118ae5d
# ╟─540ab5f3-226b-45ab-956a-2b371ec975af
# ╠═d9a797ed-54f5-4a10-8ab9-8c7ad6aff38a
# ╠═3bbffd81-268a-408d-a150-c9a65f074e5a
# ╠═2cf059d7-a966-4dd8-b5a2-bb3b01752e08
# ╟─89c62802-6242-4f7a-ab6e-cba6851d5ff5
# ╟─7c898b59-e6af-4197-ba06-9ed474373eb9
# ╟─589ebba1-0acf-4f62-a1c2-0a8ab081016b
# ╠═28b0b3b6-e5c5-42ec-b38b-50354bc724b4
