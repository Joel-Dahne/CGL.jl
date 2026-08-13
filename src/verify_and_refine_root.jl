"""
    newton_step(f, df, x::AbstractVector{T}) where {T<:Union{Arb,Acb}}

Perform one Krawczyk internval Newton iteration for the function `f`
on the input `x`. The function `df` should compute the Jacobian of
`f`.
"""
function newton_step(f, df, x::AbstractVector{Arb})
    mid = midpoint.(Arb, x)

    y = f(mid)

    all(isfinite, y) || return indeterminate.(x)

    dy = df(x)

    all(isfinite, dy) || return indeterminate.(x)

    # Compute preconditioner in Float64, it only needs to be an
    # approximate inverse.
    Y = try
        Arb.(inv(Float64.(dy)))
    catch e
        # Inversion of Float64 matrices throws a SingularException on
        # failure. We catch this and return an indeterminate value in
        # that case.
        e isa SingularException || rethrow(e)
        return indeterminate.(x)
    end

    return mid - Y * y + (one(dy) - Y * dy) * (x - mid)
end

"""
    verify_and_refine_root(
        f,
        df,
        root::Union{Vector{Arb},SVector{<:Any,Arb}};
        atol = 0,
        rtol = 4eps(one(first(root))),
        min_iterations = 1,
        max_iterations = 20,
        verbose::Bool = false,
    )

Verify that `root` contains a root of the function `f` and refine the
enclosure. The function `df` should compute the Jacobian of `f`.

If succesfull, it returns an enclosure of existence, that is proved to
contain a unique root of the function. If unsuccesful, it returns an
indeterminate enclosure.

The verification and refinement is done using successive interval
Newton iterations.

At each iteration it checks if the required tolerance is met according
to [`ArbExtras.check_tolerance`](@ref). The default tolerances are
`atol = 0` and `rtol = 4eps(one(root))`.

If the absolute error does not improve by at least a factor 1.1
between two iterations then it stops early since it's unlikely that
further iterations would improve the result much. This can for example
happen when the function is computed with too low precision. This
check is only done if `min_iterations` have been performed, this is to
avoid stopping early due to slow convergence in the beginning. To
avoid this check entirely `min_iterations` can be put to the same as
`max_iterations`.

If `verbose = true` then print the enclosure at each iteration and
some more information in the end.
"""
function verify_and_refine_root(
    f,
    df,
    root::Union{Vector{Arb},SVector{<:Any,Arb}};
    atol = 0,
    rtol = 4eps(one(first(root))),
    min_iterations = 1,
    max_iterations = 10,
    verbose::Bool = false,
)
    original_root = root
    root = copy(root)

    verbose && @info "Original enclosure" root

    if any(!isfinite, root)
        verbose && @warn "Non-finite input"
        return indeterminate.(root)
    end

    # Only used for heuristically stopping, so fine to do in Float64
    error_previous = radius.(Float64, root)
    isproved = false
    for i = 1:max_iterations
        new_root = newton_step(f, df, root)

        # Note that since Arblib.intersection! only returns an
        # enclosure of the result it's not enough to check that the
        # new enclosure is contained in the interior of previous
        # enclosure. It could happen that the previous enclosure
        # contains numbers not contained in the original enclosure and
        # in that case we are not guaranteed that the root is
        # contained in the original enclosure, only in the previous
        # one.
        if !all(isfinite, new_root)
            verbose && @warn "Newton step failed" new_root
            break
        end
        if !all(Arblib.overlaps.(root, new_root))
            verbose && @warn "New iteration doesn't overlap" new_root
            break
        end
        if all(Arblib.contains.(new_root, root))
            if verbose && isproved
                @info "New iteration contains previous - stopping early" new_root
            elseif verbose
                @warn "New iteration contains previous - stopping early" new_root
            end
            break
        end

        # Checking all(Arblib.contains_interior.(root, new_root)) proves the
        # existence of a root in new_root.
        # Checking Arblib.contains.(original_root, new_root) proves
        # that the found root is in the original enclosure. (This is
        # otherwise not guaranteed due to rounding in
        # Arblib.intersection).
        if !isproved &&
           all(Arblib.contains_interior.(root, new_root)) &&
           all(Arblib.contains.(original_root, new_root))

            verbose && @info "Proved root"
            isproved = true
        end

        root = Arblib.intersection.(root, new_root)

        if verbose && isproved
            @info "Enclosure" root
        elseif verbose
            @info "Iteration" new_root
        end

        # If the result satisfies the required tolerance - break
        if isproved && all(ArbExtras.check_tolerance.(root; atol, rtol))
            verbose && @info "Tolerance satisfied"
            break
        end

        # If the result did not improve compared to the last iteration
        # and we have performed the minimum number of iterations -
        # break
        error = radius.(Float64, root)
        min_improvement_factor = isproved ? 1.1 : 1.001
        if i >= min_iterations && all(min_improvement_factor * error .> error_previous)
            verbose &&
                @info "Diameter only improved by less than a factor $min_improvement_factor - stopping early" error_previous error
            break
        end
        error_previous = error
    end

    if isproved
        return root
    else
        verbose && @warn "Could not prove root"
        return indeterminate.(root)
    end
end

"""
    verify_root_from_approximation(
        f,
        df,
        root::Union{Vector{Arb},SVector{<:Any,Arb}};
        expansion_rate = 0.05
        max_iterations = 10,
        verbose::Bool = false,
    )

Given an approximation `root` of a root of the function `f` this
method attempts to prove the existence of a nearby root. The function
`df` should compute the Jacobian of `f`.

If succesfull, it returns an enclosure of existence and an enclosure
of uniqueness. The first enclosure is proved to contain a root of the
function, and that root is proved to be unique in the second
enclosure. If unsuccesful both return values are set to indeterminate
balls.

The method works by applying interval Newton iterations, but without
the usual intersection with the original enclosure. After each
iteration the radius of the ball is expanded by a factor determined by
`expansion_rate`. A larger `expansion_rate` generally means faster
convergence, but increases the risk of failute.

It returns immediately upon successfully proving the existence of a
root, it does not attempt to refine the enclosure further.

If `verbose = true` then print the result at each iteration and some
more information in the end.
"""
function verify_root_from_approximation(
    f,
    df,
    root::Union{Vector{Arb},SVector{<:Any,Arb}};
    expansion_rate = 0.05,
    max_iterations = 10,
    verbose::Bool = false,
)
    verbose && @info "Original approximation" root

    if any(!isfinite, root)
        verbose && @warn "Non-finite input"
        return indeterminate.(root), indeterminate.(root)
    end

    for i = 1:max_iterations
        new_root = newton_step(f, df, root)

        if !all(isfinite, new_root)
            verbose && @warn "Newton step failed" new_root
            return indeterminate.(root), indeterminate.(root)
        end

        if all(Arblib.contains_interior.(root, new_root))
            root_uniqueness = root
            root = new_root
            verbose && @info "Success" root root_uniqueness
            return root, root_uniqueness
        end

        verbose && @info "Iteration $i" new_root

        root = add_error.(new_root, expansion_rate * radius.(new_root))
    end

    verbose && @warn "Reached maximum number of iterations"
    return indeterminate.(root), indeterminate.(root)
end

"""
    expand_uniqueness(
        f,
        df,
        root_uniqueness::Union{Vector{Arb},SVector{<:Any,Arb}};
        verbose::Bool = false,
    )

Given an vector `root_uniqueness` proved to contain a unique root of a
function `f` (whose Jacobian is computed by `df`), attempt to widen
the region of uniqueness by widening `root_uniqueness`. The returned
result is guaranteed to also contain a unique root of `f`.

It works by expanding the initial enclosure as much as possible but so
that the interval enclosure of the Jacobian is still invertible. For
this it searches for a maximal factor `ρ` so that when expanding the
region by this factor the Newtin iteration still succeeds.

First it finds the smallest `i` so that `ρ = 2^i` gives a
non-invertible enclosure. It then does a binary search between `2^(i -
1)` and `2^i` to find a maximal `ρ` for which it succeeds.

Note that `f` is never actually used in the computation. It is kept as
an argument so that the function signature is similar to
[`verify_root_from_approximation`](@ref).
"""
function expand_uniqueness(
    f,
    df,
    root_uniqueness::Union{Vector{Arb},SVector{<:Any,Arb}};
    verbose::Bool = false,
)
    root_mid = midpoint.(root_uniqueness)
    root_radius = radius.(root_uniqueness)

    is_ok(ρ::Mag) =
        if isone(ρ)
            true # By assumption the root is unique in the initial enclosure
        else
            root = setball.(Arb, root_mid, ρ * root_radius)
            # Check if enclosure of Jacobian is invertible by checking
            # if determinant is non-zero.
            !Arblib.contains_zero(det(ArbMatrix(df(root))))
        end
    # Needed for searchsortedfirst to work correctly. It applies
    # is_ok also to the true argument.
    is_ok(b::Bool) = b

    # Find smallest i such that increasing radius by the factor 2^i
    # gives a non-invertible Jacobian enclosure.
    is = 1:64
    i = findfirst(i -> !is_ok(Mag(1, i)), is)

    verbose && @info "Smallest i for failure" i

    if isnothing(i)
        # All tested i values worked. Take the largest one.
        return setball.(Arb, root_mid, Mag(1, is[end]) * root_radius)
    end

    # Find largest 2^(i - 1) <= ρ < 2^i such that increasing the
    # radius by ρ gives a succesfull Newton step.

    # Take ρ values strictly between 2^(i - 1) and 2^i. Take them in
    # reverse order so that the first one that works is the largest
    # one.
    ρs = reverse(Mag.(range(1, 2, 7)[2:(end-1)])) * Mag(1, i - 1)

    ρs_idx = searchsortedfirst(ρs, true, by = is_ok)

    ρ = if ρs_idx > lastindex(ρs)
        # None of the values strictly between 2^(i - 1) and 2^i
        # worked, take 2^(i - 1).
        Mag(1, i - 1)
    else
        ρs[ρs_idx]
    end

    verbose && @info "Largest ρ for success" ρ

    if isone(ρ)
        return root_uniqueness # No improvement
    else
        return setball.(Arb, root_mid, ρ * root_radius)
    end
end
