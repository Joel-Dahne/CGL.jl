"""
    newton_step(f, df, x; verbose = false)

Perform one internval Newton iteration for the function `f` on the
input `x`. The function `df` should compute the Jacobian of `f`.
"""
function newton_step(f, df, x; verbose = false)
    mid = midpoint.(Arb, x)

    y = ArbMatrix(f(mid))

    isfinite(y) || return indeterminate.(x)

    J = ArbMatrix(df(x))
    J_div_y = similar(y)

    success = !iszero(Arblib.solve!(J_div_y, J, y))

    if success
        return mid - convert(typeof(mid), J_div_y[:])
    else
        verbose && @warn "Could not compute J \\ y" J
        return indeterminate.(x)
    end
end

"""
    verify_root(
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

The verification and refinement is done using successive interval
Newton iterations. For the method to succeed the Jacobian most be
non-zero on the enclosure of the root.

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
        new_root = newton_step(f, df, root; verbose)

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

        root = Arblib.intersection.(root, new_root)

        if !isproved && all(Arblib.contains_interior.(original_root, new_root))
            verbose && @info "Proved root"
            isproved = true
        end

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

If succesfull it returns a vector of existence and a vector of
uniqueness. The first vector is proved to contain a root of the
function, and that root is proved to be unique in the second vector.
If unsuccesful both return values are indeterminate vectors.

The method works by applying interval Newton iterations, but without
the usual intersection with the original enclosure. After each
iteration the radius of the ball is expanded by a factor determined by
`expansion_rate`. A larger `expansion_rate` generally means faster
convergence, but increases the risk of failute.

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
    original_root = root
    root = deepcopy(root)

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
        extra_verbose::Bool = false,
    )

Given an vector `root_uniqueness` proved to contain a unique root of a
function `f` (whose Jacobian is computed by `df`), attempt to widen
the region of uniqueness by widening `root_uniqueness`. The returned
result is guaranteed to also contain a unique root of `f`.

It works by expanding the initial enclosure as much as possible but so
that one Newton iteration is still enough to prove existence. For this
it searches for a maximal factor `ρ` so that when expanding the region
by this factor the Newtin iteration still succeeds.

First it finds the smallest `i` so that `ρ = 2^i` gives a failed
Newton iteration. It then does a binary search between `2^(i - 1)` and
`2^i` to find a maximal `ρ` for which it succeeds.
"""
function expand_uniqueness(
    f,
    df,
    root_uniqueness::Union{Vector{Arb},SVector{<:Any,Arb}};
    verbose::Bool = false,
    extra_verbose::Bool = false,
)
    root_mid = midpoint.(Arb, root_uniqueness)
    y = ArbMatrix(f(root_mid))

    is_ok(ρ::Mag) =
        if isone(ρ)
            true # By assumption this holds
        else
            root = add_error.(root_uniqueness, ρ * radius.(root_uniqueness))

            J = ArbMatrix(df(root))
            J_div_y = similar(y)

            success = !iszero(Arblib.solve!(J_div_y, J, y))

            if success
                new_root = root_mid - convert(typeof(root_mid), J_div_y[:])

                ok = all(Arblib.contains.(root, new_root))

                extra_verbose && @info "Newton step" ρ ok root new_root

                return ok
            else
                extra_verbose && @info "Newton step - failed inverting J" ρ
                return false
            end
        end
    # Needed for searchsortedfirst to work correctly. It applies
    # is_ok also to the true argument.
    is_ok(b::Bool) = b

    # Find smallest i such that increasing radius by the factor 2^i
    # gives a failed Newton step.
    i = findfirst(i -> !is_ok(Mag(1, i)), 1:64)

    @assert !isnothing(i) # Should never happen in practice

    verbose && @info "Smallest i for failure" i

    # Find largest 2^(i - 1) <= ρ < 2^i such that increasing the
    # radius by ρ gives a succesfull Newton step.

    # The ρ values we consider. Take them in reverse order so that the
    # first one that works is the largest one.
    ρs = reverse(Mag(1, i - 1) * Mag.(range(1, 2, 7)[2:end-1]))

    ρs_idx = searchsortedfirst(ρs, true, by = is_ok)

    ρ = if ρs_idx > lastindex(ρs)
        Mag(1, i - 1)
    else
        ρs[ρs_idx]
    end

    verbose && @info "Largest ρ for success" ρ

    return add_error.(root_uniqueness, ρ * radius.(root_uniqueness))
end
