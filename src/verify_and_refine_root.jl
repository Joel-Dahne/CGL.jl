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

    # Only used for heuristically stopping, so fine to do in Float64
    error_previous = radius.(Float64, root)
    isproved = false
    for i = 1:max_iterations
        mid = midpoint.(Arb, root)
        y = f(mid)
        dy = df(root)

        y_ArbMatrix = ArbMatrix(y)
        dy_ArbMatrix = ArbMatrix(dy)
        dy_div_y = let dy_div_y = similar(y_ArbMatrix)
            success = !iszero(Arblib.solve!(dy_div_y, dy_ArbMatrix, y_ArbMatrix))

            if !success
                verbose && @warn "Could not compute dy \\ y" dy
                break
            end

            convert(typeof(y), dy_div_y[:])
        end

        new_root = mid - dy_div_y

        # Note that since Arblib.intersection! only returns an
        # enclosure of the result it's not enough to check that the
        # new enclosure is contained in the interior of previous
        # enclosure. It could happen that the previous enclosure
        # contains numbers not contained in the original enclosure and
        # in that case we are not guaranteed that the root is
        # contained in the original enclosure, only in the previous
        # one.
        if any(isnan, new_root)
            verbose && @warn "Non-finite values" new_root
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
