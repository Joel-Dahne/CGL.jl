function initial_branches_helper(j::Integer, d::Integer)
    br = CGLBranch.branch_epsilon(CGLBranch.sverak_initial(j, d)...)

    μs = Arb.(br.μ)
    κs = Arb.(br.κ)
    ϵs = Arb.(br.param)

    _, _, _, _, ξ₁, λ = sverak_params(Arb, j, d)
    ξ₁s = Arb[ξ₁ for _ = 1:length(br)]
    λs = [λ for _ = 1:length(br)]

    μs, κs, ϵs, ξ₁s, λs
end

function initial_branches(pool, parameters, scaling)
    if length(pool) > 1
        # Distribute over workers
        tasks = map(parameters) do (j, d)
            @async Distributed.remotecall_fetch(initial_branches_helper, pool, j, d)
        end

        values = CGL.fetch_with_progress(tasks)
    else
        # Distribute over threads
        values = tmap(parameters) do (j, d)
            initial_branches_helper(j, d)
        end
    end

    endpoints = [0; cumsum(length.(getindex.(values, 1)))]

    parameter_indices =
        Dict(parameters[i] => endpoints[i]+1:endpoints[i+1] for i in eachindex(parameters))

    μs = foldl(vcat, getindex.(values, 1))
    κs = foldl(vcat, getindex.(values, 2))
    ϵs = foldl(vcat, getindex.(values, 3))
    ξ₁s = foldl(vcat, getindex.(values, 4))
    λs = foldl(vcat, getindex.(values, 5))

    if !isone(scaling)
        # Apply the scaling. No γ to scale.
        for i in eachindex(μs, κs, ϵs, ξ₁s, λs)
            μs[i], _, κs[i], ϵs[i], ξ₁s[i], λs[i] = scale_params(
                μs[i],
                indeterminate(Acb),
                κs[i],
                ϵs[i],
                ξ₁s[i],
                λs[i];
                scaling,
            )
        end
    end

    return parameter_indices, μs, κs, ϵs, ξ₁s, λs
end

function run_branch_points(
    d::Integer = 1;
    fix_kappa::Bool = false,
    scaling = 1,
    ξ₁_strategy = :default,
    ξ₁_strategy_value = nothing,
    N = 0,
    batch_size::Integer = Threads.nthreads(),
    pool = Distributed.WorkerPool(Distributed.workers()),
    directory = nothing,
    save_results::Bool = true,
    log_progress::Bool = true,
    verbose::Bool = true,
)
    js = if d == 1
        1:8
    elseif d == 3
        1:5
    else
        throw(ArgumentError("only supports d = 1 and d = 3"))
    end

    parameters = tuple.(js, d)

    verbose &&
        @info "Computing for d = $d" scaling fix_kappa ξ₁_strategy ξ₁_strategy_value N

    verbose && @info "Computing initial branches"

    parameter_indices, μ₀s, κ₀s, ϵ₀s, ξ₁_defaults, λs =
        initial_branches(pool, parameters, scaling)

    @assert allequal(λs)
    λ = λs[1]

    verbose && @info "Got $(length(μ₀s)) branch points"

    if ξ₁_strategy == :default
        @assert isnothing(ξ₁_strategy_value)
        ξ₁s = ξ₁_defaults
    elseif ξ₁_strategy == :automatic
        @assert ξ₁_strategy_value isa AbstractVector
        ξ₁s = fill(Arb.(ξ₁_strategy_value), length(ξ₁_defaults))
    elseif ξ₁_strategy == :fixed
        @assert ξ₁_strategy_value isa Real
        ξ₁s = fill(Arb(ξ₁_strategy_value), length(ξ₁_defaults))
    elseif ξ₁_strategy == :perturbed
        @assert ξ₁_strategy_value isa Real
        ξ₁s = map(ξ₁ -> Arb(ξ₁_strategy_value * ξ₁), ξ₁_defaults)
    else
        throw(ArgumentError("unkown ξ₁ strategy $ξ₁_strategy"))
    end

    # TODO: Handle this better
    if d == 3
        ξ₁s = fill(Arb.(10:5:80), length(ξ₁s))
    end

    if !iszero(N) && N < length(μ₀s)
        N >= 2 || throw(ArgumentError("Must have N >= 2, got N = $N"))

        verbose && @info "Limiting to $N branch points"

        idxs = round.(Int, range(1, length(μ₀s), N))
        μ₀s, κ₀s, ϵ₀s, ξ₁s, λs = μ₀s[idxs], κ₀s[idxs], ϵ₀s[idxs], ξ₁s[idxs], λs[idxs]

        # Make sure parameter_indices only contains valid indices
        parameter_indices = let
            new_lengths = map(parameters) do parameter
                length(intersect(parameter_indices[parameter], idxs))
            end
            new_endpoints = [0; cumsum(new_lengths)]

            Dict(
                parameters[i] => new_endpoints[i]+1:new_endpoints[i+1] for
                i in eachindex(parameters)
            )
        end
    end

    verbose && @info "Verifying branch points"

    verified_points, ξ₁s_used = branch_points(
        μ₀s,
        κ₀s,
        ϵ₀s,
        ξ₁s,
        λs,
        log_progress = true;
        batch_size,
        fix_kappa,
        verbose,
    )

    if !fix_kappa
        dfs = map(parameters) do parameter
            let idxs = parameter_indices[parameter]
                branch_points_dataframe_fix_epsilon(
                    verified_points[idxs],
                    μ₀s[idxs],
                    κ₀s[idxs],
                    ϵ₀s[idxs],
                    ξ₁s_used[idxs],
                )
            end
        end
    else
        dfs = map(parameters) do parameter
            let idxs = parameter_indices[parameter]
                branch_points_dataframe_fix_kappa(
                    verified_points[idxs],
                    μ₀s[idxs],
                    κ₀s[idxs],
                    ϵ₀s[idxs],
                    ξ₁s_used[idxs],
                )
            end
        end
    end

    if save_results
        if isnothing(directory)
            date_string = round(Dates.now(), Dates.Second)
            commit_string = readchomp(`git rev-parse --short HEAD`)
            directory = relpath(
                joinpath(
                    dirname(pathof(@__MODULE__)),
                    "../Dardel/output/branch_points_d=$d",
                    "$(date_string)_$(commit_string)",
                ),
            )
        end

        verbose && @info "Writing data" directory

        mkpath(directory)
        CGL.write_parameters(
            joinpath(directory, "parameters.csv"),
            indeterminate(Arb), # ξ₁ varies, write an indeterminate value
            λ;
            fix_kappa,
            ξ₁_strategy,
            ξ₁_strategy_value = repr(ξ₁_strategy_value),
        )
        for ((j, d), df) in zip(parameters, dfs)
            filename = "branch_points_j=$(j)_d=$(d).csv.gz"
            write_branch_points_csv(joinpath(directory, filename), df)
        end
    else
        verbose && @info "Not writing data"
    end

    return dfs, λ
end
