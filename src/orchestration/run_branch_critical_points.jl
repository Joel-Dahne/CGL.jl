function run_branch_critical_points(
    j::Integer = 1,
    d::Integer = 1,
    part = "top";
    use_midpoint::Bool = false,
    directory_existence::Union{Nothing,AbstractString} = nothing,
    N::Union{Nothing,Integer} = nothing,
    batch_size::Integer = 32,
    pool::Distributed.WorkerPool = Distributed.WorkerPool(Distributed.workers()),
    save_results::Bool = true,
    directory::Union{Nothing,AbstractString} = nothing,
    log_progress::Bool = true,
    verbose::Bool = true,
)
    verbose && @info "Computing for j = $j,  d = $d, part = $part" N directory_existence

    fix_kappa = part == "turn"

    df_existence, parameters =
        run_branch_continuation_load_data(j, d, part, directory_existence; verbose)

    (; ξ₁, λ) = parameters

    if !isnothing(N) && N < size(df_existence, 1)
        verbose && @info "Limiting to $N subintervals"
        df_existence = df_existence[1:N, :]
    end

    μs = df_existence.μ_exists
    γs = df_existence.γ_exists
    if fix_kappa
        κs = Arb.(tuple.(df_existence.κ_lower, df_existence.κ_upper))
        ϵs = df_existence.ϵ_exists
    else
        κs = df_existence.κ_exists
        ϵs = Arb.(tuple.(df_existence.ϵ_lower, df_existence.ϵ_upper))
    end

    if use_midpoint
        μs = midpoint.(Arb, μs)
        γs = midpoint.(Acb, γs)
        κs = midpoint.(Arb, κs)
        ϵs = midpoint.(Arb, ϵs)
    end

    res = branch_critical_points(μs, γs, κs, ϵs, ξ₁, λ; batch_size, pool, verbose)

    return res
end
