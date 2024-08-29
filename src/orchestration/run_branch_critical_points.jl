function run_branch_critical_points_load_existence_data(
    j::Integer,
    d::Integer,
    part,
    directory;
    N::Union{Nothing,Integer} = nothing,
    verbose::Bool = false,
)
    df, parameters = run_branch_continuation_load_data(j, d, part, directory; verbose)

    if !isnothing(N) && N < size(df, 1)
        verbose && @info "Limiting to $N subintervals"
        df = df[1:N, :]
    end

    (; ξ₁, λ) = parameters

    μs = df.μ_exists
    γs = df.γ_exists
    if part == "turn"
        κs = Arb.(tuple.(df.κ_lower, df.κ_upper))
        ϵs = df.ϵ_exists
    else
        κs = df.κ_exists
        ϵs = Arb.(tuple.(df.ϵ_lower, df.ϵ_upper))
    end
    ξ₁s = fill(ξ₁, length(μs))

    return μs, γs, κs, ϵs, ξ₁s, λ
end

function run_branch_critical_points_load_points_data(
    j::Integer,
    d::Integer,
    part,
    directory::Nothing;
    N::Union{Nothing,Integer} = nothing,
    verbose::Bool = false,
)
    verbose && @info "No directory for data given"

    base_directory = relpath(
        joinpath(
            dirname(pathof(@__MODULE__)),
            "../Dardel/output/branch_points_d=$(d)_fix_kappa",
        ),
    )

    verbose && @info "Searching for most recent data in" base_directory

    directory = joinpath(base_directory, maximum(readdir(base_directory)))

    return run_branch_critical_points_load_points_data(j, d, part, directory; N, verbose)
end

function run_branch_critical_points_load_points_data(
    j::Integer,
    d::Integer,
    part,
    directory::AbstractString;
    N::Union{Nothing,Integer} = nothing,
    verbose::Bool = false,
)
    filename = "branch_points_j=$(j)_d=$(d).csv.gz"

    verbose && @info "Loading data for branch" directory filename

    df = CGL.read_branch_points_csv(joinpath(directory, filename))

    parameters = CGL.read_parameters(joinpath(directory, "parameters.csv"))

    verbose && @info "Succesfully loaded data with $(size(df, 1)) points"

    if !isnothing(N) && N < size(df, 1)
        verbose && @info "Limiting to $N subintervals"
        df = df[1:N, :]
    end

    return df.μ, df.γ, df.κ, df.ϵ, df.ξ₁, parameters.λ
end

function run_branch_critical_points(
    j::Integer = 1,
    d::Integer = 1,
    part = "top";
    use_midpoint::Bool = false,
    use_points_data::Bool = false,
    directory_load::Union{Nothing,AbstractString} = nothing,
    N::Union{Nothing,Integer} = nothing,
    batch_size::Integer = 32,
    pool::Distributed.WorkerPool = Distributed.WorkerPool(Distributed.workers()),
    save_results::Bool = true,
    directory::Union{Nothing,AbstractString} = nothing,
    log_progress::Bool = true,
    verbose::Bool = true,
)
    verbose && @info "Computing for j = $j,  d = $d, part = $part" N directory_load

    if !use_points_data
        μs, γs, κs, ϵs, ξ₁s, λ = run_branch_critical_points_load_existence_data(
            j,
            d,
            part,
            directory_load;
            N,
            verbose,
        )
    else
        μs, γs, κs, ϵs, ξ₁s, λ = run_branch_critical_points_load_points_data(
            j,
            d,
            part,
            directory_load;
            N,
            verbose,
        )
    end

    if use_midpoint
        μs = midpoint.(Arb, μs)
        γs = midpoint.(Acb, γs)
        κs = midpoint.(Arb, κs)
        ϵs = midpoint.(Arb, ϵs)
    end

    res = branch_critical_points(μs, γs, κs, ϵs, ξ₁s, λ; batch_size, pool, verbose)

    return res
end
