function _run_branch_critical_points_load_data(
    j::Integer,
    d::Integer,
    part,
    directory::Union{Nothing,AbstractString};
    N::Union{Nothing,Integer} = nothing,
    verbose::Bool = false,
)
    if isnothing(directory)
        directory = locate_most_recent("continuation", j, d, part; verbose)

        isnothing(directory) && return Arb[], Acb[], Arb[], Arb[], Arb[], missing
    end

    filename = "branch_continuation_j=$(j)_d=$(d)_part=$part.csv.gz"

    verbose && @info "Loading data for branch" directory filename

    df = CGL.read_branch_continuation_csv(joinpath(directory, filename))

    parameters = CGL.read_parameters(joinpath(directory, "parameters.csv"))

    verbose && @info "Succesfully loaded data with $(size(df, 1)) subintervals"

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

function _run_branch_critical_points_load_data(
    j::Integer,
    d::Integer,
    directory;
    N::Union{Nothing,Integer} = nothing,
    verbose::Bool = false,
)
    μs_top, γs_top, κs_top, ϵs_top, ξ₁s_top, λ_top =
        _run_branch_critical_points_load_data(j, d, "top", directory; verbose)

    μs_turn, γs_turn, κs_turn, ϵs_turn, ξ₁s_turn, λ_turn =
        _run_branch_critical_points_load_data(j, d, "turn", directory; verbose)

    μs_bottom, γs_bottom, κs_bottom, ϵs_bottom, ξ₁s_bottom, λ_bottom =
        _run_branch_critical_points_load_data(j, d, "bottom", directory; verbose)

    μs = vcat(μs_top, μs_turn, μs_bottom)
    γs = vcat(γs_top, γs_turn, γs_bottom)
    κs = vcat(κs_top, κs_turn, κs_bottom)
    ϵs = vcat(ϵs_top, ϵs_turn, ϵs_bottom)
    ξ₁s = vcat(ξ₁s_top, ξ₁s_turn, ξ₁s_bottom)

    if !isnothing(N) && N < length(μs)
        verbose && @info "Limiting to $N subintervals"
        μs = μs[1:N]
        γs = γs[1:N]
        κs = κs[1:N]
        ϵs = ϵs[1:N]
        ξ₁s = ξ₁s[1:N]
    end

    @assert !ismissing(λ_top)
    @assert ismissing(λ_turn) || isequal(λ_top, λ_turn)
    @assert ismissing(λ_bottom) || isequal(λ_top, λ_bottom)

    return μs, γs, κs, ϵs, ξ₁s, λ_top
end

function _run_branch_critical_points_load_points_data(
    j::Integer,
    d::Integer,
    directory::Union{Nothing,AbstractString};
    N::Union{Nothing,Integer} = nothing,
    verbose::Bool = false,
)
    if isnothing(directory)
        verbose && @info "No directory for data given"

        base_directory = relpath(
            joinpath(
                dirname(pathof(@__MODULE__)),
                "../Dardel/output/branch_points_d=$(d)_fix_kappa",
            ),
        )

        verbose && @info "Searching for most recent data in" base_directory

        directory = maximum(readdir(base_directory, join = true))
    end

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

    if part == "full"
        μs, γs, κs, ϵs, ξ₁s, λ =
            _run_branch_critical_points_load_data(j, d, directory_load; N, verbose)
    elseif part == "points"
        μs, γs, κs, ϵs, ξ₁s, λ =
            _run_branch_critical_points_load_points_data(j, d, directory_load; N, verbose)
    else
        μs, γs, κs, ϵs, ξ₁s, λ =
            _run_branch_critical_points_load_data(j, d, part, directory_load; N, verbose)
    end

    if use_midpoint
        μs = midpoint.(Arb, μs)
        γs = midpoint.(Acb, γs)
        κs = midpoint.(Arb, κs)
        ϵs = midpoint.(Arb, ϵs)
    end

    num_critical_points =
        branch_critical_points(μs, γs, κs, ϵs, ξ₁s, λ; batch_size, pool, verbose)

    df = branch_critical_points_dataframe(μs, γs, κs, ϵs, ξ₁s, num_critical_points)

    if save_results
        if isnothing(directory)
            # Avoid colon (:) in date string since that is not allowed
            # on Windows
            date_string = replace(string(round(Dates.now(), Dates.Second)), ":" => "")
            commit_string = readchomp(`git rev-parse --short HEAD`)
            directory = relpath(
                joinpath(
                    dirname(pathof(@__MODULE__)),
                    "../Dardel/output/branch_critical_points_j=$(j)_d=$(d)_part=$(part)",
                    "$(date_string)_$(commit_string)",
                ),
            )
        end

        verbose && @info "Writing data" directory

        mkpath(directory)
        CGL.write_parameters(
            joinpath(directory, "parameters.csv"),
            ifelse(allequal(ξ₁s), ξ₁s[1], indeterminate(ξ₁s[1])),
            λ;
            use_midpoint,
        )
        CGL.write_branch_existence_csv(
            joinpath(directory, "branch_critical_points_j=$(j)_d=$(d)_$(part).csv.gz"),
            df,
        )
    else
        verbose && @info "Not writing data"
    end

    return df
end
