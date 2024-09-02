function _run_branch_continuation_load_data(
    j::Integer,
    d::Integer,
    part,
    directory::Union{Nothing,AbstractString};
    N::Union{Nothing,Integer} = nothing,
    verbose::Bool = false,
)
    if isnothing(directory)
        verbose && @info "No directory for data given"

        base_directory = relpath(
            joinpath(
                dirname(pathof(@__MODULE__)),
                "../Dardel/output/branch_existence_j=$(j)_d=$(d)_part=$(part)",
            ),
        )

        verbose && @info "Searching for most recent data in" base_directory

        directory = maximum(readdir(base_directory, join = true))
    end

    filename = "branch_existence_j=$(j)_d=$(d)_part=$part.csv.gz"

    verbose && @info "Loading data for branch" directory filename

    df = CGL.read_branch_existence_csv(joinpath(directory, filename))

    parameters = CGL.read_parameters(joinpath(directory, "parameters.csv"))

    verbose && @info "Succesfully loaded data with $(size(df, 1)) subintervals"

    if !isnothing(N) && N < size(df, 1)
        verbose && @info "Limiting to $N subintervals"
        df = df[1:N, :]
    end

    if part != "turn"
        ϵs_or_κs = collect(zip(df.ϵ_lower, df.ϵ_upper))
        exists =
            CGL.SVector.(df.μ_exists, real.(df.γ_exists), imag.(df.γ_exists), df.κ_exists)
        uniqs = CGL.SVector.(df.μ_uniq, real.(df.γ_uniq), imag.(df.γ_uniq), df.κ_uniq)
        approxs =
            CGL.SVector.(df.μ_approx, real.(df.γ_approx), imag.(df.γ_approx), df.κ_approx)
    else
        ϵs_or_κs = collect(zip(df.κ_lower, df.κ_upper))
        exists =
            CGL.SVector.(df.μ_exists, real.(df.γ_exists), imag.(df.γ_exists), df.ϵ_exists)
        uniqs = CGL.SVector.(df.μ_uniq, real.(df.γ_uniq), imag.(df.γ_uniq), df.ϵ_uniq)
        approxs =
            CGL.SVector.(df.μ_approx, real.(df.γ_approx), imag.(df.γ_approx), df.ϵ_approx)
    end

    return ϵs_or_κs, exists, uniqs, approxs, parameters
end

function run_branch_continuation(
    j::Integer = 1,
    d::Integer = 1,
    part = "top";
    directory_existence::Union{Nothing,AbstractString} = nothing,
    N::Union{Nothing,Integer} = nothing,
    maxevals::Integer = 50000,
    depth::Integer = 20,
    batch_size::Integer = 32,
    pool::Distributed.WorkerPool = Distributed.WorkerPool(Distributed.workers()),
    save_results::Bool = true,
    directory::Union{Nothing,AbstractString} = nothing,
    log_progress::Bool = true,
    verbose::Bool = true,
    verbose_segments::Bool = true,
)
    verbose && @info "Computing for j = $j,  d = $d, part = $part" N directory_existence

    fix_kappa = part == "turn"

    ϵs_or_κs, exists, uniqs, approxs, parameters =
        _run_branch_continuation_load_data(j, d, part, directory_existence; N, verbose)

    (; ξ₁, λ) = parameters

    runtime_continuation = @elapsed left_continuation, ϵs_or_κs, exists, uniqs, approxs =
        CGL.branch_continuation(
            ϵs_or_κs,
            exists,
            uniqs,
            approxs,
            ξ₁,
            λ;
            batch_size,
            fix_kappa,
            verbose,
        )

    if !fix_kappa
        df = CGL.branch_continuation_dataframe_fix_epsilon(
            left_continuation,
            ϵs_or_κs,
            uniqs,
            exists,
            approxs,
        )
    else
        df = CGL.branch_continuation_dataframe_fix_kappa(
            left_continuation,
            ϵs_or_κs,
            uniqs,
            exists,
            approxs,
        )
    end

    if save_results
        if isnothing(directory)
            # Avoid colon (:) in date string since that is not allowed
            # on Windows
            date_string = replace(string(round(Dates.now(), Dates.Second)), ":" => "")
            commit_string = readchomp(`git rev-parse --short HEAD`)
            directory = relpath(
                joinpath(
                    dirname(pathof(@__MODULE__)),
                    "../Dardel/output/branch_continuation_j=$(j)_d=$(d)_part=$(part)",
                    "$(date_string)_$(commit_string)",
                ),
            )
        end

        verbose && @info "Writing data" directory

        mkpath(directory)
        CGL.write_parameters(
            joinpath(directory, "parameters.csv"),
            ξ₁,
            λ;
            runtime_continuation,
            runtime_existence = parameters.runtime_existence,
        )
        CGL.write_branch_continuation_csv(
            joinpath(directory, "branch_continuation_j=$(j)_d=$(d)_part=$part.csv.gz"),
            df,
        )
    else
        verbose && @info "Not writing data"
    end

    return df
end
