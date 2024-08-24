function run_branch_continuation_load_data(
    j::Integer,
    d::Integer,
    part,
    directory_existence::Nothing;
    verbose::Bool = false,
)
    verbose && @info "No directory for existence data given"

    base_directory_existence = relpath(
        joinpath(
            dirname(pathof(@__MODULE__)),
            "../Dardel/output/branch_existence_j=$(j)_d=$(d)_part=$(part)",
        ),
    )

    verbose && @info "Searching for most recent data in" base_directory_existence

    directory_existence =
        joinpath(base_directory_existence, maximum(readdir(base_directory_existence)))

    return run_branch_continuation_load_data(j, d, part, directory_existence; verbose)
end

function run_branch_continuation_load_data(
    j::Integer,
    d::Integer,
    part,
    directory_existence::AbstractString;
    verbose::Bool = false,
)
    fix_kappa = part == "turn"

    filename_existence = "branch_existence_j=$(j)_d=$(d)_part=$part.csv.gz"

    verbose &&
        @info "Loading existence data for branch" directory_existence filename_existence

    df_existence =
        CGL.read_branch_existence_csv(joinpath(directory_existence, filename_existence))

    parameters = CGL.read_parameters(joinpath(directory_existence, "parameters.csv"))

    verbose &&
        @info "Succesfully loaded existence data with $(size(df_existence, 1)) subintervals"

    return df_existence, parameters
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

    df_existence, parameters =
        run_branch_continuation_load_data(j, d, part, directory_existence; verbose)

    (; ξ₁, λ) = parameters
    runtime_existence = parameters.runtime # Only for storing in output

    if !isnothing(N) && N < size(df_existence, 1)
        verbose && @info "Limiting to $N subintervals"
        df_existence = df_existence[1:N, :]
    end

    if !fix_kappa
        ϵs_or_κs = collect(zip(df_existence.ϵ_lower, df_existence.ϵ_upper))
        exists =
            CGL.SVector.(
                df_existence.μ_exists,
                real.(df_existence.γ_exists),
                imag.(df_existence.γ_exists),
                df_existence.κ_exists,
            )
        uniqs =
            CGL.SVector.(
                df_existence.μ_uniq,
                real.(df_existence.γ_uniq),
                imag.(df_existence.γ_uniq),
                df_existence.κ_uniq,
            )
        approxs =
            CGL.SVector.(
                df_existence.μ_approx,
                real.(df_existence.γ_approx),
                imag.(df_existence.γ_approx),
                df_existence.κ_approx,
            )
    else
        ϵs_or_κs = collect(zip(df_existence.κ_lower, df_existence.κ_upper))
        exists =
            CGL.SVector.(
                df_existence.μ_exists,
                real.(df_existence.γ_exists),
                imag.(df_existence.γ_exists),
                df_existence.ϵ_exists,
            )
        uniqs =
            CGL.SVector.(
                df_existence.μ_uniq,
                real.(df_existence.γ_uniq),
                imag.(df_existence.γ_uniq),
                df_existence.ϵ_uniq,
            )
        approxs =
            CGL.SVector.(
                df_existence.μ_approx,
                real.(df_existence.γ_approx),
                imag.(df_existence.γ_approx),
                df_existence.ϵ_approx,
            )
    end

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
            runtime_existence,
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
