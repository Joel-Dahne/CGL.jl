using Dates

include("helper.jl")

pool, num_threads = create_workers(verbose = true)

@everywhere begin
    using Arblib, CGL
    setprecision(Arb, 128)

    # Set logging to always flush
    using Logging: global_logger
    using TerminalLoggers: TerminalLogger
    global_logger(TerminalLogger(always_flush = true))
end

verbose = true

# Read arguments
j, d, part, N = read_args()

dirname_existence = nothing

verbose && @info "Determined arguments" j d part N dirname_existence

fix_kappa = part == "turn"

filename_existence = "branch_existence_j=$(j)_d=$(d)_part=$part.csv.gz"

if isnothing(dirname_existence)
    verbose && @info "No directory for existence data given"

    base_dirname = "Dardel/output/branch_existence/"
    verbose && @info "Searching for most recent data in" base_dirname filename_existence

    dirnames = sort(readdir(base_dirname))

    i = findlast(dirnames) do dirname
        in(filename_existence, readdir(joinpath(base_dirname, dirname)))
    end

    dirname_existence = joinpath(base_dirname, dirnames[i])
end

verbose && @info "Loading existence data for branch" dirname_existence filename_existence

df = CGL.read_branch_existence_csv(joinpath(dirname_existence, filename_existence))

verbose && @info "Succesfully loaded existence data with $(size(df, 1)) subintervals"

if !isnothing(N) && N < size(df, 1)
    verbose && @info "Limiting to $N subintervals"
    df = df[1:N, :]
end

_, _, _, _, ξ₁, λ = CGL.sverak_params(Arb, j, d)

if !fix_kappa
    ϵs_or_κs = collect(zip(df.ϵ_lower, df.ϵ_upper))
    exists = CGL.SVector.(df.μ_exists, real.(df.γ_exists), imag.(df.γ_exists), df.κ_exists)
    uniqs = CGL.SVector.(df.μ_uniq, real.(df.γ_uniq), imag.(df.γ_uniq), df.κ_uniq)
    approxs = CGL.SVector.(df.μ_approx, real.(df.γ_approx), imag.(df.γ_approx), df.κ_approx)
else
    ϵs_or_κs = collect(zip(df.κ_lower, df.κ_upper))
    exists = CGL.SVector.(df.μ_exists, real.(df.γ_exists), imag.(df.γ_exists), df.ϵ_exists)
    uniqs = CGL.SVector.(df.μ_uniq, real.(df.γ_uniq), imag.(df.γ_uniq), df.ϵ_uniq)
    approxs = CGL.SVector.(df.μ_approx, real.(df.γ_approx), imag.(df.γ_approx), df.ϵ_approx)
end

runtime =
    @elapsed left_continuation, ϵs_or_κs, exists, uniqs, approxs = CGL.branch_continuation(
        ϵs_or_κs,
        exists,
        uniqs,
        approxs,
        ξ₁,
        λ,
        batch_size = 8num_threads; # IMPROVE: How to pick this?
        fix_kappa,
        verbose,
    )

dirname = "Dardel/output/branch_continuation/$(round(Dates.now(), Second))"
filename = "branch_continuation_j=$(j)_d=$(d)_part=$part.csv.gz"

verbose && @info "Writing data" dirname filename

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

mkpath(dirname)
CGL.write_parameters(joinpath(dirname, "parameters.csv"), ξ₁, λ; runtime)
CGL.write_branch_continuation_csv(joinpath(dirname, filename), df)
