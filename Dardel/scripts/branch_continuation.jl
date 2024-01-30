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

filename_existence = "branch_existence_j=$(j)_d=$(d)_part=$part.csv"

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

ϵs = collect(zip(df.ϵ_lower, df.ϵ_upper))
exists = CGL.SVector.(df.μ_exists, real.(df.γ_exists), imag.(df.γ_exists), df.κ_exists)
uniqs = CGL.SVector.(df.μ_uniq, real.(df.γ_uniq), imag.(df.γ_uniq), df.κ_uniq)
approxs = CGL.SVector.(df.μ_approx, real.(df.γ_approx), imag.(df.γ_approx), df.κ_approx)

if !isnothing(N) && N < size(df, 1)
    verbose && @info "Limiting to $N subintervals"
    ϵs = ϵs[1:N]
    exists = exists[1:N]
    uniqs = uniqs[1:N]
    approxs = approxs[1:N]
end

_, _, _, ξ₁, λ = CGL.sverak_params(Arb, j, d)

left_continuation, ϵs, exists, uniqs, approxs =
    CGL.verify_branch_continuation(ϵs, exists, uniqs, approxs, ξ₁, λ; verbose)

dirname = "Dardel/output/branch_continuation/$(round(Dates.now(), Second))"
filename = "branch_continuation_j=$(j)_d=$(d)_part=$part.csv"
mkpath(dirname)

verbose && @info "Writing data" dirname filename

df = CGL.branch_continuation_dataframe(
    left_continuation,
    ϵs,
    getindex.(uniqs, 1),
    Acb.(getindex.(uniqs, 2), getindex.(uniqs, 3)),
    getindex.(uniqs, 4),
    getindex.(exists, 1),
    Acb.(getindex.(exists, 2), getindex.(exists, 3)),
    getindex.(exists, 4),
    getindex.(approxs, 1),
    Acb.(getindex.(approxs, 2), getindex.(approxs, 3)),
    getindex.(approxs, 4),
    fill(ξ₁, length(ϵs)),
)

CGL.write_branch_continuation_csv(joinpath(dirname, filename), df)
