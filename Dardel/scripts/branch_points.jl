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

    function initial_branches_helper(j, d)
        (μ, γ, κ, ξ₁, λ) = CGL.sverak_params(Arb, j, d)

        # We always want to use ξ₁ = 30 here
        br = let λ = CGL.CGLBranch.Params(λ.ϵ, λ.d, λ.ω, λ.σ, λ.δ, 30.0)
            CGL.CGLBranch.branch(Float64(μ), Float64(κ), λ)
        end

        μs = Arb.(br.branch.μ)
        κs = Arb.(br.branch.κ)
        ξ₁s = Arb[ξ₁ for _ = 1:length(br.branch)]
        λs = [CGLParams(λ; ϵ) for ϵ in br.branch.param]

        μs, κs, ξ₁s, λs
    end
end

function initial_branches(pool, parameters)
    tasks = map(parameters) do (j, d)
        @async Distributed.remotecall_fetch(initial_branches_helper, pool, j, d)
    end

    values = CGL.fetch_with_progress(tasks)

    endpoints = [0; cumsum(length.(getindex.(values, 1)))]

    parameter_indices =
        Dict(parameters[i] => endpoints[i]+1:endpoints[i+1] for i in eachindex(parameters))

    μs = foldl(vcat, getindex.(values, 1))
    κs = foldl(vcat, getindex.(values, 2))
    ξ₁s = foldl(vcat, getindex.(values, 3))
    λs = foldl(vcat, getindex.(values, 4))

    return parameter_indices, μs, κs, ξ₁s, λs
end

parameters = [
    (j = 1, d = 1),
    (j = 2, d = 1),
    (j = 3, d = 1),
    (j = 4, d = 1),
    (j = 5, d = 1),
    (j = 6, d = 1),
    (j = 7, d = 1),
    (j = 8, d = 1),
    (j = 1, d = 3),
    (j = 2, d = 3),
    (j = 3, d = 3),
    (j = 4, d = 3),
    (j = 5, d = 3),
]

verbose = true

# Read arguments
d = if length(ARGS) > 0
    parse(Int, ARGS[1])
else
    1
end

@assert d == 1 || d == 3
filter!(p -> p.d == d, parameters)

N = if length(ARGS) > 1
    parse(Int, ARGS[2])
else
    nothing
end

@assert isnothing(N) || N > 1

verbose && @info "Computing for d = $d" N

verbose && @info "Computing initial branches"

parameter_indices, μs_approx, κs_approx, ξ₁s, λs = initial_branches(pool, parameters)

verbose && @info "Got $(length(μs_approx)) branch points"

if !isnothing(N) && N < length(μs_approx)
    verbose && @info "Limiting to $N branch points"

    idxs = round.(Int, range(1, length(μs_approx), N))
    μs_approx, κs_approx, ξ₁s, λs = μs_approx[idxs], κs_approx[idxs], ξ₁s[idxs], λs[idxs]

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

verified_points = CGL.verify_branch_points(
    μs_approx,
    κs_approx,
    ξ₁s,
    λs,
    batch_size = num_threads, # IMPROVE: Should this be higher?
    log_progress = true;
    verbose,
)

dirname = "Dardel/output/branch_points/$(round(Dates.now(), Second))"
mkpath(dirname)

verbose && @info "Writing data" dirname

dfs = map(parameters) do parameter
    let idxs = parameter_indices[parameter]
        CGL.branch_points_dataframe(
            λs[idxs],
            verified_points[idxs],
            μs_approx[idxs],
            κs_approx[idxs],
            ξ₁s[idxs],
        )
    end
end

for ((j, d), df) in zip(parameters, dfs)
    filename = "branch_points_j=$(j)_d=$d.csv"
    CGL.write_branch_points_csv(joinpath(dirname, filename), df)
end
