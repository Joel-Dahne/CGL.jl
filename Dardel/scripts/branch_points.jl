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
        br = CGL.CGLBranch.branch_epsilon(CGL.CGLBranch.sverak_initial(j, d)...)

        μs = Arb.(br.μ)
        κs = Arb.(br.κ)
        ϵs = Arb.(br.param)

        _, _, _, _, ξ₁, λ = CGL.sverak_params(Arb, j, d)
        ξ₁s = Arb[ξ₁ for _ = 1:length(br)]
        λs = [λ for _ = 1:length(br)]

        μs, κs, ϵs, ξ₁s, λs
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
    ϵs = foldl(vcat, getindex.(values, 3))
    ξ₁s = foldl(vcat, getindex.(values, 4))
    λs = foldl(vcat, getindex.(values, 5))

    return parameter_indices, μs, κs, ϵs, ξ₁s, λs
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

fix_kappa = if length(ARGS) > 1
    parse(Bool, ARGS[2])
else
    false
end

N = if length(ARGS) > 2
    parse(Int, ARGS[3])
else
    nothing
end

@assert isnothing(N) || N >= 2 # Need N >= 2 for range to work below

verbose && @info "Computing for d = $d" N fix_kappa

verbose && @info "Computing initial branches"

parameter_indices, μ₀s, κ₀s, ϵ₀s, ξ₁s, λs = initial_branches(pool, parameters)

verbose && @info "Got $(length(μ₀s)) branch points"

if !isnothing(N) && N < length(μ₀s)
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

verified_points = CGL.branch_points(
    μ₀s,
    κ₀s,
    ϵ₀s,
    ξ₁s,
    λs,
    batch_size = num_threads, # IMPROVE: How to pick this?
    log_progress = true;
    fix_kappa,
    verbose,
)

dirname = "Dardel/output/branch_points/$(round(Dates.now(), Second))"

verbose && @info "Writing data" dirname

if !fix_kappa
    dfs = map(parameters) do parameter
        let idxs = parameter_indices[parameter]
            CGL.branch_points_dataframe_fix_epsilon(
                verified_points[idxs],
                μ₀s[idxs],
                κ₀s[idxs],
                ϵ₀s[idxs],
                ξ₁s[idxs],
            )
        end
    end
else
    dfs = map(parameters) do parameter
        let idxs = parameter_indices[parameter]
            CGL.branch_points_dataframe_fix_kappa(
                verified_points[idxs],
                μ₀s[idxs],
                κ₀s[idxs],
                ϵ₀s[idxs],
                ξ₁s[idxs],
            )
        end
    end
end

mkpath(dirname)
for ((j, d), df) in zip(parameters, dfs)
    fix_kappa_str = fix_kappa ? "_fix_kappa" : ""
    filename = "branch_points_j=$(j)_d=$(d)$(fix_kappa_str).csv"
    CGL.write_branch_points_csv(joinpath(dirname, filename), df)
end
