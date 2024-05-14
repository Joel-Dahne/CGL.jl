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

# Read arguments
d = length(ARGS) > 0 ? parse(Int, ARGS[1]) : 1
fix_kappa = length(ARGS) > 1 ? parse(Bool, ARGS[2]) : false
N = length(ARGS) > 2 ? parse(Int, ARGS[3]) : 0
scaling = length(ARGS) > 3 ? parse(Float64, ARGS[4]) : 1.0

CGL.run_branch_points(
    d;
    fix_kappa,
    scaling,
    ξ₁_strategy = :default,
    ξ₁_strategy_value = nothing,
    N,
    batch_size = num_threads,
    pool,
    save_results = true,
    log_progress = true,
    verbose = true,
)
