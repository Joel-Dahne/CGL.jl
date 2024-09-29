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

dataframes = map([0.9, 1.0, 1.1]) do scaling
    df_1 = CGL.run_branch_points(
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

    df_2 = CGL.run_branch_points(
        d;
        fix_kappa,
        scaling,
        ξ₁_strategy = :automatic,
        ξ₁_strategy_value = [15, 20, 25, 30, 35, 40, 50, 60, 70, 80],
        N,
        batch_size = num_threads,
        pool,
        save_results = true,
        log_progress = true,
        verbose = true,
    )

    df_3 = CGL.run_branch_points(
        d;
        fix_kappa,
        scaling,
        ξ₁_strategy = :fixed,
        ξ₁_strategy_value = 30,
        N,
        batch_size = num_threads,
        pool,
        save_results = true,
        log_progress = true,
        verbose = true,
    )

    df_4 = CGL.run_branch_points(
        d;
        fix_kappa,
        scaling,
        ξ₁_strategy = :perturbed,
        ξ₁_strategy_value = 1.05,
        N,
        batch_size = num_threads,
        pool,
        save_results = true,
        log_progress = true,
        verbose = true,
    )

    [df_1, df_2, df_3, df_4]
end

dataframes = reduce(vcat, dataframes)
