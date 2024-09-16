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
j, d, part, N, _ = read_args()

CGL.run_branch_continuation(
    j,
    d,
    part;
    N,
    batch_size = 8num_threads, # IMPROVE: How to pick this?;
    pool,
    save_results = true,
    log_progress = true,
    verbose = true,
)
