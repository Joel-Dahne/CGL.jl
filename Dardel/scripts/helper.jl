using ClusterManagers, Distributed, ProgressLogging

@assert haskey(ENV, "LD_LIBRARY_PATH")

"""
    create_workers(num_workers, num_threads; use_slurm, verbose = false)

Create `num_workers` workers, each using `num_threads` threads.
"""
function create_workers(
    num_workers = parse(Int, get(ENV, "CGL_WORKERS", get(ENV, "SLURM_NTASKS", "1"))),
    num_threads = parse(Int, get(ENV, "CGL_THREADS", get(ENV, "SLURM_CPUS_PER_TASK", "1")));
    use_slurm = haskey(ENV, "SLURM_JOB_ID"),
    verbose = false,
)
    if verbose
        @info "Preparing to set up workers" num_workers num_threads
    end

    # Launch worker processes
    ENV["JULIA_PROJECT"] = Base.active_project()
    ENV["JULIA_NUM_THREADS"] = num_threads
    ENV["JULIA_WORKER_TIMEOUT"] = 300

    if use_slurm
        # Give the current sysimage explicitly. The SlurmManager
        # doesn't handle it by itself.
        #sysimage = unsafe_string(Base.JLOptions().image_file)
        #exeflags = "--sysimage=$sysimage"
        #addprocs(SlurmManager(num_workers); exeflags)
        addprocs(SlurmManager(num_workers))
    else
        addprocs(num_workers)
    end

    if verbose
        @info "Finished setting up workers" nworkers()

        # For each worker gets its id, process id, hostname and number of threads
        calls = Dict(
            w => @spawnat(w, (myid(), getpid(), gethostname(), Threads.nthreads())) for
            w in workers()
        )
        for w in workers()
            id, pid, host, threads = fetch(calls[w])
            @info "Worker $w" id pid host threads
        end
    end

    return Distributed.WorkerPool(workers()), num_threads
end

function read_args()
    j = if length(ARGS) > 0
        parse(Int, ARGS[1])
    else
        3
    end

    d = if length(ARGS) > 1
        parse(Int, ARGS[2])
    else
        1
    end

    start = if length(ARGS) > 2
        parse(Int, ARGS[3])
    else
        1
    end

    stop = if length(ARGS) > 3
        parse(Int, ARGS[4])
    else
        nothing
    end

    return j, d, start, stop
end
