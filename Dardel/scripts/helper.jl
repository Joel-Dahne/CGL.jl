using ClusterManagers, Distributed

@assert haskey(ENV, "LD_LIBRARY_PATH")
@assert contains(ENV["LD_LIBRARY_PATH"], "capd")

"""
    create_workers(num_workers, num_threads; use_slurm, verbose = false)

Create `num_workers` workers, each using `num_threads` threads.
"""
function create_workers(
    num_workers = parse(Int, get(ENV, "CGL_WORKERS", get(ENV, "SLURM_NTASKS", "1"))),
    num_threads = parse(Int, get(ENV, "CGL_THREADS", get(ENV, "SLURM_CPUS_PER_TASK", "1")));
    heap_size_hint_G = nothing,
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

    if isnothing(heap_size_hint_G)
        if haskey(ENV, "CGL_HEAP_SIZE_HINT")
            heap_size_hint_G = parse(Int, ENV["CGL_HEAP_SIZE_HINT"])
        elseif use_slurm && haskey(ENV, "CGL_SLURM_MEM_PER_NODE")
            mem_per_node_G = parse(Int, ENV["CGL_SLURM_MEM_PER_NODE"])
            # Leave some room for other things, hence the -1
            heap_size_hint_G = mem_per_node_G ÷ num_workers - 1
        end
    end

    exeflags = String[]

    if !isnothing(heap_size_hint_G)
        verbose && @info "Using heap size hint" heap_size_hint_G
        push!(exeflags, "--heap-size-hint=$(heap_size_hint_G)G")
    end

    if use_slurm
        addprocs(SlurmManager(num_workers); exeflags)
    else
        addprocs(num_workers; exeflags)
    end

    if verbose
        @info "Finished setting up workers" nworkers()

        # For each worker gets its id, process id, hostname and number of threads
        calls = Dict(
            w => @spawnat(
                w,
                (myid(), getpid(), gethostname(), Threads.nthreads(), Threads.ngcthreads())
            ) for w in workers()
        )
        for w in workers()
            id, pid, host, threads, gcthreads = fetch(calls[w])
            @info "Worker $w" id pid host threads gcthreads
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

    part = if length(ARGS) > 2
        ARGS[3]
    else
        "top"
    end

    part == "top" ||
        part == "turn" ||
        part == "bottom" ||
        throw(ArgumentError("unknown part type $part"))

    N = if length(ARGS) > 3
        parse(Int, ARGS[4])
    else
        nothing
    end

    return j, d, part, N
end
