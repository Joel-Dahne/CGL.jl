function fetch_with_progress(
    tasks::Vector{Task},
    with_progress = true;
    timeout = 5,
    pollint = 1,
    debug = true,
    debug_interval = 60,
)
    values = similar(tasks, Any)

    if with_progress
        last_debug_time = zero(time())
        remaining = fill(true, length(tasks))
        ndone = 0
        @withprogress while any(remaining)
            first_remaining = findfirst(remaining)
            timedwait(() -> istaskdone(tasks[first_remaining]), timeout; pollint)
            for i in findall(remaining)
                if istaskdone(tasks[i])
                    values[i] = fetch(tasks[i])
                    remaining[i] = false
                    ndone += 1
                    @logprogress ndone / length(values)
                end
            end

            if debug &&
               ((time() - last_debug_time) > debug_interval || ndone == length(values))

                debug_info = readchomp(`vmstat -S M -t`)
                println(stderr, debug_info)
                flush(stderr)
                last_debug_time = time()
            end
        end
    else
        for i in eachindex(tasks, values)
            values[i] = fetch(tasks[i])
        end
    end

    # Broadcasting the identity function makes it so that the eltype
    # of the vector is narrowed to the best possible.
    return identity.(values)
end
