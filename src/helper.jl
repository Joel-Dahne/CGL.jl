function fetch_with_progress(
    tasks::Vector{Task},
    with_progress = true;
    timeout = 5,
    pollint = 1,
    debug = parse(Bool, get(ENV, "CGL_DEBUG", "0")),
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

_complex(re::Real, im::Real) = complex(re, im)
_complex(re::Arb, im::Real) = Acb(re, im)
_complex(re::Real, im::Arb) = Acb(re, im)
_complex(re::Arb, im::Arb) = Acb(re, im)
_complex(::Type{T}) where {T<:Real} = complex(T)
_complex(::Type{Arb}) = Acb

"""
    args_to_complex(μ, γ_real, γ_imag, κ_or_ϵ)

Helper function to covert from real to complex representation of `γ`.
It returns `μ, γ, κ_or_ϵ` where `γ = γ_real + im * γ_imag`.

See also [`args_to_real(μ, γ, κ_or_ϵ)`](@ref) for the inverse.
"""
function _args_to_complex(μ::T, γ_real::T, γ_imag::T, κ_or_ϵ::T) where {T<:Real}
    return μ, _complex(γ_real, γ_imag), κ_or_ϵ
end
