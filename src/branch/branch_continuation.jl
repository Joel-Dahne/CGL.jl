function verify_branch_continuation(
    ϵs::Vector{NTuple{2,Arf}},
    μs::Vector{Arb},
    γs::Vector{Acb},
    κs::Vector{Arb},
    ξ₁::Arb,
    λ::CGLParams{Arb};
    pool = Distributed.WorkerPool(Distributed.workers()),
    verbose = false,
    log_progress = false,
)
    @assert length(ϵs) == length(μs) == length(γs) == length(κs)

    # TODO: Implement this

    return ϵs, exists, uniqs
end
