export @assert_proof

"""
    ProofError <: Exception

Thrown when [`@assert_proof`](@ref) could not verify the given
condition.
"""
struct ProofError <: Exception end

Base.showerror(io::IO, e::ProofError) = print(io, "ProofError: could not verify condition")

"""
    @assert_proof cond

Return `true` if `cond` is `true`, otherwise throw a
[`ProofError`](@ref). This is used to indicate conditions that we want
to prove are true.

The main purpose of having a separate macro for this is to make it
more explicit which conditions directly correspond to a statement in
the accompanying paper.
"""
macro assert_proof(ex)
    return :($(esc(ex)) ? $(true) : throw(ProofError()))
end
