function construct_proof_witness(
    filename_top::AbstractString,
    filename_turn::AbstractString,
    filename_bottom::AbstractString,
    ξ₁::Arb,
    λ::CGLParams{Arb},
)
    data_top = read_branch_continuation_csv(filename_top)
    data_turn = read_branch_continuation_csv(filename_turn)
    data_bottom = read_branch_continuation_csv(filename_bottom)

    ## Sanity check data

    # Continuation should have succeeded for all of them
    all(data_top.left_continuation) || error("left continuation failed for top")
    all(data_turn.left_continuation) || error("left continuation failed for turn")
    all(data_bottom.left_continuation) || @warn "left continuation failed for bottom"

    iszero(data_top.ϵ_lower[1]) || error("expected ϵ to start at 0 for top")

    DataFrames.select!(
        data_top,
        [
            "ϵ_lower",
            "ϵ_upper",
            "μ_uniq",
            "γ_uniq",
            "κ_uniq",
            "μ_exists",
            "γ_exists",
            "κ_exists",
        ],
    )

    DataFrames.select!(
        data_turn,
        [
            "κ_lower",
            "κ_upper",
            "μ_uniq",
            "γ_uniq",
            "ϵ_uniq",
            "μ_exists",
            "γ_exists",
            "ϵ_exists",
        ],
    )

    DataFrames.select!(
        data_bottom,
        [
            "ϵ_lower",
            "ϵ_upper",
            "μ_uniq",
            "γ_uniq",
            "κ_uniq",
            "μ_exists",
            "γ_exists",
            "κ_exists",
        ],
    )

    parameters = (; ξ₁, λ)

    return parameters, data_top, data_turn, data_bottom
end

function write_proof_witness(
    directory,
    parameters,
    data_top::DataFrame,
    data_turn::DataFrame,
    data_bottom::DataFrame,
)
    mkpath(directory)

    write_parameters(joinpath(directory, "parameters.csv"), parameters.ξ₁, parameters.λ)

    write_branch_generic(joinpath(directory, "top.csv.gz"), data_top, compress = true)
    write_branch_generic(joinpath(directory, "turn.csv.gz"), data_turn, compress = true)
    write_branch_generic(joinpath(directory, "bottom.csv.gz"), data_bottom, compress = true)

    return
end

function read_proof_witness(directory::AbstractString)
    parameters = read_parameters(joinpath(directory, "parameters.csv"))

    types = [String, String, String, String, String, String, String, String, String, String]

    data_top = let
        data_top_dump = CSV.read(joinpath(directory, "top.csv.gz"), DataFrame; types)

        data_top = DataFrame()

        data_top.ϵ_lower = Arblib.load_string.(Arf, data_top_dump.ϵ_lower_dump)
        data_top.ϵ_upper = Arblib.load_string.(Arf, data_top_dump.ϵ_upper_dump)

        data_top.μ_uniq = Arblib.load_string.(Arb, data_top_dump.μ_uniq_dump)
        data_top.γ_uniq =
            Acb.(
                Arblib.load_string.(Arb, data_top_dump.γ_uniq_dump_real),
                Arblib.load_string.(Arb, data_top_dump.γ_uniq_dump_imag),
            )
        data_top.κ_uniq = Arblib.load_string.(Arb, data_top_dump.κ_uniq_dump)

        data_top.μ_exists = Arblib.load_string.(Arb, data_top_dump.μ_exists_dump)
        data_top.γ_exists =
            Acb.(
                Arblib.load_string.(Arb, data_top_dump.γ_exists_dump_real),
                Arblib.load_string.(Arb, data_top_dump.γ_exists_dump_imag),
            )
        data_top.κ_exists = Arblib.load_string.(Arb, data_top_dump.κ_exists_dump)

        data_top
    end

    data_turn = let
        data_turn_dump = CSV.read(joinpath(directory, "turn.csv.gz"), DataFrame; types)

        data_turn = DataFrame()

        data_turn.κ_lower = Arblib.load_string.(Arf, data_turn_dump.κ_lower_dump)
        data_turn.κ_upper = Arblib.load_string.(Arf, data_turn_dump.κ_upper_dump)

        data_turn.μ_uniq = Arblib.load_string.(Arb, data_turn_dump.μ_uniq_dump)
        data_turn.γ_uniq =
            Acb.(
                Arblib.load_string.(Arb, data_turn_dump.γ_uniq_dump_real),
                Arblib.load_string.(Arb, data_turn_dump.γ_uniq_dump_imag),
            )
        data_turn.ϵ_uniq = Arblib.load_string.(Arb, data_turn_dump.ϵ_uniq_dump)

        data_turn.μ_exists = Arblib.load_string.(Arb, data_turn_dump.μ_exists_dump)
        data_turn.γ_exists =
            Acb.(
                Arblib.load_string.(Arb, data_turn_dump.γ_exists_dump_real),
                Arblib.load_string.(Arb, data_turn_dump.γ_exists_dump_imag),
            )
        data_turn.ϵ_exists = Arblib.load_string.(Arb, data_turn_dump.ϵ_exists_dump)

        data_turn
    end

    data_bottom = let
        data_bottom_dump = CSV.read(joinpath(directory, "bottom.csv.gz"), DataFrame; types)

        data_bottom = DataFrame()

        data_bottom.ϵ_lower = Arblib.load_string.(Arf, data_bottom_dump.ϵ_lower_dump)
        data_bottom.ϵ_upper = Arblib.load_string.(Arf, data_bottom_dump.ϵ_upper_dump)

        data_bottom.μ_uniq = Arblib.load_string.(Arb, data_bottom_dump.μ_uniq_dump)
        data_bottom.γ_uniq =
            Acb.(
                Arblib.load_string.(Arb, data_bottom_dump.γ_uniq_dump_real),
                Arblib.load_string.(Arb, data_bottom_dump.γ_uniq_dump_imag),
            )
        data_bottom.κ_uniq = Arblib.load_string.(Arb, data_bottom_dump.κ_uniq_dump)

        data_bottom.μ_exists = Arblib.load_string.(Arb, data_bottom_dump.μ_exists_dump)
        data_bottom.γ_exists =
            Acb.(
                Arblib.load_string.(Arb, data_bottom_dump.γ_exists_dump_real),
                Arblib.load_string.(Arb, data_bottom_dump.γ_exists_dump_imag),
            )
        data_bottom.κ_exists = Arblib.load_string.(Arb, data_bottom_dump.κ_exists_dump)

        data_bottom
    end

    return parameters, data_top, data_turn, data_bottom
end

function check_proof_witness_part(data, reversed)
    param_lower, param_upper = names(data)[1:2]

    # Check that the subintervals are consecutive
    if !reversed
        for i = 1:nrow(data)-1
            isequal(data[i, param_upper], data[i+1, param_lower]) ||
                error("intervals not consecutive")
        end
    else
        for i = 1:nrow(data)-1
            isequal(data[i, param_lower], data[i+1, param_upper]) ||
                error("intervals not consecutive")
        end
    end

    # Read uniqueness and existence enclosures
    if "κ_uniq" in names(data)
        uniqs = SVector.(data.μ_uniq, real.(data.γ_uniq), imag.(data.γ_uniq), data.κ_uniq)
        exists =
            SVector.(
                data.μ_exists,
                real.(data.γ_exists),
                imag.(data.γ_exists),
                data.κ_exists,
            )
    else
        uniqs = SVector.(data.μ_uniq, real.(data.γ_uniq), imag.(data.γ_uniq), data.ϵ_uniq)
        exists =
            SVector.(
                data.μ_exists,
                real.(data.γ_exists),
                imag.(data.γ_exists),
                data.ϵ_exists,
            )
    end

    # Check that all enclosures are finite and that the existence is
    # strictly contained in the uniqueness
    for i in eachindex(uniqs, exists)
        all(isfinite, uniqs[i]) || error("non-finite uniqueness")
        all(isfinite, exists[i]) || error("non-finite existence")
        all(Arblib.contains_interior.(uniqs[i], exists[i])) ||
            error("exists not conained in uniqs")
    end

    # Check that we have unique continuation
    # TODO: Turn into error once we have good data
    failed = 0
    for i = 2:length(uniqs)
        failed += !all(Arblib.contains.(uniqs[i-1], exists[i]))
    end
    iszero(failed) || @warn "no unique continuation for $failed subintervals"

    return
end

# TODO: Turn warning into errors once we have good data
function check_proof_witness_top_turn(parameters, data_top, data_turn)
    # Check that the two branches reach far enough in ϵ and κ

    # The last subinterval for the top should be strictly larger than
    # the enclosure in ϵ of the first subinterval for the turn
    data_top.ϵ_lower[end] > data_turn.ϵ_exists[1] ||
        @warn "top and turn not overlapping in ϵ"

    # The first subinterval for the turn should be strictly smaller
    # than the enclosure in κ of the last subinterval for the top
    data_turn.κ_upper[1] > data_top.κ_exists[end] ||
        @warn "top and turn not overlapping in κ"

    # Check that the two branches are the same

    # For the last subinterval in the top find the first subinterval
    # in the turn so that their enclosures overlaps.
    exists_top = SVector(
        data_top.μ_exists[end],
        real(data_top.γ_exists[end]),
        imag(data_top.γ_exists[end]),
        data_top.κ_exists[end],
        Arb((data_top.ϵ_lower[end], data_top.ϵ_upper[end])),
    )

    i = findfirst(eachindex(eachrow(data_turn))) do i
        all(
            Arblib.overlaps.(
                SVector(
                    data_turn.μ_exists[i],
                    real(data_turn.γ_exists[i]),
                    imag(data_turn.γ_exists[i]),
                    Arb.(zip(data_turn.κ_lower[i], data_turn.κ_upper[i])),
                    data_turn.ϵ_exists[i],
                ),
                exists_top,
            ),
        )
    end

    exists_turn = SVector(
        data_turn.μ_exists[i],
        real(data_turn.γ_exists[i]),
        imag(data_turn.γ_exists[i]),
        Arb.(zip(data_turn.κ_lower[i], data_turn.κ_upper[i])),
        data_turn.ϵ_exists[i],
    )

    ϵ = midpoint(Arb, Arblib.intersection.(exists_top[5], exists_turn[5]))

    (; ξ₁, λ) = parameters

    μ, γ, κ = refine_approximation(exists_top[1], exists_top[4], ϵ, ξ₁, λ)

    all(Arblib.contains_interior.(exists_top, SVector(μ, real(γ), imag(γ), κ, ϵ))) ||
        error("picked point not contained in top")

    all(Arblib.contains_interior.(exists_turn, SVector(μ, real(γ), imag(γ), κ, ϵ))) ||
        @warn "picked point not contained in turn"

    exists = G_solve(μ, real(γ), imag(γ), κ, ϵ, ξ₁, λ, verbose = true)

    all(isfinite.(exists)) || error("non-finite enclosure for picked point")

    all(Arblib.contains_interior.(exists_top, SVector(exists..., ϵ))) ||
        error("enclosure for picked point not contained in top")
    all(Arblib.contains_interior.(exists_turn, SVector(exists..., ϵ))) ||
        @warn "enclosure for picked point not contained in top"
end

# TODO: Turn warning into errors once we have good data
function check_proof_witness_turn_bottom(parameters, data_turn, data_bottom)
    # Check that the two branches reach far enough in ϵ and κ

    # The last subinterval for the turn should be strictly smaller
    # than the enclosure in κ of the first subinterval for the bottom
    data_turn.κ_lower[end] < data_bottom.κ_exists[1] ||
        @warn "turn and bottom not overlapping in κ"

    # The first subinterval for the bottom should be strictly larger
    # than the enclosure in ϵ of the last subinterval for the turn
    data_bottom.ϵ_lower[1] > data_turn.ϵ_exists[end] ||
        @warn "turn and bottom not overlapping in ϵ"

    # Check that the two branches are the same

    # For the first subinterval in the botom find the last subinterval
    # in the turn so that their enclosures overlaps.
    exists_bottom = SVector(
        data_bottom.μ_exists[1],
        real(data_bottom.γ_exists[1]),
        imag(data_bottom.γ_exists[1]),
        data_bottom.κ_exists[1],
        Arb((data_bottom.ϵ_lower[1], data_bottom.ϵ_upper[1])),
    )

    i = findlast(eachindex(eachrow(data_turn))) do i
        all(
            Arblib.overlaps.(
                SVector(
                    data_turn.μ_exists[i],
                    real(data_turn.γ_exists[i]),
                    imag(data_turn.γ_exists[i]),
                    Arb.(zip(data_turn.κ_lower[i], data_turn.κ_upper[i])),
                    data_turn.ϵ_exists[i],
                ),
                exists_bottom,
            ),
        )
    end

    exists_turn = SVector(
        data_turn.μ_exists[i],
        real(data_turn.γ_exists[i]),
        imag(data_turn.γ_exists[i]),
        Arb.(zip(data_turn.κ_lower[i], data_turn.κ_upper[i])),
        data_turn.ϵ_exists[i],
    )

    ϵ = midpoint(Arb, Arblib.intersection.(exists_bottom[5], exists_turn[5]))

    (; ξ₁, λ) = parameters

    μ, γ, κ = refine_approximation(exists_bottom[1], exists_bottom[4], ϵ, ξ₁, λ)

    all(Arblib.contains_interior.(exists_bottom, SVector(μ, real(γ), imag(γ), κ, ϵ))) ||
        error("picked point not contained in bottom")

    all(Arblib.contains_interior.(exists_turn, SVector(μ, real(γ), imag(γ), κ, ϵ))) ||
        @warn "picked point not contained in turn"

    exists = G_solve(μ, real(γ), imag(γ), κ, ϵ, ξ₁, λ)

    all(isfinite.(exists)) || error("non-finite enclosure for picked point")

    all(Arblib.contains_interior.(exists_bottom, SVector(exists..., ϵ))) ||
        error("enclosure for picked point not contained in bottom")
    all(Arblib.contains_interior.(exists_turn, SVector(exists..., ϵ))) ||
        @warn "enclosure for picked point not contained in turn"
end

function check_proof_witness(parameters, data_top, data_turn, data_bottom)
    @info "Checking top part"
    check_proof_witness_part(data_top, false)

    @info "Checking turn part"
    check_proof_witness_part(data_turn, true)

    @info "Checking bottom part"
    #check_proof_witness_part(data_bottom, true)

    @info "Checking connection between top and turn"
    check_proof_witness_top_turn(parameters, data_top, data_turn)

    @info "Checking connection between turn and bottom"
    check_proof_witness_turn_bottom(parameters, data_turn, data_bottom)

    return true
end
