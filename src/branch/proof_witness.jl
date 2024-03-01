function construct_proof_witness(
    filename_top::AbstractString,
    filename_turn::AbstractString,
    filename_bottom::AbstractString,
)
    @info "Reading data"
    data_top = read_branch_continuation_csv(filename_top)
    data_turn = read_branch_continuation_csv(filename_turn)
    data_bottom = read_branch_continuation_csv(filename_bottom)

    parameters_top = read_parameters(joinpath(dirname(filename_top), "parameters.csv"))
    parameters_turn = read_parameters(joinpath(dirname(filename_turn), "parameters.csv"))
    parameters_bottom =
        read_parameters(joinpath(dirname(filename_bottom), "parameters.csv"))

    ξ₁ = parameters_top.ξ₁
    λ = parameters_top.λ

    runtime_existence_top = parameters_top.runtime_existence
    runtime_existence_turn = parameters_turn.runtime_existence
    runtime_existence_bottom = parameters_bottom.runtime_existence

    runtime_continuation_top = parameters_top.runtime_continuation
    runtime_continuation_turn = parameters_turn.runtime_continuation
    runtime_continuation_bottom = parameters_bottom.runtime_continuation

    ## Sanity check data

    # Parameters should be the same
    isequal(ξ₁, parameters_turn.ξ₁) || error("turn doesn't have same ξ₁ as top")
    isequal(ξ₁, parameters_bottom.ξ₁) || error("bottom doesn't have same ξ₁ as top")
    isequal(λ, parameters_turn.λ) || error("turn doesn't have same λ as top")
    isequal(λ, parameters_bottom.λ) || error("bottom doesn't have same λ as top")

    # Continuation should have succeeded for all of them
    all(data_top.left_continuation) || error("left continuation failed for top")
    all(data_turn.left_continuation) || error("left continuation failed for turn")
    all(data_bottom.left_continuation) || error("left continuation failed for bottom")

    iszero(data_top.ϵ_lower[1]) || error("expected ϵ to start at 0 for top")

    ## Construct output

    parameters = (;
        ξ₁,
        λ,
        runtime_existence_top,
        runtime_existence_turn,
        runtime_existence_bottom,
        runtime_continuation_top,
        runtime_continuation_turn,
        runtime_continuation_bottom,
    )

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

    @info "Checking top part"
    check_proof_witness_part(data_top, false)

    @info "Checking turn part"
    check_proof_witness_part(data_turn, true)

    @info "Checking bottom part"
    check_proof_witness_part(data_bottom, true)

    # Find points connecting the segments into one
    @info "Generating connection points"
    ϵ_top, exists_top, uniq_top =
        construct_proof_witness_connect(parameters, data_top, data_turn, true)

    ϵ_bottom, exists_bottom, uniq_bottom =
        construct_proof_witness_connect(parameters, data_bottom, data_turn, false)

    data_connection_points = DataFrame(
        ϵ = [ϵ_top, ϵ_bottom],
        μ_uniq = [uniq_top[1], uniq_bottom[1]],
        γ_uniq = [Acb(uniq_top[2], uniq_top[3]), Acb(uniq_bottom[2], uniq_bottom[3])],
        κ_uniq = [uniq_top[4], uniq_bottom[4]],
        μ_exists = [exists_top[1], exists_bottom[1]],
        γ_exists = [
            Acb(exists_top[2], exists_top[3]),
            Acb(exists_bottom[2], exists_bottom[3]),
        ],
        κ_exists = [exists_top[4], exists_bottom[4]],
    )

    return parameters, data_top, data_turn, data_bottom, data_connection_points
end

function construct_proof_witness_connect(
    parameters,
    data_top_or_bottom,
    data_turn,
    is_top;
    verbose = false,
)
    # For the top we want to look at the list interval, for the bottom
    # we want to look at the first.
    i = ifelse(is_top, nrow(data_top_or_bottom), 1)

    exists_top_or_bottom = SVector(
        data_top_or_bottom.μ_exists[i],
        real(data_top_or_bottom.γ_exists[i]),
        imag(data_top_or_bottom.γ_exists[i]),
        data_top_or_bottom.κ_exists[i],
        Arb((data_top_or_bottom.ϵ_lower[i], data_top_or_bottom.ϵ_upper[i])),
    )

    # For the top, find the find the first subinterval in the turn
    # that overlaps. For the bottom find the last.
    subinterval_overlaps(j) = all(
        Arblib.overlaps.(
            SVector(
                data_turn.μ_exists[j],
                real(data_turn.γ_exists[j]),
                imag(data_turn.γ_exists[j]),
                Arb.(zip(data_turn.κ_lower[j], data_turn.κ_upper[j])),
                data_turn.ϵ_exists[j],
            ),
            exists_top_or_bottom,
        ),
    )

    j = if is_top
        findfirst(subinterval_overlaps, eachindex(eachrow(data_turn)))
    else
        findlast(subinterval_overlaps, eachindex(eachrow(data_turn)))
    end

    exists_turn = SVector(
        data_turn.μ_exists[j],
        real(data_turn.γ_exists[j]),
        imag(data_turn.γ_exists[j]),
        Arb.(zip(data_turn.κ_lower[j], data_turn.κ_upper[j])),
        data_turn.ϵ_exists[j],
    )

    (; ξ₁, λ) = parameters
    ϵ = midpoint(Arb, Arblib.intersection.(exists_top_or_bottom[5], exists_turn[5]))
    μ, γ, κ =
        refine_approximation(exists_top_or_bottom[1], exists_top_or_bottom[4], ϵ, ξ₁, λ)

    verbose && @info "Connecting through point" μ γ κ ϵ

    all(
        Arblib.contains_interior.(exists_top_or_bottom, SVector(μ, real(γ), imag(γ), κ, ϵ)),
    ) || error("picked point not contained in top/bottom")

    all(Arblib.contains_interior.(exists_turn, SVector(μ, real(γ), imag(γ), κ, ϵ))) ||
        error("picked point not contained in turn")

    exists, uniq =
        G_solve(μ, real(γ), imag(γ), κ, ϵ, ξ₁, λ, return_uniqueness = Val{true}())

    verbose && @info "Enclosure at point" exists

    all(isfinite.(exists)) || error("non-finite enclosure for picked point")

    all(Arblib.contains_interior.(exists_top_or_bottom, SVector(exists..., ϵ))) ||
        error("enclosure for picked point not contained in top/bottom")
    all(Arblib.contains_interior.(exists_turn, SVector(exists..., ϵ))) ||
        error("enclosure for picked point not contained in top")

    return ϵ, exists, uniq
end


function write_proof_witness(
    directory,
    parameters,
    data_top::DataFrame,
    data_turn::DataFrame,
    data_bottom::DataFrame,
    data_connection_points::DataFrame,
)
    mkpath(directory)

    write_parameters(
        joinpath(directory, "parameters.csv"),
        parameters.ξ₁,
        parameters.λ;
        Base.structdiff(parameters, NamedTuple{(:ξ₁, :λ)})...,
    )

    write_branch_generic(joinpath(directory, "top.csv.gz"), data_top, compress = true)
    write_branch_generic(joinpath(directory, "turn.csv.gz"), data_turn, compress = true)
    write_branch_generic(joinpath(directory, "bottom.csv.gz"), data_bottom, compress = true)
    write_branch_generic(
        joinpath(directory, "connection_points.csv.gz"),
        data_connection_points,
        compress = true,
    )

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

    data_connection_points = let
        types_connection_points =
            [String, String, String, String, String, String, String, String, String]

        data_connection_points_dump = CSV.read(
            joinpath(directory, "connection_points.csv.gz"),
            DataFrame,
            types = types_connection_points,
        )

        data_connection_points = DataFrame()

        data_connection_points.ϵ =
            Arblib.load_string.(Arb, data_connection_points_dump.ϵ_dump)

        data_connection_points.μ_uniq =
            Arblib.load_string.(Arb, data_connection_points_dump.μ_uniq_dump)
        data_connection_points.γ_uniq =
            Acb.(
                Arblib.load_string.(Arb, data_connection_points_dump.γ_uniq_dump_real),
                Arblib.load_string.(Arb, data_connection_points_dump.γ_uniq_dump_imag),
            )
        data_connection_points.κ_uniq =
            Arblib.load_string.(Arb, data_connection_points_dump.κ_uniq_dump)

        data_connection_points.μ_exists =
            Arblib.load_string.(Arb, data_connection_points_dump.μ_exists_dump)
        data_connection_points.γ_exists =
            Acb.(
                Arblib.load_string.(Arb, data_connection_points_dump.γ_exists_dump_real),
                Arblib.load_string.(Arb, data_connection_points_dump.γ_exists_dump_imag),
            )
        data_connection_points.κ_exists =
            Arblib.load_string.(Arb, data_connection_points_dump.κ_exists_dump)

        data_connection_points
    end

    return parameters, data_top, data_turn, data_bottom, data_connection_points
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
    failed = 0
    for i = 2:length(uniqs)
        failed += !all(Arblib.contains.(uniqs[i-1], exists[i]))
    end
    iszero(failed) || error("no unique continuation for $failed subintervals")

    return
end

function check_proof_witness_top_turn(
    parameters,
    data_top,
    data_turn,
    data_connection_points,
)
    # Check that the two branches reach far enough in ϵ and κ

    # The last subinterval for the top should be strictly larger than
    # the enclosure in ϵ of the first subinterval for the turn
    data_top.ϵ_lower[end] > data_turn.ϵ_exists[1] ||
        error("top and turn not overlapping in ϵ")

    # The first subinterval for the turn should be strictly smaller
    # than the enclosure in κ of the last subinterval for the top
    data_turn.κ_upper[1] > data_top.κ_exists[end] ||
        error("top and turn not overlapping in κ")

    # Check that the first point in data_connection_points is
    # contained in both branches

    point = SVector(
        data_connection_points.μ_exists[1],
        real(data_connection_points.γ_exists[1]),
        imag(data_connection_points.γ_exists[1]),
        data_connection_points.κ_exists[1],
        data_connection_points.ϵ[1],
    )

    in_top = any(1:nrow(data_top)) do i
        all(
            Arblib.contains_interior.(
                SVector(
                    data_top.μ_uniq[i],
                    real(data_top.γ_uniq[i]),
                    imag(data_top.γ_uniq[i]),
                    data_top.κ_uniq[i],
                    Arb((data_top.ϵ_lower[i], data_top.ϵ_upper[i])),
                ),
                point,
            ),
        )
    end

    in_top || error("connection point not contained in top branch")

    in_turn = any(1:nrow(data_turn)) do i
        all(
            Arblib.contains_interior.(
                SVector(
                    data_turn.μ_uniq[i],
                    real(data_turn.γ_uniq[i]),
                    imag(data_turn.γ_uniq[i]),
                    Arb((data_turn.κ_lower[i], data_turn.κ_upper[i])),
                    data_turn.ϵ_uniq[i],
                ),
                point,
            ),
        )
    end

    in_turn || error("connection point not contained in turn branch")
end

function check_proof_witness_turn_bottom(
    parameters,
    data_turn,
    data_bottom,
    data_connection_points,
)
    # Check that the two branches reach far enough in ϵ and κ

    # The last subinterval for the turn should be strictly smaller
    # than the enclosure in κ of the first subinterval for the bottom
    data_turn.κ_lower[end] < data_bottom.κ_exists[1] ||
        error("turn and bottom not overlapping in κ")

    # The first subinterval for the bottom should be strictly larger
    # than the enclosure in ϵ of the last subinterval for the turn
    data_bottom.ϵ_lower[1] > data_turn.ϵ_exists[end] ||
        error("turn and bottom not overlapping in ϵ")

    # Check that the first point in data_connection_points is
    # contained in both branches

    point = SVector(
        data_connection_points.μ_exists[2],
        real(data_connection_points.γ_exists[2]),
        imag(data_connection_points.γ_exists[2]),
        data_connection_points.κ_exists[2],
        data_connection_points.ϵ[2],
    )

    in_bottom = any(1:nrow(data_bottom)) do i
        all(
            Arblib.contains_interior.(
                SVector(
                    data_bottom.μ_uniq[i],
                    real(data_bottom.γ_uniq[i]),
                    imag(data_bottom.γ_uniq[i]),
                    data_bottom.κ_uniq[i],
                    Arb((data_bottom.ϵ_lower[i], data_bottom.ϵ_upper[i])),
                ),
                point,
            ),
        )
    end

    in_bottom || error("connection point not contained in bottom branch")

    in_turn = any(1:nrow(data_turn)) do i
        all(
            Arblib.contains_interior.(
                SVector(
                    data_turn.μ_uniq[i],
                    real(data_turn.γ_uniq[i]),
                    imag(data_turn.γ_uniq[i]),
                    Arb((data_turn.κ_lower[i], data_turn.κ_upper[i])),
                    data_turn.ϵ_uniq[i],
                ),
                point,
            ),
        )
    end

    in_turn || error("connection point not contained in turn branch")
end

function check_proof_witness(
    parameters,
    data_top,
    data_turn,
    data_bottom,
    data_connection_points,
)
    @info "Checking top part"
    check_proof_witness_part(data_top, false)

    @info "Checking turn part"
    check_proof_witness_part(data_turn, true)

    @info "Checking bottom part"
    check_proof_witness_part(data_bottom, true)

    @info "Checking connection between top and turn"
    check_proof_witness_top_turn(parameters, data_top, data_turn, data_connection_points)

    @info "Checking connection between turn and bottom"
    check_proof_witness_turn_bottom(
        parameters,
        data_turn,
        data_bottom,
        data_connection_points,
    )

    return true
end
