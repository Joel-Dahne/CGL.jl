function _construct_proof_witness_load_data(
    j::Integer,
    d::Integer,
    directory_top::Union{AbstractString,Nothing},
    directory_turn::Union{AbstractString,Nothing},
    directory_bottom::Union{AbstractString,Nothing},
)
    @info "Reading data"

    if !isnothing(directory_top)
        @info "Reading data for top"
        data_top = read_branch_continuation_csv(
            joinpath(directory_top, "branch_continuation_d=$(d)_j=$(j)_part=top.csv.gz"),
        )
        parameters_top = read_parameters(joinpath(directory_top, "parameters.csv"))

    else
        @info "No data for top"
        data_top = nothing
        parameters_top = (
            ξ₁ = nothing,
            λ = nothing,
            runtime_existence = NaN,
            runtime_continuation = NaN,
            commit_hash = "",
        )
    end

    if !isnothing(directory_turn)
        @info "Reading data for turn"
        data_turn = read_branch_continuation_csv(
            joinpath(directory_turn, "branch_continuation_d=$(d)_j=$(j)_part=turn.csv.gz"),
        )
        parameters_turn = read_parameters(joinpath(directory_turn, "parameters.csv"))

    else
        @info "No data for turn"
        data_turn = nothing
        parameters_turn = (
            ξ₁ = nothing,
            λ = nothing,
            runtime_existence = NaN,
            runtime_continuation = NaN,
            commit_hash = "",
        )
    end

    if !isnothing(directory_bottom)
        @info "Reading data for bottom"
        data_bottom = read_branch_continuation_csv(
            joinpath(
                directory_bottom,
                "branch_continuation_d=$(d)_j=$(j)_part=bottom.csv.gz",
            ),
        )
        parameters_bottom = read_parameters(joinpath(directory_bottom, "parameters.csv"))
    else
        @info "No data for bottom"
        data_bottom = nothing
        parameters_bottom = (
            ξ₁ = nothing,
            λ = nothing,
            runtime_existence = NaN,
            runtime_continuation = NaN,
            commit_hash = "",
        )
    end

    isnothing(data_top) &&
        isnothing(data_turn) &&
        isnothing(data_bottom) &&
        error("no data read")

    ξ₁s = filter(!isnothing, [parameters_top.ξ₁, parameters_turn.ξ₁, parameters_bottom.ξ₁])
    @assert allequal(ξ₁s)
    ξ₁ = ξ₁s[1]

    λs = filter(!isnothing, [parameters_top.λ, parameters_turn.λ, parameters_bottom.λ])
    @assert allequal(λs)
    λ = λs[1]

    parameters = (;
        ξ₁,
        λ,
        runtime_existence_top = parameters_top.runtime_existence,
        runtime_existence_turn = parameters_turn.runtime_existence,
        runtime_existence_bottom = parameters_bottom.runtime_existence,
        runtime_continuation_top = parameters_top.runtime_continuation,
        runtime_continuation_turn = parameters_turn.runtime_continuation,
        runtime_continuation_bottom = parameters_bottom.runtime_continuation,
        commit_hash_top = parameters_top.commit_hash,
        commit_hash_turn = parameters_turn.commit_hash,
        commit_hash_bottom = parameters_bottom.commit_hash,
    )

    return parameters, data_top, data_turn, data_bottom
end

function _construct_proof_witness_load_data_critical_points(
    j::Integer,
    d::Integer,
    directory_top::Union{AbstractString,Nothing},
    directory_turn::Union{AbstractString,Nothing},
    directory_bottom::Union{AbstractString,Nothing},
)
    @info "Reading data about critical points"

    if !isnothing(directory_top)
        @info "Reading data for top"
        data_top = read_branch_critical_points_csv(
            joinpath(directory_top, "branch_critical_points_d=$(d)_j=$(j)_part=top.csv.gz"),
        )
        parameters_top = read_parameters(joinpath(directory_top, "parameters.csv"))

        allequal(data_top.ξ₁) || error("top doesn't have same ξ₁ for all segments")
        parameters_top.use_midpoint && error("top uses midpoint for critical points")
    else
        @info "No data for top"
        data_top = nothing
        parameters_top =
            (ξ₁ = nothing, λ = nothing, runtime_critical_points = NaN, commit_hash = "")
    end

    if !isnothing(directory_turn)
        @info "Reading data for turn"
        data_turn = read_branch_critical_points_csv(
            joinpath(
                directory_turn,
                "branch_critical_points_d=$(d)_j=$(j)_part=turn.csv.gz",
            ),
        )
        parameters_turn = read_parameters(joinpath(directory_turn, "parameters.csv"))

        allequal(data_turn.ξ₁) || error("turn doesn't have same ξ₁ for all segments")
        parameters_turn.use_midpoint && error("turn uses midpoint for critical points")
    else
        @info "No data for turn"
        data_turn = nothing
        parameters_turn =
            (ξ₁ = nothing, λ = nothing, runtime_critical_points = NaN, commit_hash = "")
    end

    if !isnothing(directory_bottom)
        @info "Reading data for bottom"
        data_bottom = read_branch_critical_points_csv(
            joinpath(
                directory_bottom,
                "branch_critical_points_d=$(d)_j=$(j)_part=bottom.csv.gz",
            ),
        )
        parameters_bottom = read_parameters(joinpath(directory_bottom, "parameters.csv"))

        allequal(data_bottom.ξ₁) || error("bottom doesn't have same ξ₁ for all segments")
        parameters_bottom.use_midpoint && error("bottom uses midpoint for critical points")
    else
        @info "No data for bottom"
        data_bottom = nothing
        parameters_bottom =
            (ξ₁ = nothing, λ = nothing, runtime_critical_points = NaN, commit_hash = "")
    end

    if isnothing(data_top) && isnothing(data_turn) && isnothing(data_bottom)
        ξ₁ = nothing
        λ = nothing
    else
        ξ₁s = filter(
            !isnothing,
            [parameters_top.ξ₁, parameters_turn.ξ₁, parameters_bottom.ξ₁],
        )
        @assert allequal(ξ₁s)
        ξ₁ = ξ₁s[1]

        λs = filter(!isnothing, [parameters_top.λ, parameters_turn.λ, parameters_bottom.λ])
        @assert allequal(λs)
        λ = λs[1]
    end

    parameters = (;
        ξ₁,
        λ,
        runtime_critical_points_top = parameters_top.runtime_critical_points,
        runtime_critical_points_turn = parameters_turn.runtime_critical_points,
        runtime_critical_points_bottom = parameters_bottom.runtime_critical_points,
        commit_hash_critical_points_top = parameters_top.commit_hash,
        commit_hash_critical_points_turn = parameters_turn.commit_hash,
        commit_hash_critical_points_bottom = parameters_bottom.commit_hash,
    )

    return parameters, data_top, data_turn, data_bottom
end

function construct_proof_witness(j::Integer, d::Integer)
    directory_top = locate_most_recent("continuation", j, d, "top")
    directory_turn = locate_most_recent("continuation", j, d, "turn")
    directory_bottom = locate_most_recent("continuation", j, d, "bottom")

    directory_critical_points_top = locate_most_recent("critical_points", j, d, "top")
    directory_critical_points_turn = locate_most_recent("critical_points", j, d, "turn")
    directory_critical_points_bottom = locate_most_recent("critical_points", j, d, "bottom")

    return CGL.construct_proof_witness(
        j,
        d,
        directory_top,
        directory_turn,
        directory_bottom,
        directory_critical_points_top,
        directory_critical_points_turn,
        directory_critical_points_bottom,
    )
end

function construct_proof_witness(
    j::Integer,
    d::Integer,
    directory_top::Union{AbstractString,Nothing},
    directory_turn::Union{AbstractString,Nothing},
    directory_bottom::Union{AbstractString,Nothing},
    directory_critical_points_top::Union{AbstractString,Nothing} = nothing,
    directory_critical_points_turn::Union{AbstractString,Nothing} = nothing,
    directory_critical_points_bottom::Union{AbstractString,Nothing} = nothing,
)
    parameters, data_top, data_turn, data_bottom = _construct_proof_witness_load_data(
        j,
        d,
        directory_top,
        directory_turn,
        directory_bottom,
    )

    parameters_critical_points,
    data_critical_points_top,
    data_critical_points_turn,
    data_critical_points_bottom = _construct_proof_witness_load_data_critical_points(
        j,
        d,
        directory_critical_points_top,
        directory_critical_points_turn,
        directory_critical_points_bottom,
    )

    parameters = (;
        parameters...,
        parameters_critical_points.runtime_critical_points_top,
        parameters_critical_points.runtime_critical_points_turn,
        parameters_critical_points.runtime_critical_points_bottom,
        parameters_critical_points.commit_hash_critical_points_top,
        parameters_critical_points.commit_hash_critical_points_turn,
        parameters_critical_points.commit_hash_critical_points_bottom,
    )

    ## Sanity check data

    # Continuation should have succeeded for all of them
    isnothing(data_top) ||
        all(data_top.left_continuation) ||
        error("left continuation failed for top")
    isnothing(data_turn) ||
        all(data_turn.left_continuation) ||
        error("left continuation failed for turn")
    isnothing(data_bottom) ||
        all(data_bottom.left_continuation) ||
        error("left continuation failed for bottom")

    if parameters.λ.d == 1
        isnothing(data_top) && error("expected top part for d = 1")
        iszero(data_top.ϵ_lower[1]) || error("expected ϵ to start at 0 for top for d = 1")
    end

    # Continuation data and critical points data should have same ξ₁ and λ
    isequal(parameters.ξ₁, parameters_critical_points.ξ₁)
    isequal(parameters.λ, parameters_critical_points.λ)

    # Continuation data and critical points data should be the same.

    # Top
    if !isnothing(data_top) && !isnothing(data_critical_points_top)
        all(isequal.(data_top.μ_exists, data_critical_points_top.μ)) ||
            error("μ not same in top for continuation and critical points")
        all(isequal.(data_top.γ_exists, data_critical_points_top.γ)) ||
            error("γ not same in top for continuation and critical points")
        all(isequal.(data_top.κ_exists, data_critical_points_top.κ)) ||
            error("κ not same in top for continuation and critical points")
        all(
            isequal.(
                Arb.(tuple.(data_top.ϵ_lower, data_top.ϵ_upper)),
                data_critical_points_top.ϵ,
            ),
        ) || error("ϵ not same in top for continuation and critical points")
    end

    # Turn
    if !isnothing(data_turn) && !isnothing(data_critical_points_turn)
        all(isequal.(data_turn.μ_exists, data_critical_points_turn.μ)) ||
            error("μ not same in turn for continuation and critical points")
        all(isequal.(data_turn.γ_exists, data_critical_points_turn.γ)) ||
            error("γ not same in turn for continuation and critical points")
        all(
            isequal.(
                Arb.(tuple.(data_turn.κ_lower, data_turn.κ_upper)),
                data_critical_points_turn.κ,
            ),
        ) || error("κ not same in turn for continuation and critical points")
        all(isequal.(data_turn.ϵ_exists, data_critical_points_turn.ϵ)) ||
            error("ϵ not same in turn for continuation and critical points")
    end

    # Bottom
    if !isnothing(data_bottom) && !isnothing(data_critical_points_bottom)
        all(isequal.(data_bottom.μ_exists, data_critical_points_bottom.μ)) ||
            error("μ not same in bottom for continuation and critical points")
        all(isequal.(data_bottom.γ_exists, data_critical_points_bottom.γ)) ||
            error("γ not same in bottom for continuation and critical points")
        all(isequal.(data_bottom.κ_exists, data_critical_points_bottom.κ)) ||
            error("κ not same in bottom for continuation and critical points")
        all(
            isequal.(
                Arb.(tuple.(data_bottom.ϵ_lower, data_bottom.ϵ_upper)),
                data_critical_points_bottom.ϵ,
            ),
        ) || error("ϵ not same in bottom for continuation and critical points")
    end

    commit_hashes = filter(
        !isempty,
        [
            parameters.commit_hash_top,
            parameters.commit_hash_turn,
            parameters.commit_hash_bottom,
            parameters.commit_hash_critical_points_top,
            parameters.commit_hash_critical_points_turn,
            parameters.commit_hash_critical_points_bottom,
        ],
    )

    if allequal(commit_hashes)
        commit_hash = commit_hashes[1]
    else
        commit_hash = ""
        @warn "All data not from same commit" commit_hashes
    end
    parameters = (; parameters..., commit_hash)

    ## Construct output

    if !isnothing(data_top)
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
    end

    if !isnothing(data_turn)
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
    end

    if !isnothing(data_bottom)
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
    end

    # Add data about critical points
    if !isnothing(data_top) && !isnothing(data_critical_points_top)
        data_top[!, :num_critical_points] = data_critical_points_top.num_critical_points
    elseif !isnothing(data_top)
        data_top[!, :num_critical_points] .= missing
    end

    if !isnothing(data_turn) && !isnothing(data_critical_points_turn)
        data_turn[!, :num_critical_points] = data_critical_points_turn.num_critical_points
    elseif !isnothing(data_turn)
        data_turn[!, :num_critical_points] .= missing
    end

    if !isnothing(data_bottom) && !isnothing(data_critical_points_bottom)
        data_bottom[!, :num_critical_points] =
            data_critical_points_bottom.num_critical_points
    elseif !isnothing(data_bottom)
        data_bottom[!, :num_critical_points] .= missing
    end

    @info "Checking top part"
    check_proof_witness_part(data_top, false)

    @info "Checking turn part"
    check_proof_witness_part(data_turn, true)

    @info "Checking bottom part"
    check_proof_witness_part(data_bottom, true)

    # Check that the number of critical points are the same along the
    # full branch
    num_critical_points_full = Int[]
    if !isnothing(data_top)
        num_critical_points_full =
            vcat(num_critical_points_full, data_top.num_critical_points)
    end
    if !isnothing(data_turn)
        num_critical_points_full =
            vcat(num_critical_points_full, data_turn.num_critical_points)
    end
    if !isnothing(data_bottom)
        num_critical_points_full =
            vcat(num_critical_points_full, data_bottom.num_critical_points)
    end

    allequal(skipmissing(num_critical_points_full)) ||
        error("not all segments have the same number of critical points")

    # Find points connecting the segments into one

    data_connection_points = DataFrame(
        ϵ = Arb[],
        μ_uniq = Arb[],
        γ_uniq = Acb[],
        κ_uniq = Arb[],
        μ_exists = Arb[],
        γ_exists = Acb[],
        κ_exists = Arb[],
    )

    if !isnothing(data_top) && !isnothing(data_turn)
        @info "Generating connection point for top and turn"
        ϵ_top, exists_top, uniq_top =
            construct_proof_witness_connect(parameters, data_top, data_turn, true)

        push!(
            data_connection_points,
            (
                ϵ_top,
                uniq_top[1],
                Acb(uniq_top[2], uniq_top[3]),
                uniq_top[4],
                exists_top[1],
                Acb(exists_top[2], exists_top[3]),
                exists_top[4],
            ),
        )
    end

    if !isnothing(data_turn) && !isnothing(data_bottom)
        @info "Generating connection point for turn and bottom"
        ϵ_bottom, exists_bottom, uniq_bottom =
            construct_proof_witness_connect(parameters, data_bottom, data_turn, false)
        push!(
            data_connection_points,
            (
                ϵ_bottom,
                uniq_bottom[1],
                Acb(uniq_bottom[2], uniq_bottom[3]),
                uniq_bottom[4],
                exists_bottom[1],
                Acb(exists_bottom[2], exists_bottom[3]),
                exists_bottom[4],
            ),
        )
    end

    return parameters, data_top, data_turn, data_bottom, data_connection_points
end

function construct_proof_witness_connect(
    parameters,
    data_top_or_bottom,
    data_turn,
    is_top;
    verbose = false,
)
    # For the top we take the connection point to be in the middle (in
    # ϵ) of the last interval. For the bottom we take the first
    # interval instead.
    i = ifelse(is_top, nrow(data_top_or_bottom), 1)

    ϵ = midpoint(Arb, Arb((data_top_or_bottom.ϵ_lower[i], data_top_or_bottom.ϵ_upper[i])))
    @assert data_top_or_bottom.ϵ_lower[i] < ϵ < data_top_or_bottom.ϵ_upper[i]

    # Compute enclosure for the connection point
    (; ξ₁, λ) = parameters
    μ, γ, κ = refine_approximation_fix_epsilon(
        data_top_or_bottom.μ_exists[i],
        data_top_or_bottom.κ_exists[i],
        ϵ,
        ξ₁,
        λ,
    )

    exists, uniq = G_solve_fix_epsilon(
        μ,
        real(γ),
        imag(γ),
        κ,
        ϵ,
        ξ₁,
        λ,
        return_uniqueness = Val{true}(),
    )

    # Check that the connection point is contained in the interior of
    # the interval for the top/bottom part
    exists_top_or_bottom = SVector(
        data_top_or_bottom.μ_exists[i],
        real(data_top_or_bottom.γ_exists[i]),
        imag(data_top_or_bottom.γ_exists[i]),
        data_top_or_bottom.κ_exists[i],
    )

    all(Arblib.contains_interior.(exists_top_or_bottom, exists)) ||
        error("enclosure for picked point not contained in top/bottom")

    # Check that the connection point is contained in some of interval
    # for the turning part

    # Find interval for which we want to check containment
    j = searchsortedfirst(data_turn.κ_lower, midpoint(exists[4]), rev = true)

    # If this fails its because the κ enclosure lies on the
    # intersection of two intervals. We currently don't handle this
    # case as it has not occurred.
    data_turn.κ_lower[j] < exists[4] < data_turn.κ_upper[j] ||
        error("κ for picked point not contained in an interval for turn")

    Arblib.contains_interior(data_turn.μ_exists[j], exists[1]) ||
        error("μ for picked point not contained in interval for turn")
    Arblib.contains_interior(real(data_turn.γ_exists[j]), exists[2]) ||
        error("real(γ) for picked point not contained in interval for turn")
    Arblib.contains_interior(imag(data_turn.γ_exists[j]), exists[3]) ||
        error("imag(γ) for picked point not contained in interval for turn")
    Arblib.contains_interior(data_turn.ϵ_exists[j], ϵ) ||
        error("ϵ for picked point not contained in interval for turn")

    return ϵ, exists, uniq
end


function write_proof_witness(
    directory,
    parameters,
    data_top::Union{DataFrame,Nothing},
    data_turn::Union{DataFrame,Nothing},
    data_bottom::Union{DataFrame,Nothing},
    data_connection_points::DataFrame,
)
    mkpath(directory)

    write_parameters(
        joinpath(directory, "parameters.csv"),
        parameters.ξ₁,
        parameters.λ;
        Base.structdiff(parameters, NamedTuple{(:ξ₁, :λ)})...,
    )

    if !isnothing(data_top)
        write_branch_csv(joinpath(directory, "top.csv.gz"), data_top)
    end
    if !isnothing(data_turn)
        write_branch_csv(joinpath(directory, "turn.csv.gz"), data_turn)
    end
    if !isnothing(data_bottom)
        write_branch_csv(joinpath(directory, "bottom.csv.gz"), data_bottom)
    end
    write_branch_csv(
        joinpath(directory, "connection_points.csv.gz"),
        data_connection_points,
    )

    return
end

function read_proof_witness(directory::AbstractString)
    parameters = read_parameters(joinpath(directory, "parameters.csv"))

    types = [
        String,
        String,
        String,
        String,
        String,
        String,
        String,
        String,
        String,
        String,
        Int,
    ]

    data_top = if ispath(joinpath(directory, "top.csv.gz"))
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

        data_top.num_critical_points = data_top_dump.num_critical_points

        data_top
    else
        nothing
    end

    data_turn = if ispath(joinpath(directory, "turn.csv.gz"))
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

        data_turn.num_critical_points = data_turn_dump.num_critical_points

        data_turn
    else
        nothing
    end

    data_bottom = if ispath(joinpath(directory, "bottom.csv.gz"))
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

        data_bottom.num_critical_points = data_bottom_dump.num_critical_points

        data_bottom
    else
        nothing
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

function check_proof_witness_part(data::Nothing, reversed)
    @info "No data to check"
    return true
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

    # Check that the number of critical points are the same for all
    # segments
    allequal(skipmissing(data.num_critical_points)) ||
        error("not all segments have the same number of critical points")

    return
end

function check_proof_witness_top_turn(
    parameters,
    data_top,
    data_turn,
    data_connection_point::DataFrameRow,
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

    # Check that connection point is contained in both branches

    point = SVector(
        data_connection_point.μ_exists,
        real(data_connection_point.γ_exists),
        imag(data_connection_point.γ_exists),
        data_connection_point.κ_exists,
        data_connection_point.ϵ,
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
    data_connection_point::DataFrameRow,
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

    # Check that the connection point is contained in both branches

    point = SVector(
        data_connection_point.μ_exists,
        real(data_connection_point.γ_exists),
        imag(data_connection_point.γ_exists),
        data_connection_point.κ_exists,
        data_connection_point.ϵ,
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

    if !isnothing(data_top) && !isnothing(data_turn)
        @info "Connection between top and turn OK"
        check_proof_witness_top_turn(
            parameters,
            data_top,
            data_turn,
            data_connection_points[1, :],
        )
    else
        @info "No connection between top and turn to check"
    end

    if !isnothing(data_turn) && !isnothing(data_bottom)
        @info "Connection between turn and bottom OK"
        check_proof_witness_turn_bottom(
            parameters,
            data_turn,
            data_bottom,
            data_connection_points[ifelse(isnothing(data_top), 1, 2), :],
        )
    else
        @info "No connection between turn and bottom to check"
    end

    return true
end
