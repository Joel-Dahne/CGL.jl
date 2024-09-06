function branch_points_dataframe_fix_epsilon(
    points::Vector{SVector{4,Arb}},
    μ₀s::Vector{Arb},
    κ₀s::Vector{Arb},
    ϵ₀s::Vector{Arb},
    ξ₁s::Vector{Arb},
)
    df = DataFrame(
        μ = getindex.(points, 1),
        γ = Acb.(getindex.(points, 2), getindex.(points, 3)),
        κ = getindex.(points, 4),
        ϵ = ϵ₀s,
        μ₀ = μ₀s,
        κ₀ = κ₀s,
        ξ₁ = ξ₁s,
    )

    return df
end

function branch_points_dataframe_fix_kappa(
    points::Vector{SVector{4,Arb}},
    μ₀s::Vector{Arb},
    κ₀s::Vector{Arb},
    ϵ₀s::Vector{Arb},
    ξ₁s::Vector{Arb},
)
    df = DataFrame(
        μ = getindex.(points, 1),
        γ = Acb.(getindex.(points, 2), getindex.(points, 3)),
        κ = κ₀s,
        ϵ = getindex.(points, 4),
        μ₀ = μ₀s,
        ϵ₀ = ϵ₀s,
        ξ₁ = ξ₁s,
    )

    return df
end

function branch_existence_dataframe_fix_epsilon(
    ϵs::Vector{NTuple{2,Arf}},
    uniqs::Vector{SVector{4,Arb}},
    exists::Vector{SVector{4,Arb}},
    approxs::Vector{SVector{4,Arb}},
)
    df = DataFrame(
        ϵ_lower = getindex.(ϵs, 1),
        ϵ_upper = getindex.(ϵs, 2),
        μ_uniq = getindex.(uniqs, 1),
        γ_uniq = Acb.(getindex.(uniqs, 2), getindex.(uniqs, 3)),
        κ_uniq = getindex.(uniqs, 4),
        μ_exists = getindex.(exists, 1),
        γ_exists = Acb.(getindex.(exists, 2), getindex.(exists, 3)),
        κ_exists = getindex.(exists, 4),
        μ_approx = getindex.(approxs, 1),
        γ_approx = Acb.(getindex.(approxs, 2), getindex.(approxs, 3)),
        κ_approx = getindex.(approxs, 4),
    )

    return df
end

function branch_existence_dataframe_fix_kappa(
    κs::Vector{NTuple{2,Arf}},
    uniqs::Vector{SVector{4,Arb}},
    exists::Vector{SVector{4,Arb}},
    approxs::Vector{SVector{4,Arb}},
)
    df = DataFrame(
        κ_lower = getindex.(κs, 1),
        κ_upper = getindex.(κs, 2),
        μ_uniq = getindex.(uniqs, 1),
        γ_uniq = Acb.(getindex.(uniqs, 2), getindex.(uniqs, 3)),
        ϵ_uniq = getindex.(uniqs, 4),
        μ_exists = getindex.(exists, 1),
        γ_exists = Acb.(getindex.(exists, 2), getindex.(exists, 3)),
        ϵ_exists = getindex.(exists, 4),
        μ_approx = getindex.(approxs, 1),
        γ_approx = Acb.(getindex.(approxs, 2), getindex.(approxs, 3)),
        ϵ_approx = getindex.(approxs, 4),
    )

    return df
end

function branch_continuation_dataframe_fix_epsilon(
    left_continuation::Union{Vector{Bool},BitVector},
    ϵs::Vector{NTuple{2,Arf}},
    uniqs::Vector{SVector{4,Arb}},
    exists::Vector{SVector{4,Arb}},
    approxs::Vector{SVector{4,Arb}},
)
    df = DataFrame(
        left_continuation = convert(Vector{Bool}, left_continuation),
        ϵ_lower = getindex.(ϵs, 1),
        ϵ_upper = getindex.(ϵs, 2),
        μ_uniq = getindex.(uniqs, 1),
        γ_uniq = Acb.(getindex.(uniqs, 2), getindex.(uniqs, 3)),
        κ_uniq = getindex.(uniqs, 4),
        μ_exists = getindex.(exists, 1),
        γ_exists = Acb.(getindex.(exists, 2), getindex.(exists, 3)),
        κ_exists = getindex.(exists, 4),
        μ_approx = getindex.(approxs, 1),
        γ_approx = Acb.(getindex.(approxs, 2), getindex.(approxs, 3)),
        κ_approx = getindex.(approxs, 4),
    )

    return df
end

function branch_continuation_dataframe_fix_kappa(
    left_continuation::Union{Vector{Bool},BitVector},
    κs::Vector{NTuple{2,Arf}},
    uniqs::Vector{SVector{4,Arb}},
    exists::Vector{SVector{4,Arb}},
    approxs::Vector{SVector{4,Arb}},
)
    df = DataFrame(
        left_continuation = convert(Vector{Bool}, left_continuation),
        κ_lower = getindex.(κs, 1),
        κ_upper = getindex.(κs, 2),
        μ_uniq = getindex.(uniqs, 1),
        γ_uniq = Acb.(getindex.(uniqs, 2), getindex.(uniqs, 3)),
        ϵ_uniq = getindex.(uniqs, 4),
        μ_exists = getindex.(exists, 1),
        γ_exists = Acb.(getindex.(exists, 2), getindex.(exists, 3)),
        ϵ_exists = getindex.(exists, 4),
        μ_approx = getindex.(approxs, 1),
        γ_approx = Acb.(getindex.(approxs, 2), getindex.(approxs, 3)),
        ϵ_approx = getindex.(approxs, 4),
    )

    return df
end

function branch_critical_points_dataframe(
    μs::Vector{Arb},
    γs::Vector{Acb},
    κs::Vector{Arb},
    ϵs::Vector{Arb},
    ξ₁s::Vector{Arb},
    num_critical_points::Vector{Union{Int,Missing}},
)
    df = DataFrame(
        μ = μs,
        γ = γs,
        κ = κs,
        ϵ = ϵs,
        ξ₁ = ξ₁s,
        num_critical_points = num_critical_points,
    )

    return df
end

function write_branch_csv(filename, data::DataFrame; compress = endswith(filename, ".gz"))
    data_dump = DataFrame()

    for col_name in names(data)
        col = data[!, col_name]
        if eltype(col) <: Union{Arf,Arb}
            insertcols!(data_dump, col_name * "_dump" => Arblib.dump_string.(col))
        elseif eltype(col) <: Acb
            insertcols!(
                data_dump,
                col_name * "_dump_real" => Arblib.dump_string.(real.(col)),
                col_name * "_dump_imag" => Arblib.dump_string.(imag.(col)),
            )
        else
            insertcols!(data_dump, col_name => col)
        end
    end

    CSV.write(filename, data_dump; compress)
end

function read_branch_points_csv(filename)
    types = [String, String, String, String, String, String, String, String]

    data_dump = CSV.read(filename, DataFrame; types)

    data = DataFrame()

    data.μ = Arblib.load_string.(Arb, data_dump.μ_dump)
    data.γ =
        Acb.(
            Arblib.load_string.(Arb, data_dump.γ_dump_real),
            Arblib.load_string.(Arb, data_dump.γ_dump_imag),
        )
    data.κ = Arblib.load_string.(Arb, data_dump.κ_dump)
    data.ϵ = Arblib.load_string.(Arb, data_dump.ϵ_dump)

    data.μ₀ = Arblib.load_string.(Arb, data_dump.μ₀_dump)
    if "κ₀_dump" in names(data_dump)
        data.κ₀ = Arblib.load_string.(Arb, data_dump.κ₀_dump)
    else
        data.ϵ₀ = Arblib.load_string.(Arb, data_dump.ϵ₀_dump)
    end

    data.ξ₁ = Arblib.load_string.(Arb, data_dump.ξ₁_dump)

    return data
end

function read_branch_csv_helper(data_dump::DataFrame)
    is_turn = "κ_lower_dump" in names(data_dump)

    data = DataFrame()

    if !is_turn
        data.ϵ_lower = Arblib.load_string.(Arf, data_dump.ϵ_lower_dump)
        data.ϵ_upper = Arblib.load_string.(Arf, data_dump.ϵ_upper_dump)
    else
        data.κ_lower = Arblib.load_string.(Arf, data_dump.κ_lower_dump)
        data.κ_upper = Arblib.load_string.(Arf, data_dump.κ_upper_dump)
    end

    data.μ_uniq = Arblib.load_string.(Arb, data_dump.μ_uniq_dump)
    data.γ_uniq =
        Acb.(
            Arblib.load_string.(Arb, data_dump.γ_uniq_dump_real),
            Arblib.load_string.(Arb, data_dump.γ_uniq_dump_imag),
        )
    if !is_turn
        data.κ_uniq = Arblib.load_string.(Arb, data_dump.κ_uniq_dump)
    else
        data.ϵ_uniq = Arblib.load_string.(Arb, data_dump.ϵ_uniq_dump)
    end

    data.μ_exists = Arblib.load_string.(Arb, data_dump.μ_exists_dump)
    data.γ_exists =
        Acb.(
            Arblib.load_string.(Arb, data_dump.γ_exists_dump_real),
            Arblib.load_string.(Arb, data_dump.γ_exists_dump_imag),
        )
    if !is_turn
        data.κ_exists = Arblib.load_string.(Arb, data_dump.κ_exists_dump)
    else
        data.ϵ_exists = Arblib.load_string.(Arb, data_dump.ϵ_exists_dump)
    end

    data.μ_approx = Arblib.load_string.(Arb, data_dump.μ_approx_dump)
    data.γ_approx =
        Acb.(
            Arblib.load_string.(Arb, data_dump.γ_approx_dump_real),
            Arblib.load_string.(Arb, data_dump.γ_approx_dump_imag),
        )
    if !is_turn
        data.κ_approx = Arblib.load_string.(Arb, data_dump.κ_approx_dump)
    else
        data.ϵ_approx = Arblib.load_string.(Arb, data_dump.ϵ_approx_dump)
    end

    return data
end


function read_branch_existence_csv(filename)
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
        String,
        String,
        String,
        String,
        String,
    ]

    data_dump = CSV.read(filename, DataFrame; types)

    return read_branch_csv_helper(data_dump)
end

function read_branch_continuation_csv(filename)
    types = [
        Bool,
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
        String,
        String,
        String,
        String,
        String,
    ]

    data_dump = CSV.read(filename, DataFrame; types)

    data = read_branch_csv_helper(data_dump)

    insertcols!(data, 1, :left_continuation => data_dump.left_continuation)

    return data
end

function read_branch_critical_points_csv(filename)
    types = [String, String, String, String, String, String, Int]

    data_dump = CSV.read(filename, DataFrame; types)

    data = DataFrame()

    data.μ = Arblib.load_string.(Arb, data_dump.μ_dump)
    data.γ =
        Acb.(
            Arblib.load_string.(Arb, data_dump.γ_dump_real),
            Arblib.load_string.(Arb, data_dump.γ_dump_imag),
        )
    data.κ = Arblib.load_string.(Arb, data_dump.κ_dump)
    data.ϵ = Arblib.load_string.(Arb, data_dump.ϵ_dump)

    data.ξ₁ = Arblib.load_string.(Arb, data_dump.ξ₁_dump)

    data.num_critical_points = data_dump.num_critical_points

    return data
end


function write_parameters(filename, ξ₁::Arb, λ::CGLParams{Arb}; kwargs...)
    parameters_raw = DataFrame(
        d = [λ.d],
        ω_dump = [Arblib.dump_string(λ.ω)],
        σ_dump = [Arblib.dump_string(λ.σ)],
        δ_dump = [Arblib.dump_string(λ.δ)],
        ξ₁_dump = [Arblib.dump_string(ξ₁)],
    )
    for (key, value) in kwargs
        insertcols!(parameters_raw, key => [value])
    end

    CSV.write(filename, parameters_raw)

    return parameters_raw
end

function read_parameters(filename)
    types = Dict(
        :d => Int,
        :ω_dump => String,
        :σ_dump => String,
        :δ_dump => String,
        :ξ₁_dump => String,
        :ξ₁_strategy => Symbol,
        :ξ₁_strategy_value => String,
    )
    parameters =
        CSV.read(filename, DataFrame, types = (_, name) -> get(types, name, nothing))

    # Extract special parameters
    ξ₁ = Arblib.load_string(Arb, only(parameters.ξ₁_dump))
    λ = CGLParams(
        only(parameters.d),
        Arblib.load_string(Arb, only(parameters.ω_dump)),
        Arblib.load_string(Arb, only(parameters.σ_dump)),
        Arblib.load_string(Arb, only(parameters.δ_dump)),
    )

    select!(parameters, Not([:d, :ω_dump, :σ_dump, :δ_dump, :ξ₁_dump]))

    if isempty(parameters)
        return (ξ₁ = ξ₁, λ = λ)
    else
        return (ξ₁ = ξ₁, λ = λ, NamedTuple(parameters[1, :])...)
    end
end

function locate_most_recent(
    data_type::AbstractString,
    j::Integer,
    d::Integer,
    part::AbstractString;
    return_filename::Bool = false,
    verbose::Bool = false,
)
    base_directory = relpath(
        joinpath(
            dirname(pathof(@__MODULE__)),
            "../Dardel/output/branch_$(data_type)_j=$(j)_d=$(d)_part=$(part)",
        ),
    )

    verbose && @info "Searching for most recent data in" base_directory

    if !isdir(base_directory)
        verbose && @warn "Directory doesn't exist"
        return nothing
    end

    directory = maximum(readdir(base_directory, join = true))

    if return_filename
        return joinpath(directory, "branch_$(data_type)_j=$(j)_d=$(d)_part=$(part).csv.gz")
    else
        return directory
    end
end
