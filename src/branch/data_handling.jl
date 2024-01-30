function branch_points_dataframe(
    λs::Vector{CGLParams{Arb}},
    μs::Vector{Arb},
    γs::Vector{Acb},
    κs::Vector{Arb},
    μs_approx::Vector{Arb},
    κs_approx::Vector{Arb},
    ξ₁s::Vector{Arb},
)
    df = DataFrame(
        ϵ = getproperty.(λs, :ϵ),
        μ = μs,
        γ = γs,
        κ = κs,
        μ_approx = μs_approx,
        κ_approx = κs_approx,
        ξ₁ = ξ₁s,
    )

    return df
end

function branch_existence_dataframe(
    ϵs::Vector{NTuple{2,Arf}},
    μs_uniq::Vector{Arb},
    γs_uniq::Vector{Acb},
    κs_uniq::Vector{Arb},
    μs_exists::Vector{Arb},
    γs_exists::Vector{Acb},
    κs_exists::Vector{Arb},
    μs_approx::Vector{Arb},
    γs_approx::Vector{Acb},
    κs_approx::Vector{Arb},
    ξ₁s::Vector{Arb},
)
    df = DataFrame(
        ϵ_lower = getindex.(ϵs, 1),
        ϵ_upper = getindex.(ϵs, 2),
        μ_uniq = μs_uniq,
        γ_uniq = γs_uniq,
        κ_uniq = κs_uniq,
        μ_exists = μs_exists,
        γ_exists = γs_exists,
        κ_exists = κs_exists,
        μ_approx = μs_approx,
        γ_approx = γs_approx,
        κ_approx = κs_approx,
        ξ₁ = ξ₁s,
    )

    return df
end

function branch_continuation_dataframe(
    left_continuation::Union{Vector{Bool},BitVector},
    ϵs::Vector{NTuple{2,Arf}},
    μs_uniq::Vector{Arb},
    γs_uniq::Vector{Acb},
    κs_uniq::Vector{Arb},
    μs_exists::Vector{Arb},
    γs_exists::Vector{Acb},
    κs_exists::Vector{Arb},
    μs_approx::Vector{Arb},
    γs_approx::Vector{Acb},
    κs_approx::Vector{Arb},
    ξ₁s::Vector{Arb},
)
    df = DataFrame(
        left_continuation = convert(Vector{Bool}, left_continuation),
        ϵ_lower = getindex.(ϵs, 1),
        ϵ_upper = getindex.(ϵs, 2),
        μ_uniq = μs_uniq,
        γ_uniq = γs_uniq,
        κ_uniq = κs_uniq,
        μ_exists = μs_exists,
        γ_exists = γs_exists,
        κ_exists = κs_exists,
        μ_approx = μs_approx,
        γ_approx = γs_approx,
        κ_approx = κs_approx,
        ξ₁ = ξ₁s,
    )

    return df
end

function write_branch_points_csv(filename, data::DataFrame)
    data_dump = DataFrame()

    for col_name in names(data)
        col = data[!, col_name]
        if eltype(col) <: Arb
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

    CSV.write(filename, data_dump)
end

function write_branch_existence_csv(filename, data::DataFrame)
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

    CSV.write(filename, data_dump)
end

function write_branch_continuation_csv(filename, data::DataFrame)
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

    CSV.write(filename, data_dump)
end

function read_branch_points_csv(filename)
    types = [String, String, String, String, String, String, String, String]

    data_dump = CSV.read(filename, DataFrame; types)

    data = DataFrame()

    data.ϵ = Arblib.load_string!.(zeros(Arb, size(data_dump, 1)), data_dump.ϵ_dump)

    data.μ = Arblib.load_string!.(zeros(Arb, size(data_dump, 1)), data_dump.μ_dump)
    data.γ =
        Acb.(
            Arblib.load_string!.(zeros(Arb, size(data_dump, 1)), data_dump.γ_dump_real),
            Arblib.load_string!.(zeros(Arb, size(data_dump, 1)), data_dump.γ_dump_imag),
        )
    data.κ = Arblib.load_string!.(zeros(Arb, size(data_dump, 1)), data_dump.κ_dump)

    data.μ_approx =
        Arblib.load_string!.(zeros(Arb, size(data_dump, 1)), data_dump.μ_approx_dump)
    data.κ_approx =
        Arblib.load_string!.(zeros(Arb, size(data_dump, 1)), data_dump.κ_approx_dump)

    data.ξ₁ = Arblib.load_string!.(zeros(Arb, size(data_dump, 1)), data_dump.ξ₁_dump)

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

    data = DataFrame()

    data.ϵ_lower =
        Arblib.load_string!.(zeros(Arf, size(data_dump, 1)), data_dump.ϵ_lower_dump)
    data.ϵ_upper =
        Arblib.load_string!.(zeros(Arf, size(data_dump, 1)), data_dump.ϵ_upper_dump)

    data.μ_uniq =
        Arblib.load_string!.(zeros(Arb, size(data_dump, 1)), data_dump.μ_uniq_dump)
    data.γ_uniq =
        Acb.(
            Arblib.load_string!.(
                zeros(Arb, size(data_dump, 1)),
                data_dump.γ_uniq_dump_real,
            ),
            Arblib.load_string!.(
                zeros(Arb, size(data_dump, 1)),
                data_dump.γ_uniq_dump_imag,
            ),
        )
    data.κ_uniq =
        Arblib.load_string!.(zeros(Arb, size(data_dump, 1)), data_dump.κ_uniq_dump)

    data.μ_exists =
        Arblib.load_string!.(zeros(Arb, size(data_dump, 1)), data_dump.μ_exists_dump)
    data.γ_exists =
        Acb.(
            Arblib.load_string!.(
                zeros(Arb, size(data_dump, 1)),
                data_dump.γ_exists_dump_real,
            ),
            Arblib.load_string!.(
                zeros(Arb, size(data_dump, 1)),
                data_dump.γ_exists_dump_imag,
            ),
        )
    data.κ_exists =
        Arblib.load_string!.(zeros(Arb, size(data_dump, 1)), data_dump.κ_exists_dump)


    data.μ_approx =
        Arblib.load_string!.(zeros(Arb, size(data_dump, 1)), data_dump.μ_approx_dump)
    data.γ_approx =
        Acb.(
            Arblib.load_string!.(
                zeros(Arb, size(data_dump, 1)),
                data_dump.γ_approx_dump_real,
            ),
            Arblib.load_string!.(
                zeros(Arb, size(data_dump, 1)),
                data_dump.γ_approx_dump_imag,
            ),
        )
    data.κ_approx =
        Arblib.load_string!.(zeros(Arb, size(data_dump, 1)), data_dump.κ_approx_dump)

    data.ξ₁ = Arblib.load_string!.(zeros(Arb, size(data_dump, 1)), data_dump.ξ₁_dump)

    return data
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

    data = DataFrame()

    data.left_continuation = data_dump.left_continuation

    data.ϵ_lower =
        Arblib.load_string!.(zeros(Arf, size(data_dump, 1)), data_dump.ϵ_lower_dump)
    data.ϵ_upper =
        Arblib.load_string!.(zeros(Arf, size(data_dump, 1)), data_dump.ϵ_upper_dump)

    data.μ_uniq =
        Arblib.load_string!.(zeros(Arb, size(data_dump, 1)), data_dump.μ_uniq_dump)
    data.γ_uniq =
        Acb.(
            Arblib.load_string!.(
                zeros(Arb, size(data_dump, 1)),
                data_dump.γ_uniq_dump_real,
            ),
            Arblib.load_string!.(
                zeros(Arb, size(data_dump, 1)),
                data_dump.γ_uniq_dump_imag,
            ),
        )
    data.κ_uniq =
        Arblib.load_string!.(zeros(Arb, size(data_dump, 1)), data_dump.κ_uniq_dump)

    data.μ_exists =
        Arblib.load_string!.(zeros(Arb, size(data_dump, 1)), data_dump.μ_exists_dump)
    data.γ_exists =
        Acb.(
            Arblib.load_string!.(
                zeros(Arb, size(data_dump, 1)),
                data_dump.γ_exists_dump_real,
            ),
            Arblib.load_string!.(
                zeros(Arb, size(data_dump, 1)),
                data_dump.γ_exists_dump_imag,
            ),
        )
    data.κ_exists =
        Arblib.load_string!.(zeros(Arb, size(data_dump, 1)), data_dump.κ_exists_dump)


    data.μ_approx =
        Arblib.load_string!.(zeros(Arb, size(data_dump, 1)), data_dump.μ_approx_dump)
    data.γ_approx =
        Acb.(
            Arblib.load_string!.(
                zeros(Arb, size(data_dump, 1)),
                data_dump.γ_approx_dump_real,
            ),
            Arblib.load_string!.(
                zeros(Arb, size(data_dump, 1)),
                data_dump.γ_approx_dump_imag,
            ),
        )
    data.κ_approx =
        Arblib.load_string!.(zeros(Arb, size(data_dump, 1)), data_dump.κ_approx_dump)

    data.ξ₁ = Arblib.load_string!.(zeros(Arb, size(data_dump, 1)), data_dump.ξ₁_dump)

    return data
end
