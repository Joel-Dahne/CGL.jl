function branch_points_dataframe(
    λs::Vector{CGLParams{Arb}},
    μs::Vector{Arb},
    μs_approx::Vector{Float64},
    γs::Vector{Acb},
    κs::Vector{Arb},
    κs_approx::Vector{Float64},
    ξ₁s::Vector{Arb},
)
    df = DataFrame(
        ϵ = getproperty.(λs, :ϵ),
        μ = μs,
        μ_approx = μs_approx,
        γ = γs,
        κ = κs,
        κ_approx = κs_approx,
        ξ₁ = ξ₁s,
    )

    return df
end

function write_branch_points_csv(filename, data::DataFrame)
    data = copy(data)

    for col_name in names(data)
        col = data[!, col_name]
        if eltype(col) == Arb
            insertcols!(
                data,
                col_name,
                col_name * "_dump" => Arblib.dump_string.(col),
                after = true,
            )
        end
        if eltype(col) == Acb
            insertcols!(
                data,
                col_name,
                col_name * "_dump_real" => Arblib.dump_string.(real.(col)),
                after = true,
            )
            insertcols!(
                data,
                col_name * "_dump_real",
                col_name * "_dump_imag" => Arblib.dump_string.(imag.(col)),
                after = true,
            )
        end
    end

    insertcols!(data, :precision => precision.(data.ϵ))

    CSV.write(filename, data)
end

function read_branch_points_csv(filename)
    # We give the types explicitly, otherwise when the Arb values
    # happen to be exact they are parsed as floating points instead of
    # as strings.
    types = [
        String,
        String,
        String,
        String,
        Float64,
        String,
        String,
        String,
        String,
        String,
        Float64,
        String,
        String,
        Int,
    ]

    data = CSV.read(filename, DataFrame; types)

    # Handle Arb dumps
    for col_name in names(data)
        if endswith(col_name, "_dump")
            target_col_name = replace(col_name, "_dump" => "")

            col = data[!, col_name]
            data[!, target_col_name] =
                Arblib.load_string!.([Arb(; prec) for prec in data.precision], col)
            select!(data, Not(col_name))
        end
    end

    # Handle Acb dump for γ
    data[!, :γ] =
        Acb.(
            Arblib.load_string!.(
                [Arb(; prec) for prec in data.precision],
                data.γ_dump_real,
            ),
            Arblib.load_string!.(
                [Arb(; prec) for prec in data.precision],
                data.γ_dump_imag,
            ),
        )
    select!(data, Not(:γ_dump_real))
    select!(data, Not(:γ_dump_imag))

    # Drop precision
    select!(data, Not(:precision))

    return data
end
