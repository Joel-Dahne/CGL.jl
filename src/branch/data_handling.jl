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

function branch_dataframe(
    ϵs::Vector{NTuple{2,Arf}},
    μs_uniq::Vector{Arb},
    γs_uniq::Vector{Acb},
    κs_uniq::Vector{Arb},
    μs_exists::Vector{Arb},
    γs_exists::Vector{Acb},
    κs_exists::Vector{Arb},
    ξ₁s::Vector{Arb},
)
    df = DataFrame(
        ϵ_lower = Arb.(getindex.(ϵs, 1)),
        ϵ_upper = Arb.(getindex.(ϵs, 2)),
        μ_uniq = μs_uniq,
        γ_uniq = γs_uniq,
        κ_uniq = κs_uniq,
        μ_exists = μs_exists,
        γ_exists = γs_exists,
        κ_exists = κs_exists,
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

function write_branch_csv(filename, data::DataFrame)
    data_raw = DataFrame()

    for col_name in names(data)
        col = data[!, col_name]
        if eltype(col) == Arb
            insertcols!(data_raw, col_name * "_dump" => Arblib.dump_string.(col))
        elseif eltype(col) == Acb
            insertcols!(
                data_raw,
                col_name * "_dump_real" => Arblib.dump_string.(real.(col)),
                col_name * "_dump_imag" => Arblib.dump_string.(imag.(col)),
            )
        else
            insertcols!(data_raw, col_name => col)
        end
    end

    insertcols!(data_raw, :precision => precision.(data.ϵ_lower))

    CSV.write(filename, data_raw)
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

function read_branch_csv(filename)
    # We give the types explicitly, otherwise when the Arb values
    # happen to be exact they are parsed as floating points instead of
    # as strings.
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
        Int,
    ]

    data_raw = CSV.read(filename, DataFrame; types)
    data = DataFrame()

    for col_name_raw in names(data_raw)
        col_raw = data_raw[!, col_name_raw]

        if endswith(col_name_raw, "_dump")
            col_name = chopsuffix(col_name_raw, "_dump")
            col =
                Arblib.load_string!.([Arb(; prec) for prec in data_raw.precision], col_raw)

            insertcols!(data, col_name => col)
        elseif endswith(col_name_raw, "_dump_real")
            col_name = chopsuffix(col_name_raw, "_dump_real")
            col_raw_imag = data_raw[!, replace(col_name_raw, "_real" => "_imag")]
            col =
                Acb.(
                    Arblib.load_string!.(
                        [Arb(; prec) for prec in data_raw.precision],
                        col_raw,
                    ),
                    Arblib.load_string!.(
                        [Arb(; prec) for prec in data_raw.precision],
                        col_raw_imag,
                    ),
                )

            insertcols!(data, col_name => col)
        end
    end

    return data

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
    for col_name in ["γ_uniq", "γ_exists"]
        col_name_dump_real = col_name * "_dump_real"
        col_name_dump_imag = col_name * "_dump_imag"

        data[!, col_name] =
            Acb.(
                Arblib.load_string!.(
                    [Arb(; prec) for prec in data.precision],
                    data[!, col_name_dump_real],
                ),
                Arblib.load_string!.(
                    [Arb(; prec) for prec in data.precision],
                    data[!, col_name_dump_imag],
                ),
            )
        select!(data, Not(col_name_dump_real))
        select!(data, Not(col_name_dump_imag))
    end

    # Drop precision
    select!(data, Not(:precision))

    return data
end
