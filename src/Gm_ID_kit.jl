module Gm_ID_kit

using MAT
using Interpolations
using StaticArrays

export ParseMAT
function ParseMAT(file_name::AbstractString, key_name::AbstractString)
    file = matopen(file_name)
    return read(file, key_name)
end

# There are three usage modes:
# (1) Simple lookup of parameters at some given (L, VGS, VDS, VSB)
# (2) Lookup of arbitrary ratios of parameters, e.g. GM_ID, GM_CGG at given (L, VGS, VDS, VSB)
# (3) Cross-lookup of one ratio against another, e.g. GM_CGG for some GM_ID
function LookUp1(
    dict::AbstractDict{String,T},
    out_var::AbstractString,
    L::SVector{N_L,Td},
    VGS::SVector{N_VGS,Td},
    VDS::SVector{N_VDS,Td},
    VSB::SVector{N_VSB,Td},
) where {T,Td<:Real,N_L,N_VGS,N_VDS,N_VSB}
    # simple interpolation in mode 1
    data = dict[out_var]
    nodes = (vec(dict["L"]), vec(dict["VGS"]), vec(dict["VDS"]), vec(dict["VSB"]))
    A = data
    return interpolate(nodes, A, Gridded(Linear()))[L, VGS, VDS, VSB]
end

function LookUp2(
    dict::AbstractDict{String,T},
    out_var_numerator::AbstractString,
    out_var_denominator::AbstractString,
    L::SVector{N_L,Td},
    VGS::SVector{N_VGS,Td},
    VDS::SVector{N_VDS,Td},
    VSB::SVector{N_VSB,Td},
) where {T,Td<:Real,N_L,N_VGS,N_VDS,N_VSB}
    # simple interpolation in mode 2
    data = dict[out_var_numerator] ./ dict[out_var_denominator]
    nodes = (vec(dict["L"]), vec(dict["VGS"]), vec(dict["VDS"]), vec(dict["VSB"]))
    A = data
    return interpolate(nodes, A, Gridded(Linear()))[L, VGS, VDS, VSB]
end

function LookUp3(
    dict::AbstractDict{String,T},
    out_var_numerator::AbstractString,
    out_var_denominator::AbstractString,
    in_var_numerator::AbstractString,
    in_var_denominator::AbstractString,
    in_desired::SVector{N_in,Td},
    L::SVector{N_L,Td},
    VGS::SVector{N_VGS,Td},
    VDS::SVector{N_VDS,Td},
    VSB::SVector{N_VSB,Td},
    Warning::Bool,
) where {T,Td<:Real,N_in,N_L,N_VGS,N_VDS,N_VSB}
    # interpolation in mode 3
    out_data = dict[out_var_numerator] ./ dict[out_var_denominator]
    in_data = dict[in_var_numerator] ./ dict[in_var_denominator]
    # assemble x and y data, then find y values at desired x
    nodes = (vec(dict["L"]), vec(dict["VGS"]), vec(dict["VDS"]), vec(dict["VSB"]))
    A = in_data
    in = interpolate(nodes, A, Gridded(Linear()))[L, VGS, VDS, VSB]

    A = out_data
    out = interpolate(nodes, A, Gridded(Linear()))[L, VGS, VDS, VSB]


    # permute so that VGS dimension always comes first
    # NOTE: swaped their indices instead
    # dropdims(in; dims = (1, 4))
    # dropdims(out; dims = (3, 4))

    dim = size(in)
    in = reshape(in, (dim[1], dim[2]))
    dim = size(out)
    out = reshape(out, (dim[1], dim[2]))

    dim = size(in)
    output = zeros(dim[1], length(in_desired))

    for i = 1:dim[1]
        for j = 1:length(in_desired)
            m, idx = findmax(in[i, :])
            if maximum(in_desired[j]) > m && Warning
                @warn(
                    "look_up warning: $(in_var_numerator)/$(in_var_denominator) input larger than maximum! (output is NaN)"
                )
            end
            # If gm/ID is the input value, find maximum and limit search range to VGS values to the RIGHT
            if in_var_numerator == "GM" && in_var_denominator == "ID"
                in_right = in[i, idx:end]
                out_right = out[i, idx:end]
                output[i, j] = extrapolate(
                    interpolate(
                        reverse!(in_right),
                        reverse!(out_right),
                        FritschButlandMonotonicInterpolation(),
                    ),
                    NaN,
                )[in_desired[j]]

                # If gm/Cgg of gm/Cgs is the input value, find maximum and limit search range to VGS values to the LEFT
            elseif in_var_numerator == "GM" &&
                   (in_var_denominator == "CGG" || in_var_denominator == "CGS")
                in_left = in[i, 1:idx]
                out_left = out[i, 1:idx]

                output[i, j] = extrapolate(
                    interpolate(
                        reverse!(in_left),
                        reverse!(out_left),
                        FritschButlandMonotonicInterpolation(),
                    ),
                    NaN,
                )[in_desired[j]]

            else
                crossings = length(
                    findall(diff(sign(in(i, :) - in_desired(j) + eps(Float64))) .!= 0),
                )
                if crossings > 1
                    output = []
                    error(
                        "*** look_up: Error! There are multiple curve intersections.\n*** Try to reduce the search range by specifying the VGS vector explicitly.\n*** Example: look_up(nch, ID_W M_GDS gm_gds, GS nch[VGS][10:end])",
                    )
                    # TODO:
                    # disp(""
                    # figure(1000)
                    # plot(1:length(x(:,i)), x(:,i), "x" 1:length(x(:,i)), ones(1, length(x(:,i)))*xdesired(j));
                    return output
                end
                k = in[i, :]
                l = out[i, :]

                output[i, j] = extrapolate(
                    interpolate(
                        reverse!(k),
                        reverse!(l),
                        FritschButlandMonotonicInterpolation(),
                    ),
                    NaN,
                )[in_desired[j]]
            end
        end
    end
    return output
end

export LookUp
function LookUp(
    data::AbstractDict{String,T},
    out_var::AbstractString,
    in_var::AbstractString = "",
    in_val = 0.0;
    # default values for parameters
    L = minimum(data["L"]),
    VGS = data["VGS"],
    VDS = maximum(data["VDS"]) / 2.0,
    VSB = 0.0,
    # TODO: add # Method = "pchip",
    Warning = true,
) where {T}

    # check if desired output is a ratio of two parameters
    out_ratio = !isnothing(findfirst('_', out_var))
    # check if first variable argument is a ratio of two parameters
    var_ratio = !isnothing(findfirst('_', in_var))

    # determine usage mode
    if (out_ratio && var_ratio)
        mode = 3
    elseif (out_ratio && !var_ratio)
        mode = 2
    elseif (!out_ratio && !var_ratio)
        mode = 1
    else
        error("Invalid syntax or usage mode!")
    end

    # input is a ratio in mode 3
    if mode == 3
        out_numerator, out_denominator = split(out_var, '_', limit = 2, keepempty = false)
        in_numerator, in_denominator = split(in_var, '_', limit = 2, keepempty = false)
        return LookUp3(
            data,
            out_numerator,
            out_denominator,
            in_numerator,
            in_denominator,
            SVector(in_val...),
            SVector(L...),
            SVector(VGS...),
            SVector(VDS...),
            SVector(VSB...),
            Warning,
        )
    elseif mode == 2
        out_numerator, out_denominator = split(out_var, '_', limit = 2, keepempty = false)
        output = LookUp2(
            data,
            out_numerator,
            out_denominator,
            SVector(L...),
            SVector(VGS...),
            SVector(VDS...),
            SVector(VSB...),
        )
    else
        output = LookUp1(
            data,
            out_var,
            SVector(L...),
            SVector(VGS...),
            SVector(VDS...),
            SVector(VSB...),
        )
    end

    # FIXME: assume that only two vectors will coexist
    dim = size(output)
    output = reshape(output, (dim[1], dim[2]))
    return output
end

# There are two basic usage scenarios:
# (1) Lookup VGS with known voltage at the source terminal
# (2) Lookup VGS with unknown source voltage, e.g. when the source of the
# transistor is the tail node of a differential pair
export LookUpVGS
function LookUpVGS(
    data::Dict{String,Any};
    # default values for parameters
    L = minimum(data["L"]),
    VGB = NaN,
    GM_ID = NaN,
    ID_W = NaN,
    VDS = maximum(data["VDS"]) / 2.0,
    VDB = NaN,
    VSB = 0.0,
    # TODO: add # Method = "pchip",
)
    # determine usage mode
    if isnan(VGB[1]) && isnan(VDB[1])
        mode = 1
    elseif !isnan(VGB[1]) && !isnan(VDB[1])
        mode = 2
    else
        error("Invalid syntax or usage mode! Please type")
        output = []
        return output
    end

    # Check whether GM_ID or ID_W was passed to function
    if !isnan(ID_W[1])
        ratio_string = "ID_W"
        ratio_data = ID_W
    elseif !isnan(GM_ID[1])
        ratio_string = "GM_ID"
        ratio_data = GM_ID
    else
        error("look_upVGS: Invalid syntax or usage mode! Please type ")
        output = []
        return output
    end

    if mode == 1
        VGS = data["VGS"]
        ratio = LookUp(data, ratio_string, "VGS", VGS; VDS = VDS, VSB = VSB, L = L)

    else # mode == 2
        step = data["VGS"][1] - data["VGS"][2]
        VSB = (maximum(data.VSB):step:minimum(data.VSB))'
        VGS = VGB - VSB
        VDS = VDB - VSB
        ratio = LookUp(
            data,
            ratio_string,
            "VGS",
            VGS;
            VDS = VDS,
            VSB = VSB,
            L = ones(length(VSB), 1) * L,
        )
        idx = isfinite(ratio)
        ratio = ratio[idx]
        VGS = VGS[idx]
    end

    ratio = permutedims(ratio, [2, 1])

    # Interpolation loop
    s = size(ratio)
    output = fill(NaN, s[2], length(ratio_data))
    for j = 1:s[2]
        ratio_range = ratio[:, j]
        VGS_range = VGS
        # If gm/ID, find maximum and limit search range to VGS values to the right
        if ratio_string == "GM_ID"
            m, idx = findmax(ratio_range)
            VGS_range = VGS_range[idx:end]
            ratio_range = ratio_range[idx:end]
            if maximum(ratio_data) > m
                @warn ("LookUpVGS: GM_ID input larger than maximum!")
            end
        end
        output[j, :] .= extrapolate(
            interpolate(
                reverse!(ratio_range),
                reverse!(VGS_range),
                FritschButlandMonotonicInterpolation(),
            ),
            NaN,
        )[ratio_data]
    end
    # Force converting one element to a number
    if length(output) == 1
        output = output[1]
    end
    return output
end

end


