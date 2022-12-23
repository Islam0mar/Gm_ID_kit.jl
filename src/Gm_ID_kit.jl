module Gm_ID_kit

using MAT
using Interpolations

export ParseMAT
function ParseMAT(file_name::AbstractString, key_name::AbstractString)
    file = matopen(file_name)
    return read(file, key_name)
end

function MakeVector(a::Matrix{T}) where {T<:AbstractFloat}
    return vec(a)
end

function MakeVector(a::Vector{T}) where {T<:AbstractFloat}
    return a
end

function MakeVector(a::T) where {T<:AbstractFloat}
    return [a]
end

# There are three usage modes:
# (1) Simple lookup of parameters at some given (L, VGS, VDS, VSB)
# (2) Lookup of arbitrary ratios of parameters, e.g. GM_ID, GM_CGG at given (L, VGS, VDS, VSB)
# (3) Cross-lookup of one ratio against another, e.g. GM_CGG for some GM_ID
# Output for the three usage modes:
# (1) Array corresponding to (L, VGS, VDS, VSB)
# (2) Array corresponding to (L, VGS, VDS, VSB)
# (3) Array corresponding to (L, VDS, VSB, IN_RATIO)
function LookUp1(
    dict::AbstractDict{String,T},
    out_var::AbstractString,
    L::Vector{Td},
    VGS::Vector{Td},
    VDS::Vector{Td},
    VSB::Vector{Td},
) where {T,Td<:AbstractFloat}
    # simple interpolation in mode 1
    data = dict[out_var]
    nodes = (vec(dict["L"]), vec(dict["VGS"]), vec(dict["VDS"]), vec(dict["VSB"]))
    A = data
    return extrapolate(interpolate(nodes, A, Gridded(Linear())), NaN)[L, VGS, VDS, VSB]
end

function LookUp2(
    dict::AbstractDict{String,T},
    out_var_numerator::AbstractString,
    out_var_denominator::AbstractString,
    L::Vector{Td},
    VGS::Vector{Td},
    VDS::Vector{Td},
    VSB::Vector{Td},
) where {T,Td<:AbstractFloat}
    # simple interpolation in mode 2
    data = dict[out_var_numerator] ./ dict[out_var_denominator]
    nodes = (vec(dict["L"]), vec(dict["VGS"]), vec(dict["VDS"]), vec(dict["VSB"]))
    A = data
    return extrapolate(interpolate(nodes, A, Gridded(Linear())), NaN)[L, VGS, VDS, VSB]
end

function LookUp3(
    dict::AbstractDict{String,T},
    out_var_numerator::AbstractString,
    out_var_denominator::AbstractString,
    in_var_numerator::AbstractString,
    in_var_denominator::AbstractString,
    in_desired::Vector{Td},
    L::Vector{Td},
    VGS::Vector{Td},
    VDS::Vector{Td},
    VSB::Vector{Td},
    Warning::Bool,
) where {T,Td<:AbstractFloat}
    # interpolation in mode 3
    out_data = dict[out_var_numerator] ./ dict[out_var_denominator]
    in_data = dict[in_var_numerator] ./ dict[in_var_denominator]
    # assemble x and y data, then find y values at desired x
    nodes = (vec(dict["L"]), vec(dict["VGS"]), vec(dict["VDS"]), vec(dict["VSB"]))
    A = in_data
    in = extrapolate(interpolate(nodes, A, Gridded(Linear())), NaN)[L, VGS, VDS, VSB]

    A = out_data
    out = extrapolate(interpolate(nodes, A, Gridded(Linear())), NaN)[L, VGS, VDS, VSB]

    dim = size(in)
    # dim[2] = VGS which is used while searching
    output = zeros(dim[1], dim[3], dim[4], length(in_desired))

    # Calculate output
    for i = 1:dim[1]
        for j = 1:dim[3]
            for k = 1:dim[4]
                for l = 1:length(in_desired)
                    m, idx = findmax(in[i, :, j, k])
                    idx = idx[1]
                    if maximum(in_desired[l]) > m && Warning
                        @warn(
                            "look_up warning: $(in_var_numerator)/$(in_var_denominator) input larger than maximum! (output is NaN)"
                        )
                    end
                    # If gm/ID is the input value, find maximum and limit search range to VGS values to the RIGHT
                    if in_var_numerator == "GM" && in_var_denominator == "ID"
                        in_right = in[i, idx:end, j, k]
                        out_right = out[i, idx:end, j, k]
                        output[i, j, k, l] = extrapolate(
                            interpolate(
                                reverse!(in_right),
                                reverse!(out_right),
                                FritschButlandMonotonicInterpolation(),
                            ),
                            NaN,
                        )[in_desired[l]]

                        # If gm/Cgg of gm/Cgs is the input value, find maximum and limit search range to VGS values to the LEFT
                    elseif in_var_numerator == "GM" &&
                           (in_var_denominator == "CGG" || in_var_denominator == "CGS")
                        in_left = in[i, 1:idx, j, k]
                        out_left = out[i, 1:idx, j, k]

                        output[i, j, k, l] = extrapolate(
                            interpolate(
                                reverse!(in_left),
                                reverse!(out_left),
                                FritschButlandMonotonicInterpolation(),
                            ),
                            NaN,
                        )[in_desired[l]]

                    else
                        crossings = length(
                            findall(
                                diff(
                                    sign(in(i, :, j, k) - in_desired(l) + eps(Float64)),
                                ) .!= 0,
                            ),
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
                        r = in[i, :, j, k]
                        s = out[i, :, j, k]

                        output[i, j, k, l] = extrapolate(
                            interpolate(
                                reverse!(r),
                                reverse!(s),
                                FritschButlandMonotonicInterpolation(),
                            ),
                            NaN,
                        )[in_desired[l]]
                    end
                end
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
    in_val::Td = 0.0;
    # default values for parameters
    L = Nothing,
    VGS = Nothing,
    VDS = Nothing,
    VSB = Nothing,
    # TODO: add # Method = "pchip",
    Warning::Bool = true,
) where {T,Td}
    NT = eltype(Td)
    in_val::Vector{NT} = MakeVector(in_val)
    L::Vector{NT} = MakeVector(L == Nothing ? minimum(data["L"]) : L)
    VGS::Vector{NT} = MakeVector(VGS == Nothing ? data["VGS"] : VGS)
    VDS::Vector{NT} = MakeVector(VDS == Nothing ? maximum(data["VDS"]) / NT(2) : VDS)
    VSB::Vector{NT} = MakeVector(VSB == Nothing ? NT(0) : VSB)

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
            MakeVector(in_val),
            MakeVector(L),
            MakeVector(VGS),
            MakeVector(VDS),
            MakeVector(VSB),
            Warning,
        )
    elseif mode == 2
        out_numerator, out_denominator = split(out_var, '_', limit = 2, keepempty = false)
        output = LookUp2(
            data,
            out_numerator,
            out_denominator,
            MakeVector(L),
            MakeVector(VGS),
            MakeVector(VDS),
            MakeVector(VSB),
        )
    else
        output = LookUp1(
            data,
            out_var,
            MakeVector(L),
            MakeVector(VGS),
            MakeVector(VDS),
            MakeVector(VSB),
        )
    end
    return output
end

# There are two basic usage scenarios:
# (1) Lookup VGS with known voltage at the source terminal
# (2) Lookup VGS with unknown source voltage, e.g. when the source of the
# transistor is the tail node of a differential pair
# Output for the two usage modes:
# (1) Array corresponding to (L, VGS, VDS, VSB)
# (2) Array corresponding to (L, VDS, VSB, IN_RATIO)
export LookUpVGS
function LookUpVGS(
    data::Dict{String,T};
    # default values for parameters
    L = minimum(data["L"]),
    VGB = NaN,
    GM_ID = NaN,
    ID_W = NaN,
    VDS = maximum(data["VDS"]) / 2.0,
    VDB = NaN,
    VSB = 0.0,
    # TODO: add # Method = "pchip",
) where {T}
    # determine usage mode
    if isnan(VGB[1]) && isnan(VDB[1])
        mode = 1
    elseif !isnan(VGB[1]) && !isnan(VDB[1])
        mode = 2
    else
        error("Invalid syntax or usage mode!")
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
        error("look_upVGS: Invalid syntax or usage mode!")
        output = []
        return output
    end

    if mode == 1
        VGS = vec(data["VGS"])
        ratio = LookUp(data, ratio_string; VGS = VGS, VDS = VDS, VSB = VSB, L = L)

    else # mode == 2
        step = data["VGS"][1] - data["VGS"][2]
        VSB = collect(maximum(data["VSB"]):step:minimum(data["VSB"]))
        VGS = VGB .- VSB
        VDS = VDB .- VSB
        ratio = Array{eltype(VSB)}(undef, (length(L),length(VGS),1,1))
        for i in 1:length(VSB)
            ratio[:,i,1,1] = LookUp(data, ratio_string; VGS = VGS[i], VDS = VDS[i], VSB = VSB[i], L = L)[:,1,1,1]
        end
        idx = vec(all(isfinite,ratio, dims=1))
        ratio = ratio[:,idx,:,:]
        VGS = vec(VGS[idx])
    end

    ratio = permutedims(ratio, [2, 1, 3 , 4])

    # Interpolation loop
    dim = size(ratio)
    output = fill(NaN, dim[2],dim[3],dim[4], length(ratio_data))
    for i = 1:dim[2]
        for j = 1:dim[3]
            for k = 1:dim[4]
                ratio_range = ratio[:, i,j,k]
                VGS_range = VGS
                # If gm/ID, find maximum and limit search range to VGS values to the right
                if ratio_string == "GM_ID"
                    m, idx = findmax(ratio_range)
                    VGS_range = reverse(VGS_range[idx:end])
                    ratio_range = reverse(ratio_range[idx:end])
                    if maximum(ratio_data) > m
                        @warn ("LookUpVGS: GM_ID input larger than maximum!")
                    end
                end
                output[i,j,k, :] .= extrapolate(
                    interpolate(
                        ratio_range,
                        VGS_range,
                        FritschButlandMonotonicInterpolation(),
                    ),
                    NaN,
        )[ratio_data]
            end
        end
    end
    return output
end

end


