module Gm_ID_kit

using MAT
using Interpolations

export ParseMAT
function ParseMAT(file_name::AbstractString, key_name::AbstractString)
    file = matopen(file_name)
    return read(file, key_name)
end

# There are three usage modes:
# (1) Simple lookup of parameters at some given (L, VGS, VDS, VSB)
# (2) Lookup of arbitrary ratios of parameters, e.g. GM_ID, GM_CGG at given (L, VGS, VDS, VSB)
# (3) Cross-lookup of one ratio against another, e.g. GM_CGG for some GM_ID
export LookUp
function LookUp(
    data::Dict{String,Any},
    out_var::String,
    in_var::String = "",
    in_val = 0;
    # default values for parameters
    L = minimum(data["L"]),
    VGS = vec(data["VGS"]),
    VDS = maximum(data["VDS"]) / 2,
    VSB = 0,
    # TODO: add # Method = "pchip",
    Warning = true,
)

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

    # output is a ratio in modes 2 and 3
    if mode == 2 || mode == 3
        numerator, denominator = split(out_var, '_', limit = 2, keepempty = false)
        ydata = data[numerator] ./ data[denominator]
    else
        # simple output in mode 1
        ydata = data[out_var]
    end

    # input is a ratio in mode 3
    if mode == 3
        numerator, denominator = split(in_var, '_', limit = 2, keepempty = false)
        xdata = data[numerator] ./ data[denominator]
        xdesired = in_val
        # assemble x and y data, then find y values at desired x
        nodes = (vec(data["L"]), vec(data["VGS"]), vec(data["VDS"]), vec(data["VSB"]))
        A = xdata
        x = interpolate(nodes, A, Gridded(Linear()))[L, VGS, VDS, VSB]

        A = ydata
        y = interpolate(nodes, A, Gridded(Linear()))[L, VGS, VDS, VSB]

        # NOTE: not needed
        # permute so that VGS dimension always comes first
        # permutedims!(x, x, [2 1 3 4])
        # permutedims!(y, y, [2 1 3 4])
        # dropdims
        # squeez to get rid of 1-dims
        # return x
        if length(size(L)) == 0 # equivalent
            x = x'
            y = y'
        end

        dim = size(x)
        output = zeros(dim[1], length(xdesired))

        for i = 1:dim[1]
            for j = 1:length(xdesired)
                m, idx = findmax(x[i, :])
                if maximum(xdesired[j]) > m && Warning
                    @warn(
                        "look_up warning: $in_var input larger than maximum! (output is NaN)"
                    )
                end
                # If gm/ID is the x value, find maximum and limit search range to VGS values to the RIGHT
                if in_var == "GM_ID"
                    x_right = x[i, idx:end]
                    y_right = y[i, idx:end]
                    output[i, j] = extrapolate(
                        interpolate(
                            reverse!(x_right),
                            reverse!(y_right),
                            FritschButlandMonotonicInterpolation(),
                        ),
                        NaN,
                    )[xdesired[j]]

                    # If gm/Cgg of gm/Cgs is the x value, find maximum and limit search range to VGS values to the LEFT
                elseif numerator == "GM" && (denominator == "CGG" || denominator == "CGS")
                    x_left = x[i, 1:idx]
                    y_left = y[i, 1:idx]

                    output[i, j] = extrapolate(
                        interpolate(
                            reverse!(x_left),
                            reverse!(y_left),
                            FritschButlandMonotonicInterpolation(),
                        ),
                        NaN,
                    )[xdesired[j]]

                else
                    crossings = length(
                        findall(diff(sign(x(i, :) - xdesired(j) + eps(Float64))) .!= 0),
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
                    k = x[i, :]
                    l = y[i, :]

                    output[i, j] = extrapolate(
                        interpolate(
                            reverse!(k),
                            reverse!(l),
                            FritschButlandMonotonicInterpolation(),
                        ),
                        NaN,
                    )[xdesired[j]]
                end
            end
        end
    else
        # simple interpolation in modes 1 and 2
        if length(data["VSB"]) > 1
            nodes = (vec(data["L"]), vec(data["VGS"]), vec(data["VDS"]), vec(data["VSB"]))
            A = ydata
            output = interpolate(nodes, A, Gridded(Linear()))[L, VGS, VDS, VSB]
        else
            nodes = (vec(data["L"]), vec(data["VGS"]), vec(data["VDS"]), vec(data["VSB"]))
            A = ydata
            output = interpolate(nodes, A, Gridded(Linear()))[L, VGS, VDS, VSB]
        end
    end

    # Force column vector to matrix if the output is one dimensional
    if length(size(output)) == 1
        output = reshape(output, :, 1)
    end
    # Force converting one element to a number
    if length(output) == 1
        output = output[1]
    end

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
    VDS = maximum(data["VDS"]) / 2,
    VDB = NaN,
    VSB = 0,
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

    # Permutation needed if L is passed as a vector
    if length(L) > 1
        permute!(ratio, [2, 1])
    end

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


