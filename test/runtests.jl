using Gm_ID_kit
using Test
# using Plots

@testset "LookUp test" begin
    nch = ParseMAT(joinpath(dirname(pathof(Gm_ID_kit)), "..", "test", "180nch.mat"), "nch")
    pch = ParseMAT(joinpath(dirname(pathof(Gm_ID_kit)), "..", "test", "180pch.mat"), "pch")

    # mode 1
    ID = LookUp(nch, "ID"; VGS = 0.6, L = 0.28)[1]
    @test ID ≈ 4.1232e-05 atol = 0.0001e-05
    ID = LookUp(pch, "ID"; VGS = 0.6, L = 0.28)[1]
    @test ID ≈ 1.5505e-05 atol = 0.0001e-05
    CGG = LookUp(nch, "CGG"; VGS = 0.6, L = 0.28)[1]
    @test CGG ≈ 1.3625e-14 atol = 0.0001e-14
    CDD = LookUp(nch, "CDD"; VGS = 0.6, L = 0.28)[1]
    @test CDD ≈ 6.0872e-15 atol = 0.0001e-15

    # mode 2
    GM_CGG = LookUp(nch, "GM_CGG"; VGS = 0.6, L = 0.28)[1]
    @test GM_CGG ≈ 4.2808e10 atol = 0.0001e10
    GM_ID = LookUp(nch, "GM_ID"; VGS = 0.6, L = 0.28)[1]
    @test GM_ID ≈ 14.1456 atol = 0.0001

    # mode 3
    GM_CDD = LookUp(nch, "GM_CDD", "GM_ID", 15.0; VDS = 0.6, L = 0.28)[1]
    @test GM_CDD ≈ 8.2656e10 atol = 0.0001e10
    GM_CDD = LookUp(pch, "GM_CDD", "GM_ID", 10.0; VDS = 0.6, L = 1.7)[1]
    @test GM_CDD ≈ 4.4941e09 atol = 0.0001e09
    GM_CGG = LookUp(nch, "GM_CGG", "GM_ID", 10.0; VDS = 0.6, L = 1.7)[1]
    @test GM_CGG ≈ 2.1313e09 atol = 0.0001e09
    GM_CGG = LookUp(pch, "GM_CGG", "GM_ID", 10.0; VDS = 0.6, L = 1.7)[1]
    @test GM_CGG ≈ 5.3293e08 atol = 0.0001e08
end

@testset "Example test" begin
    nch = ParseMAT(joinpath(dirname(pathof(Gm_ID_kit)), "..", "test", "180nch.mat"), "nch")
    pch = ParseMAT(joinpath(dirname(pathof(Gm_ID_kit)), "..", "test", "180pch.mat"), "pch")

    L = vec(nch["L"])

    mutable struct Mosfet
        L::Float64
        W::Float64
        ID::Float64
        GM_ID::Float64
        GM_GDS::Float64
        ID_W::Float64
        VDS::Float64
        VGS::Float64
        VSB::Float64
        Mosfet() = new(0, 0, 0, 0, 0, 0, 0, 0, 0)
    end

    # Compute VGS & ID
    # Givens
    M1 = Mosfet()
    VDD = 1.8
    M1.GM_ID = 10
    M1.GM_GDS = 50
    M1.VDS = VDD / 3
    M1.VSB = 0
    M1.W = 5

    # Find L that give gm/gds > given value
    # Get the gm/gds values vector corresponding to the L_vector
    gm_gds_vector = LookUp(nch, "GM_GDS", "GM_ID", M1.GM_ID; VDS = M1.VDS, L = L)
    # Get the minimum L that gives gm/gds > the given value
    # add line to get the minimum L for M1 that gives gm/gds >= M1.gm_gds
    M1.L = L[findfirst(x -> x > M1.GM_GDS, gm_gds_vector)]
    # Get the current by computing the ID/W and then multiply it by W\
    M1.ID_W = LookUp(nch, "ID_W", "GM_ID", M1.GM_ID; VDS = M1.VDS, L = M1.L)[1]
    # add line to get the current of M1
    M1.ID = M1.ID_W * M1.W
    # Get the VGS value
    # add line to get the VGS value of M1
    M1.VGS = LookUpVGS(nch; GM_ID = M1.GM_ID, VDS = M1.VDS, L = M1.L)
    @test M1.VGS ≈ 0.67 atol = 0.01
    @test M1.ID ≈ 9.52e-05 atol = 0.01e-05
end
