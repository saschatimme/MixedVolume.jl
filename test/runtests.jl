using MixedVolume
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

@testset "CayleyConfiguration" begin
    A1 = [0 1 0; 0 0 1]
    A2 = [1 0; 0 1]

    CC = CayleyConfiguration(A1, A2)
    @test CC.A == [0 1 0 1 0; 0 0 1 0 1; 1 1 1 0 0; 0 0 0 1 1]
    @test CC.offsets == [1, 4]
end
