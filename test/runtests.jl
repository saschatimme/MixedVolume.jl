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
    @test CC.offsets == [0, 3]
    @test size(CC) == (4, 5)
end

@testset "CircuitTable" begin
    A1 = [0 0 1 1; 0 2 0 1]
    A2 = [0 0 1 2; 0 1 1 0]
    A = CayleyConfiguration(A1, A2)

    M = [(2, 3), (1, 3)]

    T = CircuitTable(A, M)
    @test T isa CircuitTable{Int}

    w = [0, 0, 0, -1, 0, -3, -4, -8]

    @test dot(T, 1, w) == 8
    @test dot(T, 2, w) == 0
    @test dot(T, 3, w) == 0
    @test dot(T, 4, w) == -1
    @test dot(T, 4, [0, 0, 0, -2, 0, -3, -4, -8]) == 2
    @test dot(T, 5, w) == 0
    @test dot(T, 7, w) == 0
end

@testset "MixedCell" begin
    A1 = [0 0 1 1; 0 2 0 1]
    A2 = [0 0 1 2; 0 1 1 0]
    A = CayleyConfiguration(A1, A2)

    M = MixedCell(A, [(2, 3), (1, 3)])
    @test M isa MixedCell{Int}
end
