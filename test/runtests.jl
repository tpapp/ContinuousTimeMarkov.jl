using ContinuousTimeMarkov
using Test
using LinearAlgebra, SparseArrays

@testset "transition rate matrix checks" begin
    @test_throws ArgumentError TransitionRateMatrix(ones(3, 4))  # non-square
    @test_throws ArgumentError TransitionRateMatrix(-ones(3, 3)) # negative
end

@testset "transition matrix creation" begin
    for _ in 1:100
        A = abs.(randn(3, 3))
        Q = TransitionRateMatrix(A)
        for i in axes(Q, 1)
            @test sum(Q[i, :]) ≈ 0 atol = 10*eps()
        end
        π = stationary_distribution(Q)
        @test sum(π) ≈ 1
        @test all(π .≥ 0)
        @test Q'*π ≈ zero(π) atol = 10*eps()
    end
end
