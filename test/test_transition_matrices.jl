@testset "transition rate matrix checks" begin
    @test_throws ArgumentError TransitionRateMatrix(ones(3, 4))  # non-square
    @test_throws DomainError TransitionRateMatrix(-ones(3, 3))   # negative
end

@testset "transition matrix creation (dense)" begin
    for i in 1:100
        n = rand(5:10)
        A = abs.(randn(n, n))
        # copying constructor
        Q = TransitionRateMatrix(A)
        @test size(Q) ≡ (n, n)
        @test Base.IndexStyle(typeof(Q)) ≡ Base.IndexStyle(typeof(A))
        for i in axes(Q, 1)
            @test sum(Q[i, :]) ≈ 0 atol = n*10*eps()
        end

        # non-copying constructor
        Q2 = TransitionRateMatrix!(A)
        @test size(Q2) ≡ (n, n)
        @test Base.IndexStyle(typeof(Q2)) ≡ Base.IndexStyle(typeof(A))
        for i in axes(Q2, 1)
            @test sum(Q2[i, :]) ≈ 0 atol = n*10*eps()
        end

        π = stationary_distribution(Q)
        π2 = stationary_distribution(Q2)
        @test π == π2           # exact equivalence
        @test sum(π) ≈ 1
        @test all(π .≥ 0)
        @test Q'*π ≈ zero(π) atol = 10*eps()
    end
end

@testset "transition matrix creation (sparse)" begin
    for i in 1:100
        n = rand(5:10)
        A = abs.(sprandn(n, n, 1/n))
        # copying constructor
        Q = TransitionRateMatrix(A)
        @test size(Q) ≡ (n, n)
        @test Base.IndexStyle(typeof(Q)) ≡ Base.IndexStyle(typeof(A))
        for i in axes(Q, 1)
            @test sum(Q[i, :]) ≈ 0 atol = n*10*eps()
        end

        # non-copying constructor
        Q2 = TransitionRateMatrix!(A)
        @test size(Q2) ≡ (n, n)
        @test Base.IndexStyle(typeof(Q2)) ≡ Base.IndexStyle(typeof(A))
        for i in axes(Q2, 1)
            @test sum(Q2[i, :]) ≈ 0 atol = n*10*eps()
        end

        # FIXME tests commented out, cf
        # https://github.com/tpapp/ContinuousTimeMarkov.jl/issues/2

        # π = stationary_distribution(Q)
        # π2 = stationary_distribution(Q2)
        # @test π == π2           # exact equivalence
        # @test sum(π) ≈ 1
        # @test all(π .≥ 0)
        # @test Q'*π ≈ zero(π) atol = 10*eps()
    end
end
