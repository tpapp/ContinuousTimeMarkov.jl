@testset "competing Poisson" begin
    r = [0.2, 0.7, 3.9]
    cp = CompetingPoisson(r; events = 2:4)

    # printing
    @test repr(cp) == "Competing Poisson with total rate 4.8; P(2)=0.0417 P(3)=0.146 P(4)=0.812"

    # properties
    @test cp.rates ≈ r
    @test cp.total_rate ≈ sum(r)
    @test cp.probabilities ≈ normalize(r, 1)
    @test cp.events == 2:4
    @test Set(propertynames(cp)) == Set((:total_rate, :rates, :probabilities, :events))
    @test Set(propertynames(cp, true)) == Set((propertynames(cp)..., :index_distribution))

    # empirical count
    c = zeros(Int, 4)           # event counts
    d = zeros(10000)            # durations
    for i in axes(d, 1)
        x = @inferred rand(cp)
        d[i] = x.duration
        c[x.event] += 1
    end
    e = ecdf(d)
    p = (1:9)./10
    q = quantile.(Exponential(1 / sum(r)), p)
    @test mean(d) ≈ (1 ./ sum(r)) atol = 0.01
    @test all(isapprox.(e.(q), p; atol = 0.02))
    @test normalize(c, 1) ≈ vcat([0], normalize(r, 1)) atol = 0.02
end
