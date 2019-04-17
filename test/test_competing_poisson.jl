@testset "competing Poisson" begin
    r = [0.2, 0.7, 3.9]
    cp = CompetingPoisson(r; events = 2:4)
    c = zeros(Int, 4)
    d = zeros(10000)
    for i in axes(d, 1)
        x = rand(cp)
        d[i] = x.duration
        c[x.event] += 1
    end
    e = ecdf(d)
    p = (1:9)./10
    q = quantile.(Exponential(sum(r)), p)
    @test all(isapprox.(e.(q), p; atol = 0.02))
    @test normalize(c, 1) ≈ vcat([0], normalize(r, 1)) atol = 0.02
end
