#####
##### Competing Poisson processes
#####

export CompetingPoisson

struct CompetingPoisson{T,D,E}
    total_rate::T
    index_distribution::D
    events::E
end

"""
$(SIGNATURES)

The first arrival from competing Poisson point processes with the given `rates`.

`rand` returns a `NamedTuple{(:duration,:event)}`, where `duration` is the duration until
the first event and `event` is drawn from `events`.

Varios properties are supported, see [`propertynames`](@ref).

# Note

The distribution of `duration` is `Exponential(1/sum(rates))`, so this process can also be
thought of as competing exponentials.
"""
function CompetingPoisson(rates::AbstractVector; events = axis(rates, 1))
    total_rate = sum(rates)
    CompetingPoisson(total_rate, Categorical(rates ./ total_rate), events)
end

function Base.rand(rng::AbstractRNG, sampler::SamplerTrivial{<:CompetingPoisson})
    @unpack total_rate, index_distribution, events = sampler[]
    τ = rand(rng, Exponential(1 / total_rate))
    i = rand(rng, index_distribution)
    (duration = τ, event = events[i])
end

function Base.propertynames(cp::CompetingPoisson, private = false)
    public = (:total_rate, :rates, :probabilities, :events)
    private ? (:index_distribution, public...) : public
end

function Base.getproperty(cp::CompetingPoisson, key::Symbol)
    if key ≡ :probabilities
        probs(getfield(cp, :index_distribution))
    elseif key ≡ :rates
        cp.probabilities .* cp.total_rate
    else
        getfield(cp, key)
    end
end

function Base.show(io::IO, cp::CompetingPoisson)
    print(io, "Competing Poisson with total rate $(round(cp.total_rate; sigdigits = 3));")
    for (e, p) in zip(cp.events, cp.probabilities)
        print(io, " P($(e))=$(round(p; sigdigits = 3))")
    end
end
