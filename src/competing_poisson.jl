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

# Note

The distribution of `duration` is `Exponential(sum(rates))`, so this process can also be
thought of as competing exponentials.
"""
function CompetingPoisson(rates::AbstractVector; events = axis(rates, 1))
    total_rate = sum(rates)
    CompetingPoisson(total_rate, Categorical(rates ./ total_rate), events)
end

function Base.rand(rng::AbstractRNG, cp::CompetingPoisson)
    @unpack total_rate, index_distribution, events = cp
    τ = rand(rng, Exponential(total_rate))
    i = rand(rng, index_distribution)
    (duration = τ, event = events[i])
end
