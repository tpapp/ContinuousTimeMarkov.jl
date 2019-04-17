using ContinuousTimeMarkov, Test, LinearAlgebra, SparseArrays, Distributions, StatsBase

import Random

Random.seed!(1)

include("test_transition_matrices.jl")
include("test_competing_poisson.jl")
