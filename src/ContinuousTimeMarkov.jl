module ContinuousTimeMarkov

export TransitionRateMatrix, TransitionRateMatrix!, stationary_distribution

using AlgebraResultTypes: result_ring
using ArgCheck: @argcheck
using Distributions: Categorical, Exponential
using DocStringExtensions: SIGNATURES
using Parameters: @unpack
using LinearAlgebra: Diagonal, lu, normalize!
using Random: AbstractRNG
using SparseArrays: rowvals, nonzeros, nzrange, SparseMatrixCSC

include("transition_matrices.jl")
include("competing_poisson.jl")

end # module
