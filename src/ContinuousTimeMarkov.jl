module ContinuousTimeMarkov

export TransitionRateMatrix, TransitionRateMatrix!, stationary_distribution

using AlgebraResultTypes: result_ring
using ArgCheck: @argcheck
using DocStringExtensions: SIGNATURES
using Parameters: @unpack
using LinearAlgebra: Diagonal, lu, normalize!
using SparseArrays: rowvals, nonzeros, nzrange, SparseMatrixCSC

include("transition_matrices.jl")

end # module
