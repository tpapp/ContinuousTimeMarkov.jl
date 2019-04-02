module ContinuousTimeMarkov

export TransitionRateMatrix, TransitionRateMatrix!, stationary_distribution

using AlgebraResultTypes: result_ring
using ArgCheck: @argcheck
using DocStringExtensions: SIGNATURES
using Parameters: @unpack
using LinearAlgebra: Diagonal, lu, normalize!
using SparseArrays: rowvals, nonzeros, nzrange, SparseMatrixCSC

####
#### utilities
####

"""
$(SIGNATURES)

Assert that both axes are the same, and return a single one.

Internal, not part of the API. Works with generalized indexing.
"""
function checked_square_axis(A::AbstractMatrix)
    a1, a2 = axes(A)
    @argcheck a1 == a2 "Mismatching axes for matrix."
    a1
end

"""
$(SIGNATURES)

Check that off-diagonal elements are nonnegative, and return their sum by row.

Internal, not part of the API. Works with generalized indexing.
"""
function off_diagonal_checked_row_sum(A::AbstractMatrix{T}) where {T}
    _off_diagonal_checked_row_sum!(zeros(result_ring(T), checked_square_axis(A)), A)
end

"""
$(SIGNATURES)

Accumulate the sum of off-diagonal elements into `d` by row. Return `d`.

Internal, not part of the API. Works with generalized indexing.
"""
function _off_diagonal_checked_row_sum!(d, A::AbstractMatrix)
    @inbounds for ι in CartesianIndices(A)
        if ι[1] ≠ ι[2]
            a = A[ι]
            @argcheck a ≥ 0 DomainError(a, "Element at $(ι) is negative.")
            d[ι[1]] += a
        end
    end
    d
end

function _off_diagonal_checked_row_sum!(d, A::SparseMatrixCSC)
    rows = rowvals(A)
    vals = nonzeros(A)
    m, n = size(A)
    for j in 1:n
        for k in nzrange(A, j)
            i = rows[k]
            if i ≠ j
                a = vals[k]
                @argcheck a ≥ 0 DomainError(a, "Element at $(CartesianIndex(i,j)) is negative.")
                d[i] += a
            end
        end
    end
    d
end

"""
$(SIGNATURES)

Set the diagonal of the argument (which is modified) so that rows sum to `0`.

Internal, not part of the API. Works with generalized indexing.
"""
function normalize_diagonal!(A::AbstractMatrix)
    d = off_diagonal_checked_row_sum(A)
    @inbounds for i in checked_square_axis(A)
        A[i, i] = -d[i]
    end
    A
end

####
#### Transition rate matrices
####

struct TransitionRateMatrix{T, S <: AbstractMatrix} <: AbstractMatrix{T}
    matrix::S
    # NOTE: this is the unchecked constructor. Don't use unless you are sure that
    # 1. `matrix` is square, with the same `axis` 1 and 2,
    # 2. its rows sum to `0`,
    # 3. it will not be modified later (ie does not share structure).
    function TransitionRateMatrix(::Val{:trust}, matrix::S) where {T, S <: AbstractMatrix{T}}
        checked_square_axis(matrix)
        new{T, S}(matrix)
    end
end

Base.size(m::TransitionRateMatrix) = size(m.matrix)

Base.getindex(m::TransitionRateMatrix, I...) = Base.getindex(m.matrix, I...)

Base.IndexStyle(::Type{TransitionRateMatrix{T,S}}) where {T,S} = IndexStyle(S)

"""
$(SIGNATURES)

Create a transition rate matrix from the argument. This makes a copy, checks that
off-diagonal elements are nonnegative, and sets the diagonal so that rows sum to `0`,
striving to preserve sparsity and structure.
"""
function TransitionRateMatrix(Q::AbstractMatrix)
    TransitionRateMatrix(Val{:trust}(), normalize_diagonal!(copy(Q)))
end

"""
$(SIGNATURES)

Create a `TransitionRateMatrix` from a matrix `Q`. Modifies the argument (to normalize the
diagonal).
"""
function TransitionRateMatrix!(Q::AbstractMatrix)
    TransitionRateMatrix(Val{:trust}(), normalize_diagonal!(Q))
end

"""
$(SIGNATURES)

Return the stationary distribution of a transition rate matrix as a vector.
"""
function stationary_distribution(m::TransitionRateMatrix)
    @unpack matrix = m
    @assert !Base.has_offset_axes(matrix) "Generalized indexing version not implemented."
    LU = lu(matrix, Val(false); check = false) # not pivoted, allow singular
    L = LU.L
    x = vcat(L[1:(end-1), 1:(end-1)]' \ -L[end, 1:(end-1)], 1)
    normalize!(x, 1)
    x
end

end # module
