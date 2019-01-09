module ContinuousTimeMarkov

export TransitionRateMatrix, stationary_distribution

using ArgCheck: @argcheck
using DocStringExtensions: SIGNATURES
using Parameters: @unpack
using LinearAlgebra: Diagonal, lu, normalize!
using SparseArrays: rowvals, nonzeros, nzrange, SparseMatrixCSC

####
#### Transition rate matrices
####

"""
$(SIGNATURES)

Return true iff off-diagonal elements of `A` are `≥ 0`.

Used internally, not exported.
"""
function is_nonnegative_offdiagonal(A::SparseMatrixCSC)
    rows = rowvals(A)
    vals = nonzeros(A)
    m, n = size(A)
    for i in 1:n
        for j in nzrange(A, i)
            rows[j] ≠ i && vals[j] < 0 && return false
        end
    end
    true
end

function is_nonnegative_offdiagonal(A::AbstractMatrix)
    for j in axes(A, 2)
        for i in axes(A, 1)     # column-major traversal
            i ≠ j && A[i, j] < 0 && return false
        end
    end
    true
end

is_square_matrix(A::AbstractMatrix) = axes(A, 1) == axes(A, 2)

struct TransitionRateMatrix{T, S <: AbstractMatrix} <: AbstractMatrix{T}
    matrix::S
    # NOTE: this is the unchecked constructor. Don't use unless you are sure that
    # 1. `matrix` is square, with the same `axis` 1 and 2,
    # 2. its rows sum to `0`,
    # 3. it will not be modified later (ie does not share structure).
    function TransitionRateMatrix(::Val{:trust}, matrix::S) where {T, S <: AbstractMatrix{T}}
        new{T, S}(matrix)
    end
end

Base.size(m::TransitionRateMatrix) = size(m.matrix)

Base.getindex(m::TransitionRateMatrix, I...) = Base.getindex(m.matrix, I...)

Base.IndexStyle(::Type{TransitionRateMatrix{T,S}}) where {T,S} = IndexStyle(S)

"""
$(SIGNATURES)

Create a transition rate matrix from the argument. This makes a copy and sets the diagonal
so that rows sum to `0`, striving to preserve sparsity and structure.

When the argument is already a transition rate matrix (all non-diagonal elements
nonnegative, rows sum to `0`), use
```julia
TransitionRateMatrix(Val(:trust), matrix)
```
"""
function TransitionRateMatrix(Q)
    @argcheck is_square_matrix(Q)
    @argcheck is_nonnegative_offdiagonal(Q)
    TransitionRateMatrix(Val{:trust}(), Q - Diagonal(vec(sum(Q; dims = 2))))
end

function stationary_distribution(m::TransitionRateMatrix)
    @unpack matrix = m
    LU = lu(matrix, Val(false); check = false) # not pivoted, singular
    L = LU.L
    x = vcat(L[1:(end-1), 1:(end-1)]' \ -L[end, 1:(end-1)], 1)
    normalize!(x, 1)
    x
end

end # module
