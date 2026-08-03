## Description #############################################################################
#
# Definition of the types used in the package.
#
############################################################################################

export LegendreCoefficients

# Container used to store the packed coefficients. In Julia versions 1.11 and above, we use
# `Memory` for better performance. In earlier versions, we use a standard `Vector`.
@static if VERSION >= v"1.11-"
    const PackedStorage{T} = Memory{T}
else
    const PackedStorage{T} = Vector{T}
end

"""
    struct LegendreCoefficients{N, T<:AbstractFloat}

Precomputed recursion coefficients to compute the associated Legendre function
`P_n,m[cos(ϕ)]` and its first-order derivative with respect to `ϕ` using the normalization
`N` (`:full`, `:schmidt`, or `:unnormalized`) and the element type `T`.

The coefficients depend only on the degree and order. Hence, when the functions are
evaluated at many angles `ϕ` with the same maximum degree and order, precomputing them
avoids evaluating square roots at every call, largely improving the performance.

The coefficients related to a degree and order are stored in packed vectors containing only
the lower triangular part of the coefficient set, clamped at the maximum order (see
`_packed_index`). In Julia versions 1.11 and above, the storage uses `Memory`. In earlier
versions, it uses a standard `Vector`. The fields are internal and must not be modified or
accessed directly.

# Fields

- `n_max::Int`: Maximum degree supported by the coefficients.
- `m_max::Int`: Maximum order supported by the derivative coefficients.
- `m_max_P::Int`: Maximum order supported by the Legendre function coefficients. It is one
    order higher than `m_max` (clamped to `n_max`) so the matrix `P` required by
    [`dlegendre!`](@ref) can be computed with the same object.
- `seed::T`: Coefficient of the terms with degree 1.
- `diag::PackedStorage{T}`: Coefficients of the diagonal recursion (`n == m`), stored per
    degree.
- `a::PackedStorage{T}`: Coefficients `a_nm` of the Legendre function recursion, packed
    with maximum order `m_max_P`.
- `b::PackedStorage{T}`: Coefficients `b_nm` of the Legendre function recursion, packed
    with maximum order `m_max_P`.
- `da::PackedStorage{T}`: Coefficients `a_nm` of the derivative equation, packed with
    maximum order `m_max`.
- `db::PackedStorage{T}`: Coefficients `b_nm` of the derivative equation, packed with
    maximum order `m_max`.

# See Also

[`legendre!`](@ref), [`dlegendre!`](@ref)
"""
struct LegendreCoefficients{N, T <: AbstractFloat}
    n_max::Int
    m_max::Int
    m_max_P::Int
    seed::T
    diag::PackedStorage{T}
    a::PackedStorage{T}
    b::PackedStorage{T}
    da::PackedStorage{T}
    db::PackedStorage{T}
end
