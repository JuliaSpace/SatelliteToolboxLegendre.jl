## Description #############################################################################
#
# Definition of the types used in the package.
#
############################################################################################

export LegendreCoefficients

"""
    struct LegendreCoefficients{N, T<:AbstractFloat}

Precomputed recursion coefficients to compute the associated Legendre function
`P_n,m[cos(ϕ)]` and its first-order derivative with respect to `ϕ` using the normalization
`N` (`:full`, `:schmidt`, or `:unnormalized`) and the element type `T`.

The coefficients depend only on the degree and order. Hence, when the functions are
evaluated at many angles `ϕ` with the same maximum degree and order, precomputing them
avoids evaluating square roots at every call, largely improving the performance.

# Fields

- `n_max::Int`: Maximum degree supported by the coefficients.
- `m_max::Int`: Maximum order supported by the derivative coefficients.
- `m_max_P::Int`: Maximum order supported by the Legendre function coefficients. It is one
    order higher than `m_max` (clamped to `n_max`) so the matrix `P` required by
    [`dlegendre!`](@ref) can be computed with the same object.
- `seed::T`: Coefficient of the terms with degree 1.
- `diag::Vector{T}`: Coefficients of the diagonal recursion (`n == m`).
- `a::Matrix{T}`: Coefficients `a_nm` of the Legendre function recursion.
- `b::Matrix{T}`: Coefficients `b_nm` of the Legendre function recursion.
- `da::Matrix{T}`: Coefficients `a_nm` of the derivative equation.
- `db::Matrix{T}`: Coefficients `b_nm` of the derivative equation.

# See Also

[`legendre!`](@ref), [`dlegendre!`](@ref)
"""
struct LegendreCoefficients{N, T<:AbstractFloat}
    n_max::Int
    m_max::Int
    m_max_P::Int
    seed::T
    diag::Vector{T}
    a::Matrix{T}
    b::Matrix{T}
    da::Matrix{T}
    db::Matrix{T}
end
