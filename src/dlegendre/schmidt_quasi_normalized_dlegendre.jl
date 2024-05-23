## Description #############################################################################
#
# Compute the first-order derivative of associated Legendre functions with the Schmidt
# quasi-normalization.
#
############################################################################################

"""
    schmidt_quasi_normalized_dlegendre!(dP::AbstractMatrix, ϕ::Number, P::AbstractMatrix, n_max::Integer = -1, m_max::Integer = -1; kwargs...) -> Nothing

Compute the first-order derivative of the Schmidt quasi-normalized associated Legendre
function `P_n,m[cos(ϕ)]` with respect to `ϕ` [rad]:

    ∂P_n,m[cos(ϕ)]
    ──────────────
         ∂ϕ

The maximum degree and order that will be computed are given by the parameters `n_max` and
`m_max`. If they are negative (default), the dimensions of matrix `dP` will be used:

    maximum degree -> number of rows - 1
    maximum order  -> number of columns - 1

The derivatives will be stored in the matrix `dP`.

This algorithm needs the matrix `P` with the values of the Schmidt quasi-normalized
associated Legendre function. This can be computed using the function
[`schmidt_quase_normalized_legendre!`](@ref).

!!! warning

    The user is responsible to pass a matrix `P` with the correct values. For example, if
    `ph_term` is `true`, `P` must also be computed with `ph_term` set to `true`.

# Keywords

- `ph_term::Bool`: If `true`, the Condon-Shortley phase term `(-1)^m` will be included.
    (**Default** = `false`)
"""
function schmidt_quasi_normalized_dlegendre!(
    dP::AbstractMatrix,
    ϕ::Number,
    P::AbstractMatrix,
    n_max::Integer = -1,
    m_max::Integer = -1;
    ph_term::Bool = false
)
    # The algorithm to compute the first-order derivative using Schmidt normalization is
    # precisely the same as the one that computes using full normalization, given that `P`
    # has the correct coefficients.
    return fully_normalized_dlegendre!(dP, ϕ, P, n_max, m_max; ph_term = ph_term)
end

"""
    schmidt_quasi_normalized_dlegendre(ϕ::T, n_max::Integer, m_max::Integer = -1; kwargs...) where T<:Number -> Matrix{float(T)}, Matrix{float(T)}

Compute the first-order derivative of the Schmidt quasi-normalized associated Legendre
function `P_n,m[cos(ϕ)]` with respect to `ϕ` [rad]:

    ∂P_n,m[cos(ϕ)]
    ──────────────
         ∂ϕ

The maximum degree that will be computed is `n_max` and the maximum order is `m_max`. Notice
that if `m_max` is higher than `n_max` or negative (default), it is set to `n_max`.

# Keywords

- `ph_term::Bool`: If `true`, the Condon-Shortley phase term `(-1)^m` will be included.
    (**Default** = `false`)

# Returns

- `Matrix{float(T)}`: A matrix with the first-order derivative of the Legendre associated
    functions `P_n,m[cos(ϕ)]`.
- `Matrix{float(T)}`: A matrix with the Legendre associated functions `P_n,m[cos(ϕ)]`.
"""
function schmidt_quasi_normalized_dlegendre(
    ϕ::T,
    n_max::Integer,
    m_max::Integer = -1;
    ph_term::Bool = false
) where T<:Number
    if (m_max < 0) || (m_max > n_max)
        m_max = n_max
    end

    # Check if we need to compute and additional degree in `P` to provide the desire order
    # in `dP`.
    if n_max == m_max
        n_max_P = m_max_P = n_max
    else
        n_max_P = n_max
        m_max_P = m_max + 1
    end

    # First, compute the matrix with the associated Legendre functions.
    P = schmidt_quasi_normalized_legendre(ϕ, n_max_P, m_max_P; ph_term = ph_term)

    # Now, compute and return the derivative of the associated Legendre
    # functions.
    dP = zeros(float(T), n_max + 1, m_max + 1)
    schmidt_quasi_normalized_dlegendre!(dP, ϕ, P; ph_term = ph_term)

    return dP, P
end
