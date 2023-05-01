# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Functions related to the first-order derivative of associated Legendre
#   functions.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] Du, J., Chen, C., Lesur, V., and Wang, L (2015). Non-singular spherical harmonic
#       expressions of geomagnetic vector and gradient tensor fields in the local
#       north-oriented reference frame. Geoscientific Model Development, 8, pp. 1979-1990.
#
#   [2] Ilk, K. H.: Ein Beitrag zur Dynamik ausgedehnter Körper-Gravitationswechselwirkung,
#       Deutsche Geodätische Kommission.  Reihe C, Heft Nr. 288, München, 1983.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export dlegendre!, dlegendre

"""
    dlegendre!(dP::AbstractMatrix, ϕ::Number, P::AbstractMatrix, n_max::Integer = -1, m_max::Integer = -1; kwargs...) -> Nothing

Compute the first-order derivative of the associated Legendre function `P_n,m[x]` with
respect to `ϕ` [rad]:

    ∂P_n,m[cos(ϕ)]
    ──────────────
         ∂ϕ

The maximum degree and order that will be computed are given by the parameters `n_max` and
`m_max`. If they are negative (default), the dimensions of matrix `dP` will be used:

    maximum degree -> number of rows - 1
    maximum order  -> number of columns - 1

The derivatives will be stored in the matrix `dP`.

The parameter `N` selects the normalization. The following values are valid:

- `Val(:full)`: Compute the fully normalized associated Legendre function.
- `Val(:schmidt)`: Compute the Schmidt quasi-normalized associated Legendre function.
- `Val(:unnormalized)`: Compute the unnormalized associated Legendre function.

This algorithm needs the matrix `P` with the values of the associated Legendre function
using the same normalization `N`, which can be computed using the function
[`legendre`](@ref).

!!! warning
    The user is responsible to pass a matrix `P` with the correct values. For example, if
    `ph_term` is `true`, `P` must also be computed with `ph_term` set to `true`.

# Keywords

- `ph_term::Bool`: If `true`, then the Condon-Shortley phase term `(-1)^m` will be included.
    (**Default** = `false`)
"""
function dlegendre!(
    ::Val{:full},
    dP::AbstractMatrix,
    ϕ::Number,
    P::AbstractMatrix,
    n_max::Integer = -1,
    m_max::Integer = -1;
    ph_term::Bool = false
)
    return fully_normalized_dlegendre!(dP, ϕ, P, n_max, m_max; ph_term = ph_term)
end

function dlegendre!(
    ::Val{:schmidt},
    dP::AbstractMatrix,
    ϕ::Number,
    P::AbstractMatrix,
    n_max::Integer = -1,
    m_max::Integer = -1;
    ph_term::Bool = false
)
    return schmidt_quasi_normalized_dlegendre!(dP, ϕ, P, n_max, m_max; ph_term = ph_term)
end

function dlegendre!(
    ::Val{:unnormalized},
    dP::AbstractMatrix,
    ϕ::Number,
    P::AbstractMatrix,
    n_max::Integer = -1,
    m_max::Integer = -1;
    ph_term::Bool = false
)
    return unnormalized_dlegendre!(dP, ϕ, P, n_max, m_max; ph_term = ph_term)
end

"""
    dlegendre(N, ϕ::T, n_max::Integer, m_max::Integer = -1; kwargs...) where T<:Number -> Matrix{float(T)}, Matrix{float(T)}

Compute the first-order derivative of the associated Legendre function `P_n,m[cos(ϕ)]`
with respect to `ϕ` [rad]:

    ∂P_n,m[cos(ϕ)]
    ──────────────
         ∂ϕ

The maximum degree that will be computed is `n_max` and the maximum order is `m_max`. Notice
that if `m_max` is higher than `n_max` or negative (default), it is set to `n_max`.

The parameter `N` selects the normalization. The following values are valid:

- `Val(:full)`: Compute the fully normalized associated Legendre function.
- `Val(:schmidt)`: Compute the Schmidt quasi-normalized associated Legendre function.
- `Val(:unnormalized)`: Compute the unnormalized associated Legendre function.

# Keywords

- `ph_term::Bool`: If `true`, then the Condon-Shortley phase term `(-1)^m` will be included.
    (**Default** = `false`)

# Returns

- `Matrix{float(T)}`: A matrix with the first-order derivative of the Legendre associated
    functions `P_n,m[cos(ϕ)]`.
- `Matrix{float(T)}`: A matrix with the Legendre associated functions `P_n,m[cos(ϕ)]`.
"""
function dlegendre(
    ::Val{:full},
    ϕ::Number,
    n_max::Integer,
    m_max::Integer = -1;
    ph_term::Bool = false
)
    return fully_normalized_dlegendre(ϕ, n_max, m_max; ph_term = ph_term)
end

function dlegendre(
    ::Val{:schmidt},
    ϕ::Number,
    n_max::Integer,
    m_max::Integer = -1;
    ph_term::Bool = false
)
    return schmidt_quasi_normalized_dlegendre(ϕ, n_max, m_max; ph_term = ph_term)
end

function dlegendre(
    ::Val{:unnormalized},
    ϕ::Number,
    n_max::Integer,
    m_max::Integer = -1;
    ph_term::Bool = false
)
    return unnormalized_dlegendre(ϕ, n_max, m_max; ph_term = ph_term)
end
