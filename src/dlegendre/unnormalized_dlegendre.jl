## Description #############################################################################
#
# Compute the first-order derivative of associated Legendre functions without
#   normalization.
#
## References ##############################################################################
#
# [1] Du, J., Chen, C., Lesur, V., and Wang, L (2015). Non-singular spherical harmonic
#     expressions of geomagnetic vector and gradient tensor fields in the local
#     north-oriented reference frame. Geoscientific Model Development, 8, pp. 1979-1990.
#
# [2] Ilk, K. H.: Ein Beitrag zur Dynamik ausgedehnter Körper-Gravitationswechselwirkung,
#     Deutsche Geodätische Kommission.  Reihe C, Heft Nr. 288, München, 1983.
#
############################################################################################

"""
    unnormalized_dlegendre!(dP::AbstractMatrix{T}, ϕ::Number, P::AbstractMatrix, n_max::Integer = -1, m_max::Integer = -1; kwargs...) where T<:Number -> Nothing

Compute the first-order derivative of the unnormalized (or conventional) associated Legendre
function `P_n,m[cos(ϕ)]` with respect to `ϕ` [rad]:

    ∂P_n,m[cos(ϕ)]
    ──────────────
         ∂ϕ

The maximum degree and order that will be computed are given by the parameters `n_max` and
`m_max`. If they are negative (default), the dimensions of matrix `dP` will be used:

    maximum degree -> number of rows - 1
    maximum order  -> number of columns - 1

The derivatives will be stored in the matrix `dP`.

This algorithm needs the matrix `P` with the values of the unnormalized associated Legendre
function. This can be computed using the function [`unnormalized_legendre!`](@ref). The
algorithm accesses the terms `P[n, m + 1]` when the order `m` is lower than the degree `n`.
Hence, the matrix `P` must have at least `n_max + 1` rows and `m_max + 2` columns if
`m_max < n_max`, or `m_max + 1` columns otherwise. This function throws an `ArgumentError`
if `P` does not have the required dimensions.

!!! warning
    The user is responsible for passing a matrix `P` with the correct values. For example,
    if `ph_term` is `true`, `P` must also be computed with `ph_term` set to `true`.

# Keywords

- `ph_term::Bool`: If `true`, the Condon-Shortley phase term `(-1)^m` will be included.
    (**Default** = `false`)
"""
function unnormalized_dlegendre!(
    dP::AbstractMatrix{T},
    ϕ::Number,
    P::AbstractMatrix,
    n_max::Integer = -1,
    m_max::Integer = -1;
    ph_term::Bool = false
) where T<:Number
    # Obtain the maximum degree and order that must be computed.
    n_max, m_max = _get_degree_and_order(dP, P, n_max, m_max)

    return _dlegendre_kernel!(dP, ϕ, P, Val(:unnormalized), n_max, m_max, ph_term, T)
end

"""
    unnormalized_dlegendre(ϕ::T, n_max::Integer, m_max::Integer = -1; ph_term::Bool = false) where T<:Number -> Matrix{float(T)}, Matrix{float(T)}

Compute the first-order derivative of the unnormalized (or conventional) associated Legendre
function `P_n,m[cos(ϕ)]` with respect to `ϕ` [rad]:

    ∂P_n,m[cos(ϕ)]
    ──────────────
         ∂ϕ

The maximum degree that will be computed is `n_max` and the maximum order is `m_max`. Notice
that if `m_max` is higher than `n_max` or negative, it is set to `n_max`.

# Keywords

- `ph_term::Bool`: If `true`, the Condon-Shortley phase term `(-1)^m` will be included.
    (**Default** = `false`)

# Returns

- `Matrix{float(T)}`: A matrix with the first-order derivative of the Legendre associated
    functions `P_n,m[cos(ϕ)]`.
- `Matrix{float(T)}`: A matrix with the Legendre associated functions `P_n,m[cos(ϕ)]`.
"""
function unnormalized_dlegendre(
    ϕ::T,
    n_max::Integer,
    m_max::Integer = -1;
    ph_term::Bool = false
) where T<:Number

    if (m_max < 0) || (m_max > n_max)
        m_max = n_max
    end

    # Check if we need to compute an additional degree in `P` to provide the desired order
    # in `dP`.
    if n_max == m_max
        n_max_P = m_max_P = n_max
    else
        n_max_P = n_max
        m_max_P = m_max + 1
    end

    # First, compute the matrix with the associated Legendre functions.
    P = unnormalized_legendre(ϕ, n_max_P, m_max_P; ph_term = ph_term)

    # Now, compute and return the derivative of the associated Legendre functions.
    dP = zeros(float(T), n_max + 1, m_max + 1)
    unnormalized_dlegendre!(dP, ϕ, P; ph_term = ph_term)

    return dP, P
end
