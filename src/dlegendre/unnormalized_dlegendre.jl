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
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

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
function. This can be computed using the function [`unnormalized_legendre!`](@ref).

!!! warning
    The user is responsible to pass a matrix `P` with the correct values. For example, if
    `ph_term` is `true`, `P` must also be computed with `ph_term` set to `true`.

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

    # The derivative is computed using the following equation [1, p. 1981]:
    #
    #   ∂P(n, m)
    #   ──────── = ¹/₂ . ((n + m) . (n - m + 1) . P(n, m-1) - P(n, m+1)),
    #      ∂θ
    #
    # TODO: This algorithm is based on eq. Z.1.44 of [2]. However, it was verified that it
    # does not provide the correct sign when ϕ ∈ [π, 2π]. This makes sense because the
    # algorithm uses only the values of the coefficients, which are equal for ϕ and -ϕ.
    # However, the derivative w.r.t.  ϕ does change. This hack was used so that the values
    # are correct, but further verification is needed.
    #
    # In fact, in [2, p. 119], it is mentioned that `0 <= ϕ <= π`.  However, further
    # clarification is required.

    ϕ    = mod(ϕ, T(2π))
    fact = ϕ > π ? -1 : 1

    if ph_term
        fact *= -1
    end

    # Get the first indices in `P` to take into account offset arrays.
    i₀, j₀ = first.(axes(P))

    # Get the first indices in `dP` to take into account offset arrays.
    di₀, dj₀ = first.(axes(dP))

    dP[di₀, dj₀] = 0

    m_max < 0 && return nothing

    @inbounds for n in 1:n_max
        for m in 0:n
            if m == 0
                dP[di₀+n, dj₀+m] = -P[i₀+n, j₀+1]
            elseif n != m
                dP[di₀+n, dj₀+m] = ((n + m) * (n - m + 1) * P[i₀+n, j₀+m-1] - P[i₀+n, j₀+m+1]) / 2
            else
                dP[di₀+n, dj₀+m] = ((n + m) * (n - m + 1) * P[i₀+n, j₀+m-1]) / 2
            end

            dP[di₀+n, dj₀+m] *= fact

            # Check if the maximum desired order has been reached.
            m >= m_max && break
        end
    end

    nothing
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

    # Check if we need to compute and additional degree in `P` to provide the desire order
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
