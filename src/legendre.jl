## Description #############################################################################
#
# Functions related to the associated Legendre functions.
#
## References ##############################################################################
#
# [1] Holmes, S. A. and W. E. Featherstone, 2002. A unified approach to the Clenshaw
#     summation and the recursive computation of very high degree and order normalised
#     associated Legendre functions. Journal of Geodesy, 76(5), pp. 279-299.
#
#     For more info.: http://mitgcm.org/~mlosch/geoidcookbook/node11.html
#
# [2] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. Microcosm Press,
#     Hawthorn, CA, USA.
#
# [3] Schmidt, A (1917). Erdmagnetismus, Enzykl. Math. Wiss., 6, pp. 265–396.
#
# [4] Winch, D. E., Ivers, D. J., Turner, J. P. R., Stening R. J (2005). Geomagnetism and
#     Schmidt quasi-normalization. Geophysical Journal International, 160(2), pp. 487-504.
#
############################################################################################

export legendre!, legendre

"""
    legendre!(N, P::AbstractMatrix, ϕ::Number, n_max::Integer = -1, m_max::Integer = -1; kwargs...) -> Nothing

Compute the associated Legendre function `P_n,m[cos(ϕ)]`. The maximum degree and order that
will be computed are given by the parameters `n_max` and `m_max`. If they are negative
(default), the dimensions of matrix `P` will be used:

    maximum degree -> number of rows - 1
    maximum order  -> number of columns - 1

The result will be stored in matrix `P`.

!!! note

    This function only writes to the elements of `P` in the lower triangular part up to
    the computed degree and order. If `P` is being reused with a smaller `n_max` or
    `m_max`, the remaining elements will keep their previous (stale) values.

The parameter `N` selects the normalization. The following values are valid:

- `Val(:full)`: Compute the fully normalized associated Legendre function.
- `Val(:schmidt)`: Compute the Schmidt quasi-normalized associated Legendre function.
- `Val(:unnormalized)`: Compute the unnormalized associated Legendre function.

# Keywords

- `ph_term::Bool`: If `true`, the Condon-Shortley phase term `(-1)^m` will be included.
    (**Default** = `false`)

# Remarks

## Full normalization

This algorithm was based on **[1]**. Our definition of fully normalized associated Legendre
function can be seen in **[2, p. 546]**. The conversion is obtained by:

             ┌                       ┐
             │  (n-m)! . k . (2n+1)  │
    K_n,m = √│ ───────────────────── │,  k = (m = 0) ? 1 : 2.
             │         (n+m)!        │
             └                       ┘

    P̄_n,m = P_n,m * K_n,m,

where `P̄_n,m` is the fully normalized Legendre associated function.

## Schmidt quasi-normalization

This algorithm was based on **[3, 4]**. The conversion is obtained by:

             ┌              ┐
             │      (n-m)!  │
    K_n,m = √│ k . ──────── │,  k = (m = 0) ? 1 : 2.
             │      (n+m)!  │
             └              ┘

    P̂_n,m = P_n,m * K_n,m,

where `P̂_n,m` is the Schmidt quasi-normalized Legendre associated function.

# References

- **[1]** Holmes, S. A. and W. E. Featherstone, 2002. A unified approach to the Clenshaw
    summation and the recursive computation of very high degree and order normalised
    associated Legendre functions. Journal of Geodesy, 76(5), pp. 279-299. For more info.:
    http://mitgcm.org/~mlosch/geoidcookbook/node11.html

- **[2]** Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. Microcosm
    Press, Hawthorn, CA, USA.

- **[3]** Schmidt, A (1917). Erdmagnetismus, Enzykl. Math. Wiss., 6, pp. 265–396.

- **[4]** Winch, D. E., Ivers, D. J., Turner, J. P. R., Stening R. J (2005). Geomagnetism
    and Schmidt quasi-normalization. Geophysical Journal International, 160(2), pp. 487-504.
"""
function legendre!(
    N::Union{Val{:full}, Val{:schmidt}, Val{:unnormalized}},
    P::AbstractMatrix{T},
    ϕ::Number,
    n_max::Integer = -1,
    m_max::Integer = -1;
    ph_term::Bool = false,
) where {T <: Number}
    # Obtain the maximum degree and order that must be computed.
    n_max, m_max = _get_degree_and_order(P, n_max, m_max)

    return _legendre_kernel!(P, ϕ, N, n_max, m_max, ph_term, T)
end

"""
    legendre(N, ϕ::T, n_max::Integer, m_max::Integer = -1; ph_term::Bool = false) where T<:Number -> Matrix{float(T)}

Compute the associated Legendre function `P_n,m[cos(ϕ)]`. The maximum degree that will be
computed is `n_max` and the maximum order is `m_max`. Notice that if `m_max` is higher than
`n_max` or negative (default), it is set to `n_max`.

The parameter `N` selects the normalization. The following values are valid:

- `Val(:full)`: Compute the fully normalized associated Legendre function.
- `Val(:schmidt)`: Compute the Schmidt quasi-normalized associated Legendre function.
- `Val(:unnormalized)`: Compute the unnormalized associated Legendre function.

# Keywords

- `ph_term::Bool`: If `true`, the Condon-Shortley phase term `(-1)^m` will be included.
    (**Default** = `false`)

# Returns

- `Matrix{float(T)}`: A matrix with the Legendre associated functions `P_n,m[cos(ϕ)]`.

# Remarks

## Full normalization

This algorithm was based on **[1]**. Our definition of fully normalized associated Legendre
function can be seen in **[2, p. 546]**. The conversion is obtained by:

             ┌                       ┐
             │  (n-m)! . k . (2n+1)  │
    K_n,m = √│ ───────────────────── │,  k = (m = 0) ? 1 : 2.
             │         (n+m)!        │
             └                       ┘

    P̄_n,m = P_n,m * K_n,m,

where `P̄_n,m` is the fully normalized Legendre associated function.

## Schmidt quasi-normalization

This algorithm was based on **[3, 4]**. The conversion is obtained by:

             ┌              ┐
             │      (n-m)!  │
    K_n,m = √│ k . ──────── │,  k = (m = 0) ? 1 : 2.
             │      (n+m)!  │
             └              ┘

    P̂_n,m = P_n,m * K_n,m,

where `P̂_n,m` is the Schmidt quasi-normalized Legendre associated function.

# References

- **[1]** Holmes, S. A. and W. E. Featherstone, 2002. A unified approach to the Clenshaw
    summation and the recursive computation of very high degree and order normalised
    associated Legendre functions. Journal of Geodesy, 76(5), pp. 279-299. For more info.:
    http://mitgcm.org/~mlosch/geoidcookbook/node11.html

- **[2]** Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. Microcosm
    Press, Hawthorn, CA, USA.

- **[3]** Schmidt, A (1917). Erdmagnetismus, Enzykl. Math. Wiss., 6, pp. 265–396.

- **[4]** Winch, D. E., Ivers, D. J., Turner, J. P. R., Stening R. J (2005). Geomagnetism
    and Schmidt quasi-normalization. Geophysical Journal International, 160(2), pp. 487-504.
"""
function legendre(
    N::Union{Val{:full}, Val{:schmidt}, Val{:unnormalized}},
    ϕ::T,
    n_max::Integer,
    m_max::Integer = -1;
    ph_term::Bool = false,
) where {T <: Number}
    n_max < 0 && throw(ArgumentError("n_max must not be negative."))

    if (m_max < 0) || (m_max > n_max)
        m_max = n_max
    end

    P = zeros(float(T), n_max + 1, m_max + 1)
    legendre!(N, P, ϕ; ph_term = ph_term)

    return P
end

"""
    legendre!(P::AbstractMatrix, ϕ::Number, coefs::LegendreCoefficients{N, T}, n_max::Integer = -1, m_max::Integer = -1; kwargs...) where {N, T<:AbstractFloat} -> Nothing

Compute the associated Legendre function `P_n,m[cos(ϕ)]` using the precomputed coefficients
`coefs`, which also select the normalization (see [`LegendreCoefficients`](@ref)). The
maximum degree and order that will be computed are given by the parameters `n_max` and
`m_max`. If they are negative (default), the dimensions of matrix `P` will be used. This
function throws an `ArgumentError` if the maximum degree or order exceeds the ones
supported by `coefs`.

The result will be stored in matrix `P` and all the arithmetic operations are performed
using the element type `T` of the coefficients.

!!! note

    This function only writes to the elements of `P` in the lower triangular part up to
    the computed degree and order. If `P` is being reused with a smaller `n_max` or
    `m_max`, the remaining elements will keep their previous (stale) values.

# Keywords

- `ph_term::Bool`: If `true`, the Condon-Shortley phase term `(-1)^m` will be included.
    (**Default** = `false`)
"""
function legendre!(
    P::AbstractMatrix,
    ϕ::Number,
    coefs::LegendreCoefficients{N, T},
    n_max::Integer = -1,
    m_max::Integer = -1;
    ph_term::Bool = false,
) where {N, T <: AbstractFloat}
    # Obtain the maximum degree and order that must be computed.
    n_max, m_max = _get_degree_and_order(P, n_max, m_max)

    if (n_max > coefs.n_max) || (m_max > coefs.m_max_P)
        throw(
            ArgumentError(
                "The coefficients support the maximum degree $(coefs.n_max) and order " *
                "$(coefs.m_max_P), but the computation requires degree $n_max and order " *
                "$m_max.",
            ),
        )
    end

    return _legendre_kernel!(P, ϕ, coefs, n_max, m_max, ph_term, T)
end
