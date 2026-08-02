## Description #############################################################################
#
# Compute the associated Legendre functions with full normalization.
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
############################################################################################

"""
    fully_normalized_legendre!(P::AbstractMatrix{T}, ϕ::Number, n_max::Integer = -1, m_max::Integer = -1; kwargs...) where T<:Number -> Nothing

Compute the fully normalized associated Legendre function `P_n,m[cos(ϕ)]`. The maximum
degree and order that will be computed are given by the arguments `n_max` and `m_max`. If
they are negative (default), the dimensions of matrix `P` will be used:

    maximum degree -> number of rows - 1
    maximum order  -> number of columns - 1

The result will be stored in the matrix `P`.

# Keywords

- `ph_term::Bool`: If `true`, the Condon-Shortley phase term `(-1)^m` will be included.
    (**Default** = `false`)

# Remarks

This algorithm was based on **[1]**. Our definition of fully normalized associated Legendre
function can be seen in **[2, p. 546]**. The conversion is obtained by:

             ┌                       ┐
             │  (n-m)! . k . (2n+1)  │
    K_n,m = √│ ───────────────────── │,  k = (m = 0) ? 1 : 2.
             │         (n+m)!        │
             └                       ┘

    P̄_n,m = P_n,m * K_n,m,

where `P̄_n,m` is the fully normalized Legendre associated function.

# References

- **[1]** Holmes, S. A. and W. E. Featherstone, 2002. A unified approach to the Clenshaw
    summation and the recursive computation of very high degree and order normalised
    associated Legendre functions. Journal of Geodesy, 76(5), pp. 279-299. For more info.:
    http://mitgcm.org/~mlosch/geoidcookbook/node11.html

- **[2]** Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. Microcosm
    Press, Hawthorn, CA, USA.
"""
function fully_normalized_legendre!(
    P::AbstractMatrix{T},
    ϕ::Number,
    n_max::Integer = -1,
    m_max::Integer = -1;
    ph_term::Bool = false
) where T<:Number
    # Obtain the maximum degree and order that must be computed.
    n_max, m_max = _get_degree_and_order(P, n_max, m_max)

    return _legendre_kernel!(P, ϕ, Val(:full), n_max, m_max, ph_term, T)
end

"""
    fully_normalized_legendre(ϕ::T, n_max::Integer, m_max::Integer = -1; kwargs...) where T<:Number -> Matrix{float(T)}

Compute the fully normalized associated Legendre function `P_n,m[cos(ϕ)]`. The maximum
degree that will be computed is `n_max` and the maximum order is `m_max`. Notice that if
`m_max` is higher than `n_max` or negative, it is set to `n_max`.

# Keywords

- `ph_term::Bool`: If `true`, the Condon-Shortley phase term `(-1)^m` will be included.
    (**Default** = `false`)

# Returns

- `Matrix{float(T)}`: A matrix with the Legendre associated functions `P_n,m[cos(ϕ)]`.

# Remarks

This algorithm was based on **[1]**. Our definition of fully normalized associated Legendre
function can be seen in **[2, p. 546]**. The conversion is obtained by:

             ┌                       ┐
             │  (n-m)! . k . (2n+1)  │
    K_n,m = √│ ───────────────────── │,  k = (m = 0) ? 1 : 2.
             │         (n+m)!        │
             └                       ┘

    P̄_n,m = P_n,m * K_n,m,

where `P̄_n,m` is the fully normalized Legendre associated function.

# References

- **[1]** Holmes, S. A. and W. E. Featherstone, 2002. A unified approach to the Clenshaw
    summation and the recursive computation of very high degree and order normalised
    associated Legendre functions. Journal of Geodesy, 76(5), pp. 279-299. For more info.:
    http://mitgcm.org/~mlosch/geoidcookbook/node11.html

- **[2]** Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. Microcosm
    Press, Hawthorn, CA, USA.
"""
function fully_normalized_legendre(
    ϕ::T,
    n_max::Integer,
    m_max::Integer = -1;
    ph_term::Bool = false
) where T<:Number
    n_max < 0 && throw(ArgumentError("n_max must not be negative."))

    if (m_max < 0) || (m_max > n_max)
        m_max = n_max
    end

    P = zeros(float(T), n_max + 1, m_max + 1)
    fully_normalized_legendre!(P, ϕ; ph_term = ph_term)

    return P
end
