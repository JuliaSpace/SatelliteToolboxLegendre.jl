# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Functions related to the associated Legendre functions.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] Holmes, S. A. and W. E. Featherstone, 2002. A unified approach to the Clenshaw
#       summation and the recursive computation of very high degree and order normalised
#       associated Legendre functions. Journal of Geodesy, 76(5), pp. 279-299.
#
#       For more info.: http://mitgcm.org/~mlosch/geoidcookbook/node11.html
#
#   [2] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. Microcosm
#       Press, Hawthorn, CA, USA.
#
#   [3] Schmidt, A (1917). Erdmagnetismus, Enzykl. Math. Wiss., 6, pp. 265–396.
#
#   [4] Winch, D. E., Ivers, D. J., Turner, J. P. R., Stening R. J (2005). Geomagnetism and
#       Schmidt quasi-normalization. Geophysical Journal International, 160(2), pp. 487-504.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export legendre!, legendre

"""
    legendre!(N, P::AbstractMatrix, ϕ::Number, n_max::Integer = -1, m_max::Integer = -1; kwargs...) -> Nothing

Compute the associated Legendre function `P_n,m[cos(ϕ)]`. The maximum degree and order that
will be computed are given by the parameters `n_max` and `m_max`. If they are negative
(default), the dimensions of matrix `P` will be used:

    maximum degree -> number of rows - 1
    maximum order  -> number of columns - 1

The result will be stored at matrix `P`.

The parameter `N` selects the normalization. The following values are valid:

- `Val(:full)`: Compute the fully normalized associated Legendre function.
- `Val(:schmidt)`: Compute the Schmidt quasi-normalized associated Legendre function.
- `Val(:unnormalized)`: Compute the unnormalized associated Legendre function.

# Keywords

- `ph_term::Bool`: If `true`, then the Condon-Shortley phase term `(-1)^m` will be included.
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
    ::Val{:full},
    P::AbstractMatrix,
    ϕ::Number,
    n_max::Integer = -1,
    m_max::Integer = -1,
    ph_term::Bool = false
)
    return fully_normalized_legendre!(P, ϕ, n_max, m_max; ph_term = ph_term)
end

function legendre!(
    ::Val{:schmidt},
    P::AbstractMatrix,
    ϕ::Number,
    n_max::Integer = -1,
    m_max::Integer = -1;
    ph_term::Bool = false
)
    return schmidt_quasi_normalized_legendre!(P, ϕ, n_max, m_max; ph_term = ph_term)
end

function legendre!(
    ::Val{:unnormalized},
    P::AbstractMatrix,
    ϕ::Number,
    n_max::Integer = -1,
    m_max::Integer = -1;
    ph_term::Bool = false
)
    return unnormalized_legendre!(P, ϕ, n_max, m_max; ph_term = ph_term)
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

- `ph_term::Bool`: If `true`, then the Condon-Shortley phase term `(-1)^m` will be included.
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
    ::Val{:full},
    ϕ::Number,
    n_max::Integer,
    m_max::Integer = -1;
    ph_term::Bool = false
)
    return fully_normalized_legendre(ϕ, n_max, m_max; ph_term = ph_term)
end

function legendre(
    ::Val{:schmidt},
    ϕ::Number,
    n_max::Integer,
    m_max::Integer = -1;
    ph_term::Bool = false
)
    return schmidt_quasi_normalized_legendre(ϕ, n_max, m_max; ph_term = ph_term)
end

function legendre(
    ::Val{:unnormalized},
    ϕ::Number,
    n_max::Integer,
    m_max::Integer = -1;
    ph_term::Bool = false
)
    return unnormalized_legendre(ϕ, n_max, m_max; ph_term = ph_term)
end
