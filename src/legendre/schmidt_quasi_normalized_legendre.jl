## Description #############################################################################
#
# Compute the associated Legendre functions with the Schmidt quasi-normalization.
#
## References ##############################################################################
#
# [1] Schmidt, A (1917). Erdmagnetismus, Enzykl. Math. Wiss., 6, pp. 265–396.
#
# [2] Winch, D. E., Ivers, D. J., Turner, J. P. R., Stening R. J (2005). Geomagnetism and
#     Schmidt quasi-normalization. Geophysical Journal International, 160(2), pp. 487-504.
#
############################################################################################

"""
    schmidt_quasi_normalized_legendre!(P::AbstractMatrix{T}, ϕ::Number, n_max::Integer = -1, m_max::Integer = -1; kwargs...) where T<:Number -> Nothing

Compute the Schmidt quasi-normalized associated Legendre function `P_n,m[cos(ϕ)]`
**[1, 2]**. The maximum degree and order that will be computed are given by the parameters
`n_max` and `m_max`. If they are negative (default), the dimensions of matrix `P` will be
used:

    maximum degree -> number of rows - 1
    maximum order  -> number of columns - 1

The result will be stored in matrix `P`.

# Keywords

- `ph_term::Bool`: If `true`, the Condon-Shortley phase term `(-1)^m` will be included.
    (**Default** = `false`)

# Remarks

This algorithm was based on **[1, 2]**. The conversion is obtained by:

             ┌              ┐
             │      (n-m)!  │
    K_n,m = √│ k . ──────── │,  k = (m = 0) ? 1 : 2.
             │      (n+m)!  │
             └              ┘

    P̂_n,m = P_n,m * K_n,m,

where `P̂_n,m` is the Schmidt quasi-normalized Legendre associated function.

# References

- **[1]** Schmidt, A (1917). Erdmagnetismus, Enzykl. Math. Wiss., 6, pp. 265–396.

- **[2]** Winch, D. E., Ivers, D. J., Turner, J. P. R., Stening R. J (2005). Geomagnetism
    and Schmidt quasi-normalization. Geophysical Journal International, 160(2), pp. 487-504.
"""
function schmidt_quasi_normalized_legendre!(
    P::AbstractMatrix{T},
    ϕ::Number,
    n_max::Integer = -1,
    m_max::Integer = -1;
    ph_term::Bool = false
) where T<:Number
    # Obtain the maximum degree and order that must be computed.
    n_max, m_max = _get_degree_and_order(P, n_max, m_max)

    return _legendre_kernel!(P, ϕ, Val(:schmidt), n_max, m_max, ph_term, T)
end

"""
    schmidt_quasi_normalized_legendre(ϕ::T, n_max::Integer, m_max::Integer = -1; kwargs...) where T<:Number -> Matrix{float(T)}

Compute the Schmidt quasi-normalized associated Legendre function `P_n,m[cos(ϕ)]`. The
maximum degree that will be computed is `n_max` and the maximum order is `m_max`. Notice
that if `m_max` is higher than `n_max` or negative, it is set to `n_max`.

# Keywords

- `ph_term::Bool`: If `true`, the Condon-Shortley phase term `(-1)^m` will be included.
    (**Default** = `false`)

# Returns

- `Matrix{float(T)}`: A matrix with the Legendre associated functions `P_n,m[cos(ϕ)]`.

# Remarks

This algorithm was based on **[1, 2]**. The conversion is obtained by:

             ┌              ┐
             │      (n-m)!  │
    K_n,m = √│ k . ──────── │,  k = (m = 0) ? 1 : 2.
             │      (n+m)!  │
             └              ┘

    P̂_n,m = P_n,m * K_n,m,

where `P̂_n,m` is the Schmidt quasi-normalized Legendre associated function.

# References

- **[1]** Schmidt, A (1917). Erdmagnetismus, Enzykl. Math. Wiss., 6, pp. 265–396.

- **[2]** Winch, D. E., Ivers, D. J., Turner, J. P. R., Stening R. J (2005). Geomagnetism
    and Schmidt quasi-normalization. Geophysical Journal International, 160(2), pp. 487-504.
"""
function schmidt_quasi_normalized_legendre(
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
    schmidt_quasi_normalized_legendre!(P, ϕ; ph_term = ph_term)

    return P
end
