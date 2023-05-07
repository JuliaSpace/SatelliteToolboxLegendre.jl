# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Compute the associated Legendre functions with the Schmidt quasi-normalization.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] Schmidt, A (1917). Erdmagnetismus, Enzykl. Math. Wiss., 6, pp. 265–396.
#
#   [2] Winch, D. E., Ivers, D. J., Turner, J. P. R., Stening R. J (2005).
#       Geomagnetism and Schmidt quasi-normalization. Geophysical Journal
#       International, 160(2), pp. 487-504.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

"""
    schmidt_quasi_normalized_legendre!(P::AbstractMatrix{T}, ϕ::Number, n_max::Integer = -1, m_max::Integer = -1; kwargs...) where T<:Number -> Nothing

Compute the Schmidt quasi-normalized associated Legendre function `P_n,m[cos(ϕ)]`
**[1, 2]**. The maximum degree and order that will be computed are given by the parameters
`n_max` and `m_max`. If they are negative (default), the dimensions of matrix `P` will be
used:

    maximum degree -> number of rows - 1
    maximum order  -> number of columns - 1

The result will be stored at matrix `P`.

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

    # Auxiliary variables to improve code performance.
    s, c = sincos(T(ϕ))

    # The sine must be always positive. In fact, `s` was previously computed using
    # `√(1 - c^2)`. However, we had numerical problems for very small angles that lead to
    # `cos(ϕ) = 1`.
    if (s < 0)
        s = -s
    end

    s_fact = !ph_term ? +s : -s

    # Get the first indices in `P` to take into account offset arrays.
    i₀, j₀ = first.(axes(P))

    @inbounds for n in 0:n_max
        # Starting values.
        if n == 0
            P[i₀, j₀] = 1
            continue

        elseif n == 1
            P[i₀+1, j₀] = +c

            if m_max > 0
                P[i₀+1, j₀+1] = +s_fact
            end

            continue
        end

        aux_n = T(2n - 1) # ......................................... √((2n - 1) * (2n - 1))

        for m in 0:n

            if m == n
                P[i₀+n, j₀+n] = s_fact * √(aux_n / T(2n)) * P[i₀+n-1, j₀+n-1]

            else
                aux_nm = √(T(n - m) * T(n + m))
                a_nm   = aux_n / aux_nm * c
                b_nm   = √(T(n + m - 1) * T(n - m - 1)) / aux_nm

                # We assume that the matrix is not initialized. Hence, we must not access
                # elements on the upper triangle.
                if m != n - 1
                    P[i₀+n, j₀+m] = a_nm * P[i₀+n-1, j₀+m] - b_nm * P[i₀+n-2, j₀+m]
                else
                    P[i₀+n, j₀+m] = a_nm * P[i₀+n-1, j₀+m]
                end
            end

            # Check if the maximum desired order has been reached.
            m == m_max && break
        end
    end

    return nothing
end

"""
    schmidt_quasi_normalized_legendre(ϕ::T, n_max::Integer = -1, m_max::Integer = -1; kwargs...) where T<:Number -> Matrix{float(T)}

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
    n_max::Integer = -1,
    m_max::Integer = -1;
    ph_term::Bool = false
) where T<:Number
    n_max < 0 && throw(ArgumentError("n_max must be positive."))

    if (m_max < 0) || (m_max > n_max)
        m_max = n_max
    end

    P = zeros(float(T), n_max + 1, m_max + 1)
    schmidt_quasi_normalized_legendre!(P, ϕ; ph_term = ph_term)

    return P
end
