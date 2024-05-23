## Description #############################################################################
#
# Compute the associated Legendre functions without normalization.
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
    unnormalized_legendre!(P::AbstractMatrix{T}, ϕ::Number, n_max::Integer = -1, m_max::Integer = -1; kwargs...) where T<:Number -> Nothing

Compute the unnormalized (or conventional) associated Legendre function `P_n,m[cos(ϕ)]`. The
maximum degree and order that will be computed are given by the parameters `n_max` and
`m_max`. If they are negative (default), the dimensions of matrix `P` will be used:

    maximum degree -> number of rows - 1
    maximum order  -> number of columns - 1

The result will be stored at matrix `P`.

# Keywords

- `ph_term::Bool`: If `true`, the Condon-Shortley phase term `(-1)^m` will be included.
    (**Default** = `false`)
"""
function unnormalized_legendre!(
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
    if s < 0
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

            if n == m
                P[i₀+n, j₀+n] = s_fact * aux_n * P[i₀+n-1, j₀+n-1]
            else
                aux_nm = T(n - m)              # ...................... √((n - m) * (n - m))
                a_nm   = aux_n / aux_nm * c
                b_nm   = T(n + m - 1) / aux_nm # ..... √((n + m - 1) * (n + m - 1)) / aux_nm

                # We assume that the matrix is not initialized. Hence, we must
                # not access elements on the upper triangle.
                if m != n-1
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
    unnormalized_legendre(ϕ::T, n_max::Integer, m_max::Integer = -1; kwargs...) where T<:Number -> Matrix{float(T)}

Compute the unnormalized (or conventional) associated Legendre function `P_n,m[cos(ϕ)]`. The
maximum degree that will be computed is `n_max` and the maximum order is `m_max`. Notice
that if `m_max` is higher than `n_max` or negative, it is set to `n_max`.

# Keywords

- `ph_term::Bool`: If `true`, the Condon-Shortley phase term `(-1)^m` will be included.
    (**Default** = `false`)

# Returns

- `Matrix{float(T)}`: A matrix with the Legendre associated functions `P_n,m[cos(ϕ)]`.
"""
function unnormalized_legendre(
    ϕ::T,
    n_max::Integer,
    m_max::Integer = -1;
    ph_term::Bool = false
) where T<:Number
    n_max < 0 && throw(ArgumentError("n_max must be positive."))

    if (m_max < 0) || (m_max > n_max)
        m_max = n_max
    end

    P = zeros(float(T), n_max + 1, m_max + 1)
    unnormalized_legendre!(P, ϕ; ph_term = ph_term)

    return P
end
