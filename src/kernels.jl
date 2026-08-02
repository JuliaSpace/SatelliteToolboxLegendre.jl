## Description #############################################################################
#
# Generic kernels to compute the associated Legendre functions and their first-order
#   derivatives.
#
# The kernels are parameterized by a coefficient provider `coefs`, which can be:
#
#   1. A `Val` selecting the normalization (`Val(:full)`, `Val(:schmidt)`, or
#      `Val(:unnormalized)`), leading to coefficients computed on the fly; or
#   2. A `LegendreCoefficients` object, leading to coefficients read from precomputed
#      arrays.
#
# All the arithmetic operations are performed using the element type `T` passed to the
# kernels.
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
# [5] Du, J., Chen, C., Lesur, V., and Wang, L (2015). Non-singular spherical harmonic
#     expressions of geomagnetic vector and gradient tensor fields in the local
#     north-oriented reference frame. Geoscientific Model Development, 8, pp. 1979-1990.
#
# [6] Ilk, K. H.: Ein Beitrag zur Dynamik ausgedehnter Körper-Gravitationswechselwirkung,
#     Deutsche Geodätische Kommission.  Reihe C, Heft Nr. 288, München, 1983.
#
############################################################################################

############################################################################################
#                                     Legendre Kernel                                      #
############################################################################################

# Compute the associated Legendre function `P_n,m[cos(ϕ)]` up to the degree `n_max` and
# order `m_max` using the coefficient provider `coefs`, storing the result in `P`. All the
# arithmetic operations are performed using the type `T`.
function _legendre_kernel!(
    P::AbstractMatrix,
    ϕ::Number,
    coefs,
    n_max::Int,
    m_max::Int,
    ph_term::Bool,
    ::Type{T}
) where T
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

    seed = _kernel_seed(coefs, T)

    @inbounds for n in 0:n_max
        # Starting values.
        if n == 0
            P[i₀, j₀] = 1
            continue

        elseif n == 1
            P[i₀ + 1, j₀] = seed * c

            if m_max > 0
                P[i₀ + 1, j₀ + 1] = seed * s_fact
            end

            continue
        end

        # Auxiliary terms that depend only on the degree `n`.
        aux = _kernel_legendre_aux(coefs, T, n)

        for m in 0:n
            P_nm = zero(T)

            if n == m
                P_nm = s_fact * _kernel_legendre_diag(coefs, T, n) *
                    T(P[i₀ + n - 1, j₀ + n - 1])

            else
                a_nm, b_nm = _kernel_legendre_ab(coefs, aux, T, n, m)
                a_nm *= c

                # We assume that the matrix is not initialized. Hence, we must not access
                # elements on the upper triangle.
                if m != n - 1
                    P_nm = a_nm * T(P[i₀ + n - 1, j₀ + m]) -
                        b_nm * T(P[i₀ + n - 2, j₀ + m])
                else
                    P_nm = a_nm * T(P[i₀ + n - 1, j₀ + m])
                end
            end

            P[i₀ + n, j₀ + m] = P_nm

            # Check if the maximum desired order has been reached.
            m >= m_max && break
        end
    end

    return nothing
end

############################################################################################
#                                    Derivative Kernel                                     #
############################################################################################

# Compute the first-order derivative of the associated Legendre function `P_n,m[cos(ϕ)]`
# with respect to `ϕ` up to the degree `n_max` and order `m_max` using the coefficient
# provider `coefs`, storing the result in `dP`. The matrix `P` must contain the values of
# the associated Legendre function with the same normalization. All the arithmetic
# operations are performed using the type `T`.
function _dlegendre_kernel!(
    dP::AbstractMatrix,
    ϕ::Number,
    P::AbstractMatrix,
    coefs,
    n_max::Int,
    m_max::Int,
    ph_term::Bool,
    ::Type{T}
) where T
    # The derivative is computed using the following equation [5, p. 1981]:
    #
    #   ∂P(n, m)
    #   ──────── = a_nm . P(n,m-1) - b_nm . P(n,m+1),
    #      ∂ϕ
    #
    # where the coefficients `a_nm` and `b_nm` depend on the normalization. When `m = 0`,
    # the equation uses the terms `P(n, -1)` and `P(n, +1)`. Notice that [5, p. 1985]:
    #
    #                   m
    #   P_(n, -m) = (-1)  . P_(n, m).
    #
    # Since both coefficients are equal in this case, the equation simplifies to a single
    # product.
    #
    # NOTE: This algorithm is based on eq. Z.1.44 of [6], which is valid only for
    # ϕ ∈ [0, π] (see [6, p. 119]). In this package, the Legendre associated functions are
    # computed using the convention that `sin(ϕ)` is always positive. Under this
    # convention, the computed function is even about ϕ = 0 and ϕ = π. Hence, for
    # ϕ ∈ (π, 2π), the derivative equals the negative of the value obtained from the
    # coefficients, which is applied here using the variable `fact`. This behavior was
    # verified numerically against finite differences of the values returned by the
    # corresponding Legendre function for the entire circle, including angles outside
    # [0, 2π]. At the points ϕ ∈ {0, π, 2π}, where the convention renders the computed
    # function non-differentiable, this function returns the one-sided derivative (from
    # the right at 0 and 2π, and from the left at π).

    ϕ    = mod(ϕ, T(2π))
    fact = ϕ > T(π) ? -1 : 1

    if ph_term
        fact *= -1
    end

    # Get the first indices in `P` to take into account offset arrays.
    i₀, j₀ = first.(axes(P))

    # Get the first indices in `dP` to take into account offset arrays.
    di₀, dj₀ = first.(axes(dP))

    dP[di₀, dj₀] = 0

    @inbounds for n in 1:n_max
        for m in 0:n
            dP_nm = zero(T)

            if m == 0
                dP_nm = -_kernel_dlegendre_a(coefs, T, n, 0) * T(P[i₀ + n, j₀ + 1])
            else
                dP_nm = _kernel_dlegendre_a(coefs, T, n, m) * T(P[i₀ + n, j₀ + m - 1])

                if n != m
                    dP_nm -= _kernel_dlegendre_b(coefs, T, n, m) * T(P[i₀ + n, j₀ + m + 1])
                end
            end

            dP[di₀ + n, dj₀ + m] = fact * dP_nm

            # Check if the maximum desired order has been reached.
            m >= m_max && break
        end
    end

    return nothing
end

############################################################################################
#                                  Coefficient Providers                                   #
############################################################################################

# == Seed ==================================================================================

# Return the coefficient of the terms with degree 1.
@inline _kernel_seed(::Val{:full}, ::Type{T}) where T         = √T(3)
@inline _kernel_seed(::Val{:schmidt}, ::Type{T}) where T      = one(T)
@inline _kernel_seed(::Val{:unnormalized}, ::Type{T}) where T = one(T)
@inline _kernel_seed(coefs::LegendreCoefficients, ::Type{T}) where T = coefs.seed

# == Legendre Function Coefficients ========================================================

# Return the auxiliary terms of the Legendre function recursion that depend only on the
# degree `n`.
@inline function _kernel_legendre_aux(::Val{:full}, ::Type{T}, n::Int) where T
    # Compute the square roots of the terms that depend only on `n` outside the inner
    # loop. Using ratios of already-rooted factors avoids products that can overflow or
    # lose precision at very high degrees.
    sq_2n_p_1 = √T(2n + 1)
    return (
        √T(2n - 1) * sq_2n_p_1, # ................................ √((2n - 1) * (2n + 1))
        sq_2n_p_1 / √T(2n - 3)  # .................................. √((2n + 1) / (2n - 3))
    )
end

@inline _kernel_legendre_aux(::Val{:schmidt}, ::Type{T}, n::Int) where T      = (T(2n - 1),)
@inline _kernel_legendre_aux(::Val{:unnormalized}, ::Type{T}, n::Int) where T = (T(2n - 1),)
@inline _kernel_legendre_aux(::LegendreCoefficients, ::Type{T}, n::Int) where T = nothing

# Return the coefficient of the diagonal recursion (`n == m`) of the Legendre function.
@inline _kernel_legendre_diag(::Val{:full}, ::Type{T}, n::Int) where T    = √(T(2n + 1) / T(2n))
@inline _kernel_legendre_diag(::Val{:schmidt}, ::Type{T}, n::Int) where T = √(T(2n - 1) / T(2n))
@inline _kernel_legendre_diag(::Val{:unnormalized}, ::Type{T}, n::Int) where T = T(2n - 1)
@inline function _kernel_legendre_diag(coefs::LegendreCoefficients, ::Type{T}, n::Int) where T
    return @inbounds coefs.diag[n + 1]
end

# Return the coefficients `a_nm` and `b_nm` of the Legendre function recursion.
@inline function _kernel_legendre_ab(::Val{:full}, aux, ::Type{T}, n::Int, m::Int) where T
    aux_nm = √(T(n - m) * T(n + m))
    a_nm   = aux[1] / aux_nm
    b_nm   = √(T(n + m - 1) * T(n - m - 1)) * aux[2] / aux_nm
    return a_nm, b_nm
end

@inline function _kernel_legendre_ab(::Val{:schmidt}, aux, ::Type{T}, n::Int, m::Int) where T
    aux_nm = √(T(n - m) * T(n + m))
    a_nm   = aux[1] / aux_nm
    b_nm   = √(T(n + m - 1) * T(n - m - 1)) / aux_nm
    return a_nm, b_nm
end

@inline function _kernel_legendre_ab(::Val{:unnormalized}, aux, ::Type{T}, n::Int, m::Int) where T
    aux_nm = T(n - m)              # .......................... √((n - m) * (n - m))
    a_nm   = aux[1] / aux_nm
    b_nm   = T(n + m - 1) / aux_nm # ......... √((n + m - 1) * (n + m - 1)) / aux_nm
    return a_nm, b_nm
end

@inline function _kernel_legendre_ab(
    coefs::LegendreCoefficients,
    aux,
    ::Type{T},
    n::Int,
    m::Int
) where T
    return @inbounds coefs.a[n + 1, m + 1], coefs.b[n + 1, m + 1]
end

# == Derivative Coefficients ===============================================================

# Return the coefficient `a_nm` of the derivative equation.
#
# The coefficients for the full normalization are [5, p. 1981]:
#
#   a_nm = ¹/₂ . √(n + m) . √(n - m + 1) . √(C_{m} / C_{m-1}),
#
#   b_nm = ¹/₂ . √(n + m + 1) . √(n - m) . √(C_{m} / C_{m+1}),
#
#           ┌
#           │ 1, m  = 0
#   C_{m} = │
#           │ 2, m != 0
#           └
#
# NOTE: The conversion of the coefficients between the full normalization and the Schmidt
# quasi-normalization is performed using:
#
#   √(2n + 1),
#
# which depends only on `n`. Since the derivative equation only has terms related to the
# order `n`, the same coefficients work for both normalizations.
@inline function _kernel_dlegendre_a(
    ::Union{Val{:full}, Val{:schmidt}},
    ::Type{T},
    n::Int,
    m::Int
) where T
    (m == 0) && return √(T(n) * T(n + 1) / 2)

    # We should consider the case `m == 1` separately from the general one because of the
    # coefficient `C_{m}`.
    (m == 1) && return √(T(2n) * T(n + 1)) / 2

    return √(T(n + m) * T(n - m + 1)) / 2
end

@inline function _kernel_dlegendre_a(::Val{:unnormalized}, ::Type{T}, n::Int, m::Int) where T
    (m == 0) && return one(T)
    return T(n + m) * T(n - m + 1) / 2
end

@inline function _kernel_dlegendre_a(coefs::LegendreCoefficients, ::Type{T}, n::Int, m::Int) where T
    return @inbounds coefs.da[n + 1, m + 1]
end

# Return the coefficient `b_nm` of the derivative equation.
@inline function _kernel_dlegendre_b(
    ::Union{Val{:full}, Val{:schmidt}},
    ::Type{T},
    n::Int,
    m::Int
) where T
    # We should consider the case `m == 1` separately from the general one because of the
    # coefficient `C_{m}`.
    (m == 1) && return √(T(n + 2) * T(n - 1)) / 2

    return √(T(n + m + 1) * T(n - m)) / 2
end

@inline _kernel_dlegendre_b(::Val{:unnormalized}, ::Type{T}, n::Int, m::Int) where T = T(1) / 2

@inline function _kernel_dlegendre_b(coefs::LegendreCoefficients, ::Type{T}, n::Int, m::Int) where T
    return @inbounds coefs.db[n + 1, m + 1]
end
