## Description #############################################################################
#
# Precomputed coefficients to accelerate the computation of the associated Legendre
#   functions and their first-order derivatives.
#
############################################################################################

export LegendreCoefficients

"""
    struct LegendreCoefficients{N, T<:AbstractFloat}

Precomputed recursion coefficients to compute the associated Legendre function
`P_n,m[cos(ϕ)]` and its first-order derivative with respect to `ϕ` using the normalization
`N` (`:full`, `:schmidt`, or `:unnormalized`) and the element type `T`.

The coefficients depend only on the degree and order. Hence, when the functions are
evaluated at many angles `ϕ` with the same maximum degree and order, precomputing them
avoids evaluating square roots at every call, largely improving the performance.

# Fields

- `n_max::Int`: Maximum degree supported by the coefficients.
- `m_max::Int`: Maximum order supported by the derivative coefficients.
- `m_max_P::Int`: Maximum order supported by the Legendre function coefficients. It is one
    order higher than `m_max` (clamped to `n_max`) so the matrix `P` required by
    [`dlegendre!`](@ref) can be computed with the same object.
- `seed::T`: Coefficient of the terms with degree 1.
- `diag::Vector{T}`: Coefficients of the diagonal recursion (`n == m`).
- `a::Matrix{T}`: Coefficients `a_nm` of the Legendre function recursion.
- `b::Matrix{T}`: Coefficients `b_nm` of the Legendre function recursion.
- `da::Matrix{T}`: Coefficients `a_nm` of the derivative equation.
- `db::Matrix{T}`: Coefficients `b_nm` of the derivative equation.

# See Also

[`legendre!`](@ref), [`dlegendre!`](@ref)
"""
struct LegendreCoefficients{N, T<:AbstractFloat}
    n_max::Int
    m_max::Int
    m_max_P::Int
    seed::T
    diag::Vector{T}
    a::Matrix{T}
    b::Matrix{T}
    da::Matrix{T}
    db::Matrix{T}
end

############################################################################################
#                                       Constructor                                        #
############################################################################################

"""
    LegendreCoefficients(N, n_max::Integer, m_max::Integer = -1; T::Type{<:AbstractFloat} = Float64) -> LegendreCoefficients

Create the precomputed coefficients to compute the associated Legendre function
`P_n,m[cos(ϕ)]` and its first-order derivative with respect to `ϕ` up to the maximum degree
`n_max` and maximum order `m_max` using the element type `T`. Notice that if `m_max` is
higher than `n_max` or negative (default), it is set to `n_max`.

The parameter `N` selects the normalization. The following values are valid:

- `Val(:full)`: Compute the fully normalized associated Legendre function.
- `Val(:schmidt)`: Compute the Schmidt quasi-normalized associated Legendre function.
- `Val(:unnormalized)`: Compute the unnormalized associated Legendre function.

The returned object can be passed to [`legendre!`](@ref) and [`dlegendre!`](@ref) to
compute the functions without evaluating square roots at every call.

# Keywords

- `T::Type{<:AbstractFloat}`: Element type of the coefficients. (**Default** = `Float64`)

# Examples

```julia-repl
julia> coefs = LegendreCoefficients(Val(:full), 60);

julia> P = zeros(61, 61);

julia> legendre!(P, 0.123, coefs)

julia> dP = zeros(61, 61);

julia> dlegendre!(dP, 0.123, P, coefs)
```
"""
function LegendreCoefficients(
    N::Union{Val{:full}, Val{:schmidt}, Val{:unnormalized}},
    n_max::Integer,
    m_max::Integer = -1;
    T::Type{<:AbstractFloat} = Float64
)
    n_max = Int(n_max)
    m_max = Int(m_max)

    n_max < 0 && throw(ArgumentError("n_max must not be negative."))

    if (m_max < 0) || (m_max > n_max)
        m_max = n_max
    end

    # The Legendre function coefficients support one additional order so the matrix `P`
    # required by the derivative can be computed with the same object.
    m_max_P = min(m_max + 1, n_max)

    seed = _legendre_seed(N, T)
    diag = zeros(T, n_max + 1)
    a    = zeros(T, n_max + 1, m_max_P + 1)
    b    = zeros(T, n_max + 1, m_max_P + 1)
    da   = zeros(T, n_max + 1, m_max + 1)
    db   = zeros(T, n_max + 1, m_max + 1)

    _fill_legendre_coefficients!(N, diag, a, b, n_max, m_max_P)
    _fill_dlegendre_coefficients!(N, da, db, n_max, m_max)

    return LegendreCoefficients{_normalization(N), T}(
        n_max,
        m_max,
        m_max_P,
        seed,
        diag,
        a,
        b,
        da,
        db
    )
end

############################################################################################
#                                    Public Functions                                      #
############################################################################################

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
    ph_term::Bool = false
) where {N, T<:AbstractFloat}
    # Obtain the maximum degree and order that must be computed.
    n_max, m_max = _get_degree_and_order(P, n_max, m_max)

    if (n_max > coefs.n_max) || (m_max > coefs.m_max_P)
        throw(ArgumentError(
            "The coefficients support the maximum degree $(coefs.n_max) and order " *
            "$(coefs.m_max_P), but the computation requires degree $n_max and order " *
            "$m_max."
        ))
    end

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
            P[i₀ + 1, j₀] = coefs.seed * c

            if m_max > 0
                P[i₀ + 1, j₀ + 1] = coefs.seed * s_fact
            end

            continue
        end

        for m in 0:n
            P_nm = zero(T)

            if n == m
                P_nm = s_fact * coefs.diag[n + 1] * T(P[i₀ + n - 1, j₀ + n - 1])

            else
                a_nm = coefs.a[n + 1, m + 1] * c

                # We assume that the matrix is not initialized. Hence, we must not access
                # elements on the upper triangle.
                if m != n - 1
                    P_nm = a_nm * T(P[i₀ + n - 1, j₀ + m]) -
                        coefs.b[n + 1, m + 1] * T(P[i₀ + n - 2, j₀ + m])
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

"""
    dlegendre!(dP::AbstractMatrix, ϕ::Number, P::AbstractMatrix, coefs::LegendreCoefficients{N, T}, n_max::Integer = -1, m_max::Integer = -1; kwargs...) where {N, T<:AbstractFloat} -> Nothing

Compute the first-order derivative of the associated Legendre function `P_n,m[cos(ϕ)]`
with respect to `ϕ` [rad] using the precomputed coefficients `coefs`, which also select the
normalization (see [`LegendreCoefficients`](@ref)). The maximum degree and order that will
be computed are given by the parameters `n_max` and `m_max`. If they are negative
(default), the dimensions of matrix `dP` will be used. This function throws an
`ArgumentError` if the maximum degree or order exceeds the ones supported by `coefs`.

The derivatives will be stored in the matrix `dP` and all the arithmetic operations are
performed using the element type `T` of the coefficients.

This algorithm needs the matrix `P` with the values of the associated Legendre function
using the same normalization, which can be computed using [`legendre!`](@ref) with the same
coefficients `coefs`. The algorithm accesses the terms `P[n, m + 1]` when the order `m` is
lower than the degree `n`. Hence, the matrix `P` must have at least `n_max + 1` rows and
`m_max + 2` columns if `m_max < n_max`, or `m_max + 1` columns otherwise. This function
throws an `ArgumentError` if `P` does not have the required dimensions.

!!! note

    This function only writes to the elements of `dP` in the lower triangular part up to
    the computed degree and order. If `dP` is being reused with a smaller `n_max` or
    `m_max`, the remaining elements will keep their previous (stale) values.

!!! warning

    The user is responsible for passing a matrix `P` with the correct values. For example,
    if `ph_term` is `true`, `P` must also be computed with `ph_term` set to `true`.

# Keywords

- `ph_term::Bool`: If `true`, the Condon-Shortley phase term `(-1)^m` will be included.
    (**Default** = `false`)
"""
function dlegendre!(
    dP::AbstractMatrix,
    ϕ::Number,
    P::AbstractMatrix,
    coefs::LegendreCoefficients{N, T},
    n_max::Integer = -1,
    m_max::Integer = -1;
    ph_term::Bool = false
) where {N, T<:AbstractFloat}
    # Obtain the maximum degree and order that must be computed.
    n_max, m_max = _get_degree_and_order(dP, P, n_max, m_max)

    if (n_max > coefs.n_max) || (m_max > coefs.m_max)
        throw(ArgumentError(
            "The coefficients support the maximum degree $(coefs.n_max) and order " *
            "$(coefs.m_max), but the computation requires degree $n_max and order $m_max."
        ))
    end

    # See the kernels in `src/dlegendre/` for the definition of the variable `fact`.
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
                dP_nm = -coefs.da[n + 1, 1] * T(P[i₀ + n, j₀ + 1])
            else
                dP_nm = coefs.da[n + 1, m + 1] * T(P[i₀ + n, j₀ + m - 1])

                if n != m
                    dP_nm -= coefs.db[n + 1, m + 1] * T(P[i₀ + n, j₀ + m + 1])
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
#                                    Private Functions                                     #
############################################################################################

# Return the symbol of the normalization related to the value `N`.
_normalization(::Val{:full})         = :full
_normalization(::Val{:schmidt})      = :schmidt
_normalization(::Val{:unnormalized}) = :unnormalized

# Return the coefficient of the terms with degree 1 for the normalization `N` using the
# element type `T`.
_legendre_seed(::Val{:full}, T::Type)         = √T(3)
_legendre_seed(::Val{:schmidt}, T::Type)      = one(T)
_legendre_seed(::Val{:unnormalized}, T::Type) = one(T)

# Fill the coefficients used to compute the associated Legendre function. The values are
# precisely those computed by the kernels in `src/legendre/`.
function _fill_legendre_coefficients!(
    ::Val{:full},
    diag::Vector{T},
    a::Matrix{T},
    b::Matrix{T},
    n_max::Integer,
    m_max::Integer
) where T<:AbstractFloat
    @inbounds for n in 2:n_max
        sq_2n_p_1 = √T(2n + 1)
        aux_an    = √T(2n - 1) * sq_2n_p_1 # ..................... √((2n - 1) * (2n + 1))
        aux_bn    = sq_2n_p_1 / √T(2n - 3) # ....................... √((2n + 1) / (2n - 3))

        diag[n + 1] = √(T(2n + 1) / T(2n))

        for m in 0:min(n - 1, m_max)
            aux_nm       = √(T(n - m) * T(n + m))
            a[n + 1, m + 1] = aux_an / aux_nm
            b[n + 1, m + 1] = √(T(n + m - 1) * T(n - m - 1)) * aux_bn / aux_nm
        end
    end

    return nothing
end

function _fill_legendre_coefficients!(
    ::Val{:schmidt},
    diag::Vector{T},
    a::Matrix{T},
    b::Matrix{T},
    n_max::Integer,
    m_max::Integer
) where T<:AbstractFloat
    @inbounds for n in 2:n_max
        aux_n = T(2n - 1)

        diag[n + 1] = √(aux_n / T(2n))

        for m in 0:min(n - 1, m_max)
            aux_nm       = √(T(n - m) * T(n + m))
            a[n + 1, m + 1] = aux_n / aux_nm
            b[n + 1, m + 1] = √(T(n + m - 1) * T(n - m - 1)) / aux_nm
        end
    end

    return nothing
end

function _fill_legendre_coefficients!(
    ::Val{:unnormalized},
    diag::Vector{T},
    a::Matrix{T},
    b::Matrix{T},
    n_max::Integer,
    m_max::Integer
) where T<:AbstractFloat
    @inbounds for n in 2:n_max
        aux_n = T(2n - 1)

        diag[n + 1] = aux_n

        for m in 0:min(n - 1, m_max)
            aux_nm       = T(n - m)
            a[n + 1, m + 1] = aux_n / aux_nm
            b[n + 1, m + 1] = T(n + m - 1) / aux_nm
        end
    end

    return nothing
end

# Fill the coefficients used to compute the first-order derivative of the associated
# Legendre function. The values are precisely those computed by the kernels in
# `src/dlegendre/`.
function _fill_dlegendre_coefficients!(
    ::Union{Val{:full}, Val{:schmidt}},
    da::Matrix{T},
    db::Matrix{T},
    n_max::Integer,
    m_max::Integer
) where T<:AbstractFloat
    @inbounds for n in 1:n_max
        for m in 0:min(n, m_max)
            if m == 0
                da[n + 1, 1] = √(T(n) * T(n + 1) / 2)
            elseif m == 1
                da[n + 1, 2] = √(T(2n) * T(n + 1)) / 2
                (n > 1) && (db[n + 1, 2] = √(T(n + 2) * T(n - 1)) / 2)
            else
                da[n + 1, m + 1] = √(T(n + m) * T(n - m + 1)) / 2
                (n != m) && (db[n + 1, m + 1] = √(T(n + m + 1) * T(n - m)) / 2)
            end
        end
    end

    return nothing
end

function _fill_dlegendre_coefficients!(
    ::Val{:unnormalized},
    da::Matrix{T},
    db::Matrix{T},
    n_max::Integer,
    m_max::Integer
) where T<:AbstractFloat
    @inbounds for n in 1:n_max
        for m in 0:min(n, m_max)
            if m == 0
                da[n + 1, 1] = 1
            else
                da[n + 1, m + 1] = T(n + m) * T(n - m + 1) / 2
                (n != m) && (db[n + 1, m + 1] = T(1) / 2)
            end
        end
    end

    return nothing
end
