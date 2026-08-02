## Description #############################################################################
#
# Precomputed coefficients to accelerate the computation of the associated Legendre
#   functions and their first-order derivatives.
#
############################################################################################

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
