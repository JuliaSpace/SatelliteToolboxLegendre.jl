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
    T::Type{<:AbstractFloat} = Float64,
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

    len_P = _packed_length(n_max, m_max_P)
    len_d = _packed_length(n_max, m_max)

    seed = _kernel_seed(N, T)
    diag = _zeros_storage(T, n_max + 1)
    a    = _zeros_storage(T, len_P)
    b    = _zeros_storage(T, len_P)
    da   = _zeros_storage(T, len_d)
    db   = _zeros_storage(T, len_d)

    _fill_legendre_coefficients!(N, diag, a, b, n_max, m_max_P)
    _fill_dlegendre_coefficients!(N, da, db, n_max, m_max)

    return LegendreCoefficients{_normalization(N), T}(
        n_max, m_max, m_max_P, seed, diag, a, b, da, db
    )
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _normalization(N::Val) -> Symbol

Return the symbol of the normalization related to the value `N`.
"""
_normalization(::Val{:full}) = :full

_normalization(::Val{:schmidt}) = :schmidt

_normalization(::Val{:unnormalized}) = :unnormalized

"""
    _fill_legendre_coefficients!(N::Val, diag::_PackedStorage{T}, a::_PackedStorage{T}, b::_PackedStorage{T}, n_max::Integer, m_max::Integer) where T<:AbstractFloat -> Nothing

Fill the packed storages `diag`, `a`, and `b` with the coefficients used to compute the
associated Legendre function with the normalization `N` up to the degree `n_max` and order
`m_max`. The values are precisely those computed on the fly by the kernel in
`src/kernels.jl`, and `a` and `b` are laid out as described in `_packed_index`.
"""
function _fill_legendre_coefficients!(
    N::Union{Val{:full}, Val{:schmidt}, Val{:unnormalized}},
    diag::_PackedStorage{T},
    a::_PackedStorage{T},
    b::_PackedStorage{T},
    n_max::Integer,
    m_max::Integer,
) where {T <: AbstractFloat}
    @inbounds for n in 2:n_max
        aux = _kernel_legendre_aux(N, T, n)

        diag[n + 1] = _kernel_legendre_diag(N, T, n)

        for m in 0:min(n - 1, m_max)
            i = _packed_index(n, m, m_max)
            a[i], b[i] = _kernel_legendre_ab(N, aux, T, n, m)
        end
    end

    return nothing
end

"""
    _fill_dlegendre_coefficients!(N::Val, da::_PackedStorage{T}, db::_PackedStorage{T}, n_max::Integer, m_max::Integer) where T<:AbstractFloat -> Nothing

Fill the packed storages `da` and `db` with the coefficients used to compute the
first-order derivative of the associated Legendre function with the normalization `N` up
to the degree `n_max` and order `m_max`. The values are precisely those computed on the
fly by the kernel in `src/kernels.jl`, and `da` and `db` are laid out as described in
`_packed_index`.
"""
function _fill_dlegendre_coefficients!(
    N::Union{Val{:full}, Val{:schmidt}, Val{:unnormalized}},
    da::_PackedStorage{T},
    db::_PackedStorage{T},
    n_max::Integer,
    m_max::Integer,
) where {T <: AbstractFloat}
    @inbounds for n in 1:n_max
        for m in 0:min(n, m_max)
            i = _packed_index(n, m, m_max)

            da[i] = _kernel_dlegendre_a(N, T, n, m)

            # The coefficient `b_nm` is not used when `n == m`.
            if (m > 0) && (n != m)
                db[i] = _kernel_dlegendre_b(N, T, n, m)
            end
        end
    end

    return nothing
end
