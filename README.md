SatelliteToolboxLegendre.jl
===========================

[![CI](https://github.com/JuliaSpace/SatelliteToolboxLegendre.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JuliaSpace/SatelliteToolboxLegendre.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/JuliaSpace/SatelliteToolboxLegendre.jl/branch/main/graph/badge.svg?token=AUE8ZZ5IXJ)](https://codecov.io/gh/JuliaSpace/SatelliteToolboxLegendre.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

This package contains function to compute the Legendre associated functions and its
time-derivatives for the models in the **SatelliteToolbox.jl** ecosystem.

## Installation

```julia
julia> using Pkg
julia> Pkg.install("SatelliteToolboxLegendre")
```

## Usage

### Legendre Associated Functions

This package exports two methods to compute the Legendre associated functions: `legendre`
and `legendre!`.

    legendre(N, ϕ::T, n_max::Integer, m_max::Integer = -1; ph_term::Bool = false) where T<:Number -> Matrix{float(T)}

Compute the associated Legendre function `P_n,m[cos(ϕ)]`. The maximum degree that will be
computed is `n_max` and the maximum order is `m_max`. Notice that if `m_max` is higher than
`n_max` or negative (default), it is set to `n_max`.

The parameter `N` selects the normalization. The following values are valid:

- `Val(:full)`: Compute the fully normalized associated Legendre function.
- `Val(:schmidt)`: Compute the Schmidt quasi-normalized associated Legendre function.
- `Val(:unnormalized)`: Compute the unnormalized associated Legendre function.

This function has the following keywords:

- `ph_term::Bool`: If `true`, then the Condon-Shortley phase term `(-1)^m` will be included.
    (**Default** = `false`)

It returns a `Matrix{float(T)}` with the Legendre associated functions `P_n,m[cos(ϕ)]`.

```julia
julia> legendre(Val(:unnormalized), 0.45, 4)
5×5 Matrix{Float64}:
 1.0       0.0       0.0       0.0      0.0
 0.900447  0.434966  0.0       0.0      0.0
 0.716207  1.17499   0.567585  0.0      0.0
 0.474547  1.99259   2.5554    1.2344   0.0
 0.210627  2.61987   6.63455   7.78058  3.75845

julia> legendre(Val(:schmidt), 0.45, 4, 3; ph_term = true)
5×4 Matrix{Float64}:
 1.0        0.0       0.0        0.0
 0.900447  -0.434966  0.0        0.0
 0.716207  -0.678381  0.163848   0.0
 0.474547  -0.813473  0.329901  -0.0650586
 0.210627  -0.828476  0.49451   -0.154993

julia> legendre(Val(:full), 0.45, 4, 1)
5×2 Matrix{Float64}:
 1.0       0.0
 1.55962   0.753382
 1.60149   1.51691
 1.25553   2.15225
 0.631881  2.48543
```

    legendre!(N, P::AbstractMatrix, ϕ::Number, n_max::Integer = -1, m_max::Integer = -1; kwargs...) -> Nothing

Compute the associated Legendre function `P_n,m[cos(ϕ)]`. The maximum degree and order that
will be computed are given by the parameters `n_max` and `m_max`. If they are negative
(default), the dimensions of matrix `P` will be used:

    maximum degree -> number of rows - 1
    maximum order  -> number of columns - 1

The result will be stored at matrix `P`. Hence, this function can be used to reduce
allocations.

The parameter `N` selects the normalization. The following values are valid:

- `Val(:full)`: Compute the fully normalized associated Legendre function.
- `Val(:schmidt)`: Compute the Schmidt quasi-normalized associated Legendre function.
- `Val(:unnormalized)`: Compute the unnormalized associated Legendre function.

This function has the following keywords:

- `ph_term::Bool`: If `true`, then the Condon-Shortley phase term `(-1)^m` will be included.
    (**Default** = `false`)
    
### Time-Derivative of the Legendre Associated Functions
    
This package exports two methods to compute the time-derivative of the Legendre associated
functions: `dlegendre` and `dlegendre!`.
    
    dlegendre(N, ϕ::T, n_max::Integer, m_max::Integer = -1; kwargs...) where T<:Number -> Matrix{float(T)}, Matrix{float(T)}

Compute the first-order derivative of the associated Legendre function `P_n,m[cos(ϕ)]`
with respect to `ϕ` [rad]:

    ∂P_n,m[cos(ϕ)]
    ──────────────
         ∂ϕ

The maximum degree that will be computed is `n_max` and the maximum order is `m_max`. Notice
that if `m_max` is higher than `n_max` or negative (default), it is set to `n_max`.

The parameter `N` selects the normalization. The following values are valid:

- `Val(:full)`: Compute the fully normalized associated Legendre function.
- `Val(:schmidt)`: Compute the Schmidt quasi-normalized associated Legendre function.
- `Val(:unnormalized)`: Compute the unnormalized associated Legendre function.

This function has the following keywords:

- `ph_term::Bool`: If `true`, then the Condon-Shortley phase term `(-1)^m` will be included.
    (**Default** = `false`)

It returns the following objects:

- `Matrix{float(T)}`: A matrix with the first-order derivative of the Legendre associated
    functions `P_n,m[cos(ϕ)]`.
- `Matrix{float(T)}`: A matrix with the Legendre associated functions `P_n,m[cos(ϕ)]`.

```julia
julia> dP, P = dlegendre(Val(:unnormalized), 0.45, 4)

julia> dP
5×5 Matrix{Float64}:
  0.0        0.0        0.0       0.0      0.0
 -0.434966   0.900447   0.0       0.0      0.0
 -1.17499    1.86483    2.34998   0.0      0.0
 -1.99259    1.56958    9.34577   7.6662   0.0
 -2.61987   -1.21101   19.6885   44.5626  31.1223

julia> P
5×5 Matrix{Float64}:
 1.0       0.0       0.0       0.0      0.0
 0.900447  0.434966  0.0       0.0      0.0
 0.716207  1.17499   0.567585  0.0      0.0
 0.474547  1.99259   2.5554    1.2344   0.0
 0.210627  2.61987   6.63455   7.78058  3.75845

julia> dP, P = dlegendre(Val(:schmidt), 0.45, 4, 3; ph_term = true)

julia> dP
5×4 Matrix{Float64}:
  0.0        0.0       0.0        0.0
 -0.434966  -0.900447  0.0        0.0
 -1.17499   -1.07666   0.678381   0.0
 -1.99259   -0.640778  1.20653   -0.404044
 -2.61987    0.382954  1.4675    -0.887709

julia> P
5×5 Matrix{Float64}:
 1.0        0.0       0.0        0.0        0.0
 0.900447  -0.434966  0.0        0.0        0.0
 0.716207  -0.678381  0.163848   0.0        0.0
 0.474547  -0.813473  0.329901  -0.0650586  0.0
 0.210627  -0.828476  0.49451   -0.154993   0.0264706

julia> dP, P = dlegendre(Val(:full), 0.45, 4, 1)

julia> dP
5×2 Matrix{Float64}:
  0.0        0.0
 -0.753382   1.55962
 -2.62736    2.40749
 -5.27191    1.69534
 -7.85961   -1.14886

julia> P
5×3 Matrix{Float64}:
 1.0       0.0       0.0
 1.55962   0.753382  0.0
 1.60149   1.51691   0.366375
 1.25553   2.15225   0.872836
 0.631881  2.48543   1.48353
```

    dlegendre!(dP::AbstractMatrix, ϕ::Number, P::AbstractMatrix, n_max::Integer = -1, m_max::Integer = -1; kwargs...) -> Nothing

Compute the first-order derivative of the associated Legendre function `P_n,m[x]` with
respect to `ϕ` [rad]:

    ∂P_n,m[cos(ϕ)]
    ──────────────
         ∂ϕ

The maximum degree and order that will be computed are given by the parameters `n_max` and
`m_max`. If they are negative (default), the dimensions of matrix `dP` will be used:

    maximum degree -> number of rows - 1
    maximum order  -> number of columns - 1

The derivatives will be stored in the matrix `dP`. Hence, this function can be used to reduce
allocations.

The parameter `N` selects the normalization. The following values are valid:

- `Val(:full)`: Compute the fully normalized associated Legendre function.
- `Val(:schmidt)`: Compute the Schmidt quasi-normalized associated Legendre function.
- `Val(:unnormalized)`: Compute the unnormalized associated Legendre function.

This algorithm needs the matrix `P` with the values of the associated Legendre function
using the same normalization `N`, which can be computed using the function `legendre`.

> **Warning**
> The user is responsible to pass a matrix `P` with the correct values. For example, if
> `ph_term` is `true`, `P` must also be computed with `ph_term` set to `true`.

This function has the following keywords:

- `ph_term::Bool`: If `true`, then the Condon-Shortley phase term `(-1)^m` will be included.
    (**Default** = `false`)
    
## Normalizations

### Full normalization

This algorithm was based on **[1]**. Our definition of fully normalized associated Legendre
function can be seen in **[2, p. 546]**. The conversion is obtained by:

             ┌                       ┐
             │  (n-m)! . k . (2n+1)  │
    K_n,m = √│ ───────────────────── │,  k = (m = 0) ? 1 : 2.
             │         (n+m)!        │
             └                       ┘

    P̄_n,m = P_n,m * K_n,m,

where `P̄_n,m` is the fully normalized Legendre associated function.

### Schmidt quasi-normalization

This algorithm was based on **[3, 4]**. The conversion is obtained by:

             ┌              ┐
             │      (n-m)!  │
    K_n,m = √│ k . ──────── │,  k = (m = 0) ? 1 : 2.
             │      (n+m)!  │
             └              ┘

    P̂_n,m = P_n,m * K_n,m,

where `P̂_n,m` is the Schmidt quasi-normalized Legendre associated function.

## References

- **[1]** Holmes, S. A. and W. E. Featherstone, 2002. A unified approach to the Clenshaw
    summation and the recursive computation of very high degree and order normalised
    associated Legendre functions. Journal of Geodesy, 76(5), pp. 279-299. For more info.:
    http://mitgcm.org/~mlosch/geoidcookbook/node11.html

- **[2]** Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. Microcosm
    Press, Hawthorn, CA, USA.

- **[3]** Schmidt, A (1917). Erdmagnetismus, Enzykl. Math. Wiss., 6, pp. 265–396.

- **[4]** Winch, D. E., Ivers, D. J., Turner, J. P. R., Stening R. J (2005). Geomagnetism
    and Schmidt quasi-normalization. Geophysical Journal International, 160(2), pp. 487-504.
