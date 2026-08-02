SatelliteToolboxLegendre.jl Changelog
=====================================

Version 1.2.0
-------------

- ![Feature][badge-feature] Add the structure `LegendreCoefficients`, which stores
  precomputed recursion coefficients. The new methods `legendre!(P, ϕ, coefs, ...)` and
  `dlegendre!(dP, ϕ, P, coefs, ...)` use it to avoid evaluating square roots at every
  call, largely improving the performance when the functions are evaluated at many angles
  with the same maximum degree and order. The coefficients are stored in packed vectors
  containing only the lower triangular part of the coefficient set, halving the memory
  usage compared to dense matrices.
- ![Enhancement][badge-enhancement] Remove type instabilities when the configuration
  integers are not `Int` and when the element types of `dP` and `P` differ in the
  derivative computation.
- ![Enhancement][badge-enhancement] Improve the accuracy of the recursion coefficients of
  the fully normalized associated Legendre function at very high degrees.
- ![Enhancement][badge-enhancement] The sign adjustment applied to the derivative when
  `ϕ ∈ (π, 2π)` was numerically verified against finite differences for the entire circle,
  and the new tests cover this behavior.
- ![Bugfix][badge-bugfix] The derivative kernels wrote the result using the first indices
  of the matrix `P` instead of those of the matrix `dP`. Hence, the result was stored in
  the wrong cells if matrices with offset axes and different origins were used.
- ![Bugfix][badge-bugfix] `dlegendre!` now throws an `ArgumentError` if the matrix `P`
  does not have the dimensions required by the selected maximum degree and order. Before,
  the algorithm silently read memory outside the matrix bounds in this case, leading to
  wrong results.
- ![Info][badge-info] Many documentation errors were fixed in the docstrings and in the
  README.md file. Additionally, all the internal functions are now documented.
- ![Info][badge-info] The source code was restructured. The recursions are now implemented
  in a single generic kernel per operation, shared by all the normalizations and by the
  methods that use the precomputed coefficients. Hence, the internal functions
  `fully_normalized_legendre!`, `schmidt_quasi_normalized_legendre!`,
  `unnormalized_legendre!`, their non-mutating and derivative counterparts were removed.
  Those functions were never exported or declared public.

Version 1.1.3
-------------

- ![Enhancement][badge-enhancement] Fix allocation in Julia v1.12. ([#2][gh-pr-2])

Version 1.1.2
-------------

- ![Enhancement][badge-enhancement] Huge performance increase for lower triangular storage.

Version 1.1.1
-------------

- ![Enhancement][badge-enhancement] Performance improvements.

Version 1.1.0
-------------

- ![Info][badge-info] We dropped support for Julia 1.6. This version only supports the
  current Julia version and v1.10 (LTS).

Version 1.0.3
-------------

- ![Enhancement][badge-enhancement] Minor source-code updates.

Version 1.0.2
-------------

- ![Enhancement][badge-enhancement] Documentation update.

Version 1.0.1
-------------

- ![Bugfix][badge-bugfix] Fix keyword argument when computing Legendre associated functions
  with full normalization using the API.

Version 1.0.0
-------------

- ![Enhancement][badge-enhancement] The documentation received some minor improvements.
- ![Info][badge-info] This algorithm has been stable in **SatelliteToolbox.jl** for over six
  years. Hence, we are setting the current version as 1.0.0. The API modifications from the
  current version in **SatelliteToolbox.jl** are breaking but subtle. Therefore, we do not
  expect any problems.

Version 0.1.0
-------------

- Initial version.
  - This version was based on the code in **SatelliteToolbox.jl**.

[badge-breaking]: https://img.shields.io/badge/Breaking-DC2626?style=flat-square
[badge-deprecation]: https://img.shields.io/badge/Deprecation-D97706?style=flat-square
[badge-feature]: https://img.shields.io/badge/Feature-16A34A?style=flat-square
[badge-enhancement]: https://img.shields.io/badge/Enhancement-0284C7?style=flat-square
[badge-bugfix]: https://img.shields.io/badge/Bugfix-DB2777?style=flat-square
[badge-info]: https://img.shields.io/badge/Info-475569?style=flat-square

[gh-pr-2]: https://github.com/JuliaSpace/SatelliteToolboxLegendre.jl/pull/2
