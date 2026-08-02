# Repository Guide

## Package Structure

- Single Julia package with zero runtime dependencies. Requires Julia 1.10 or newer (`[compat]` in `Project.toml`).
- `src/SatelliteToolboxLegendre.jl` is the module entrypoint and controls the include order: `types.jl`, `misc.jl`, `kernels.jl`, `coefficients.jl`, `dlegendre.jl`, `legendre.jl`. All feature files are plain includes, so a file can only reference symbols defined in files included before it.
- Public API (the only exported names): `legendre`, `legendre!`, `dlegendre`, `dlegendre!`, and `LegendreCoefficients`. Everything prefixed with `_` is internal.
- The recursions are implemented once in `src/kernels.jl` as generic kernels (`_legendre_kernel!`, `_dlegendre_kernel!`) parameterized by a coefficient provider: a `Val` normalization (`Val(:full)`, `Val(:schmidt)`, `Val(:unnormalized)`) computes coefficients on the fly, while a `LegendreCoefficients` object reads precomputed arrays. Never duplicate a recursion loop; add or change coefficient providers instead.
- Tests are included unconditionally from `test/runtests.jl` into four verbose testsets: `test/legendre.jl`, `test/dlegendre.jl`, `test/offset_arrays.jl`, and `test/coefficients.jl`. There is no per-`src`-file mirroring.
- Test-only dependencies (`OffsetArrays`, `Test`) are declared in `[extras]` + `[targets]` in `Project.toml`; there is no `test/Project.toml`.
- `Manifest.toml` is not committed; do not add it to git.

## Commands

- Instantiate: `julia --project=. -e 'using Pkg; Pkg.instantiate()'`
- Full test suite: `julia --project=. -e 'using Pkg; Pkg.test()'` — the first run precompiles and can take a while; use generous timeouts.
- Focused test file: `julia --project=. -e 'using SatelliteToolboxLegendre, Test; include("test/legendre.jl")'` (same for `test/dlegendre.jl` and `test/coefficients.jl`). `test/offset_arrays.jl` additionally needs the test-only dependency `OffsetArrays`, which a plain `--project=.` session may not have; run the full suite for it. There is no test-name selector.
- CI builds the package (`julia-actions/julia-buildpkg`) before testing and covers Julia 1.10 and the latest stable release on Linux x64, macOS arm64, and Windows x64, with coverage uploaded to Codecov. A separate workflow runs the same matrix on Julia nightly.

## Code Style

- `.JuliaFormatter.toml` at the repo root is the source of truth: Blue style with alignment options enabled (`align_assignment`, `align_pair_arrow`, etc.) and several transformations disabled (`short_to_long_function_def`, `conditional_to_if`, ...). CI does not run a format check.
- Wrap code and comments at 92 characters. If a function declaration does not fit in 92 characters, break it with one argument per line; otherwise keep it on a single line.
- End every comment with a period. Use `;` to separate positional from keyword arguments in calls and definitions, and spaces around `=` in keyword arguments.
- Every function and structure has a docstring, except methods extending an already documented function. Docstring signature lines follow the `name(args; kwargs...) -> ReturnType` convention used throughout `src/`.
- Commit messages start with a gitmoji short code (e.g. `:bug:`, `:sparkles:`), use a capitalized imperative summary of at most 50 characters, and wrap the body at 72 characters.

## Behavioral Constraints

- Numerical invariant: `legendre!`/`dlegendre!` with a `LegendreCoefficients` object must produce exactly the same values as the `Val`-based methods. The tests in `test/coefficients.jl` compare with `==`, so the fill functions in `src/coefficients.jl` must keep using the same `_kernel_*` provider functions as the on-the-fly path.
- The in-place methods must remain allocation-free after warm-up; keep the `_kernel_*` provider functions `@inline` so the compiler folds them into the kernel loops.
- The kernels use the convention that `sin(ϕ)` is always positive, and the derivative applies a verified sign flip for `ϕ ∈ (π, 2π)` (see the note in `src/kernels.jl`). The finite-difference tests in `test/dlegendre.jl` lock in this behavior, including one-sided derivatives at `ϕ ∈ {0, π, 2π}`.
- The kernels do not validate inputs (they run under `@inbounds`); all dimension checks live in `_get_degree_and_order` (`src/misc.jl`) and in the public methods. Any new entry point must validate before calling a kernel.
- The kernels only write the lower triangular part of the destination matrices and support matrices with offset axes (`test/offset_arrays.jl`); always index through the `first.(axes(...))` offsets.
- New tests follow the `@testset "Name" begin ... end` pattern with hard-coded reference values, matching the existing files.

## Not Configured

- No documentation build exists: `docs/` contains only the logo asset, and there is no `docs/make.jl` (the README docs badges point to a site that is not built from this repo). Do not create one unless asked.
- No linter, Aqua/JET checks, pre-commit hooks, or CI format check are configured; do not invent them.
