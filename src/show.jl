## Description #############################################################################
#
# Methods to show the types defined in this package.
#
############################################################################################

# Escape sequences of the ANSI codes used to print the types. Since this package has no
# dependencies, we use the codes directly instead of the package Crayons.jl, which is used
# by the other packages in the SatelliteToolbox.jl ecosystem.
const _B = "\e[1m" # ........................................................ Bold sequence
const _D = "\e[0m" # ....................................................... Reset sequence

function Base.show(io::IO, coefs::LegendreCoefficients{N, T}) where {N, T}
    print(
        io,
        "LegendreCoefficients{:",
        N,
        ", ",
        T,
        "}: n_max = ",
        coefs.n_max,
        ", m_max = ",
        coefs.m_max
    )

    return nothing
end

function Base.show(
    io::IO,
    ::MIME"text/plain",
    coefs::LegendreCoefficients{N, T}
) where {N, T}
    # Check for color support in the `io`.
    color = get(io, :color, false)
    b = color ? _B : ""
    d = color ? _D : ""

    println(io, "LegendreCoefficients{:", N, ", ", T, "}:")
    println(io, "$(b)   Normalization :$(d) ", _normalization_name(Val(N)))
    println(io, "$(b)  Maximum degree :$(d) ", coefs.n_max)
    print(io,   "$(b)   Maximum order :$(d) ", coefs.m_max)

    return nothing
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _normalization_name(N::Val) -> String

Return the human-readable name of the normalization related to the value `N`.
"""
_normalization_name(::Val{:full}) = "Fully normalized"

_normalization_name(::Val{:schmidt}) = "Schmidt quasi-normalized"

_normalization_name(::Val{:unnormalized}) = "Unnormalized"
