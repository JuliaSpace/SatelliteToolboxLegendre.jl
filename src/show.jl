## Description #############################################################################
#
# Methods to show the types defined in this package.
#
############################################################################################

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
    println(io, "LegendreCoefficients{:", N, ", ", T, "}:")
    println(io, styled"{bold:   Normalization :} $(_normalization_name(Val(N)))")
    println(io, styled"{bold:  Maximum degree :} $(coefs.n_max)")
    print(io,   styled"{bold:   Maximum order :} $(coefs.m_max)")

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
