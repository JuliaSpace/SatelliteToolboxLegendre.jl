## Description #############################################################################
#
# Tests related to the computation of the Legendre associated functions and their
#   derivatives using precomputed coefficients.
#
############################################################################################

# == Structure: LegendreCoefficients =======================================================

############################################################################################
#                                       Test Results                                       #
############################################################################################
#
# The results obtained using the precomputed coefficients must match those obtained using
# the functions that compute the coefficients at every call.
#
############################################################################################

@testset "Legendre Associated Functions" begin
    for T in (Float64, Float32)
        ϕ = T(0.123)

        for norm in (Val(:unnormalized), Val(:schmidt), Val(:full))
            for (n_max, m_max) in ((5, 5), (5, 3))
                coefs = LegendreCoefficients(norm, n_max, m_max; T = T)

                for ph_term in (false, true)
                    P_ref = legendre(norm, ϕ, n_max, m_max; ph_term = ph_term)

                    P = zeros(T, n_max + 1, m_max + 1)
                    legendre!(P, ϕ, coefs; ph_term = ph_term)

                    @test P == P_ref
                    @test eltype(P) == T
                end
            end
        end
    end
end

@testset "Derivative of the Legendre Associated Functions" begin
    for T in (Float64, Float32)
        ϕ = T(0.123)

        for norm in (Val(:unnormalized), Val(:schmidt), Val(:full))
            for (n_max, m_max) in ((5, 5), (5, 3))
                coefs = LegendreCoefficients(norm, n_max, m_max; T = T)

                for ph_term in (false, true)
                    dP_ref, P_ref = dlegendre(norm, ϕ, n_max, m_max; ph_term = ph_term)

                    P = zeros(T, size(P_ref)...)
                    legendre!(P, ϕ, coefs; ph_term = ph_term)
                    @test P == P_ref

                    dP = zeros(T, n_max + 1, m_max + 1)
                    dlegendre!(dP, ϕ, P, coefs; ph_term = ph_term)

                    @test dP ≈ dP_ref
                    @test eltype(dP) == T
                end
            end
        end
    end
end

@testset "Errors" begin
    for norm in (Val(:unnormalized), Val(:schmidt), Val(:full))
        @test_throws ArgumentError LegendreCoefficients(norm, -2)

        coefs = LegendreCoefficients(norm, 3)

        # The matrices require a degree or order higher than the ones supported by the
        # coefficients.
        P = zeros(6, 6)
        @test_throws ArgumentError legendre!(P, 0.123, coefs)

        dP = zeros(6, 6)
        @test_throws ArgumentError dlegendre!(dP, 0.123, P, coefs)

        # The matrix `P` does not have the required dimensions to compute the derivative.
        coefs = LegendreCoefficients(norm, 3, 1)
        P = zeros(4, 2)
        dP = zeros(4, 2)
        @test_throws ArgumentError dlegendre!(dP, 0.123, P, coefs)
    end
end

@testset "Show" begin
    coefs = LegendreCoefficients(Val(:full), 3, 2)

    expected = "LegendreCoefficients{:full, Float64}: n_max = 3, m_max = 2"
    @test sprint(show, coefs) == expected

    expected = """
        LegendreCoefficients{:full, Float64}:
           Normalization : Fully normalized
          Maximum degree : 3
           Maximum order : 2"""
    @test sprint(show, MIME("text/plain"), coefs) == expected

    coefs = LegendreCoefficients(Val(:schmidt), 3; T = Float32)

    expected = "LegendreCoefficients{:schmidt, Float32}: n_max = 3, m_max = 3"
    @test sprint(show, coefs) == expected

    expected = """
        LegendreCoefficients{:schmidt, Float32}:
           Normalization : Schmidt quasi-normalized
          Maximum degree : 3
           Maximum order : 3"""
    @test sprint(show, MIME("text/plain"), coefs) == expected

    coefs = LegendreCoefficients(Val(:unnormalized), 3)

    expected = """
        LegendreCoefficients{:unnormalized, Float64}:
           Normalization : Unnormalized
          Maximum degree : 3
           Maximum order : 3"""
    @test sprint(show, MIME("text/plain"), coefs) == expected
end
