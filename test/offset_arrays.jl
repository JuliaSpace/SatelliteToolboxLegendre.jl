## Description #############################################################################
#
# Tests related to the computation of the Legendre associated functions and their
#   derivatives using matrices with offset axes.
#
############################################################################################

# == Functions: legendre! and dlegendre! ===================================================

############################################################################################
#                                       Test Results                                       #
############################################################################################
#
# The results obtained using matrices with offset axes must match those obtained using
# ordinary matrices.
#
############################################################################################

@testset "Offset Arrays" begin
    ϕ = 0.123

    for norm in (Val(:unnormalized), Val(:schmidt), Val(:full))
        for ph_term in (false, true)
            # == Reference Using Ordinary Matrices =========================================

            P_ref  = zeros(4, 4)
            dP_ref = zeros(4, 4)

            legendre!(norm, P_ref, ϕ; ph_term = ph_term)
            dlegendre!(norm, dP_ref, ϕ, P_ref; ph_term = ph_term)

            # == Computation Using Matrices With Offset Axes ===============================

            # Notice that `P` and `dP` have different origins to make sure the functions
            # treat the indices of each matrix independently.
            P  = OffsetArray(zeros(4, 4), -2:1, 0:3)
            dP = OffsetArray(zeros(4, 4), 0:3, -1:2)

            legendre!(norm, P, ϕ; ph_term = ph_term)
            dlegendre!(norm, dP, ϕ, P; ph_term = ph_term)

            @test parent(P) ≈ P_ref
            @test parent(dP) ≈ dP_ref
        end
    end
end
