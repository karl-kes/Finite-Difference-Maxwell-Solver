#include "test_framework.hpp"
#include "../src/Classes/Config/config.hpp"
#include "../src/Classes/PML/PML.hpp"
#include <cmath>

// ============================================================================
// PML Unit Tests
// ============================================================================

TEST(PML, ActiveWhenConfigured) {
    Simulation_Config cfg{};
    PML pml{ cfg };

    // Default config has use_pml=true, thickness=8:
    ASSERT_TRUE( pml.is_active() );
    ASSERT_EQ( pml.thickness(), std::size_t{8} );
}

TEST(PML, CoefficientsInValidRange) {
    Simulation_Config cfg{};
    PML pml{ cfg };

    // b coefficients should be in (0, 1] — they are exponential decay factors:
    for ( std::size_t i = 0; i < pml.thickness(); ++i ) {
        ASSERT_GT( pml.b_Ex_ptr()[i], 0.0 );
        ASSERT_LT( pml.b_Ex_ptr()[i], 1.0 + 1e-10 );

        ASSERT_GT( pml.b_Bx_ptr()[i], 0.0 );
        ASSERT_LT( pml.b_Bx_ptr()[i], 1.0 + 1e-10 );
    }

    // c coefficients should be non-positive (sigma*(b-1)/denom, where b<1):
    for ( std::size_t i = 0; i < pml.thickness(); ++i ) {
        ASSERT_LT( pml.c_Ex_ptr()[i], 1e-15 );  // <= 0
        ASSERT_LT( pml.c_Bx_ptr()[i], 1e-15 );
    }

    // kappa should be >= 1.0:
    for ( std::size_t i = 0; i < pml.thickness(); ++i ) {
        ASSERT_GT( pml.kappa_Ex_ptr()[i], 1.0 - 1e-10 );
        ASSERT_GT( pml.kappa_Bx_ptr()[i], 1.0 - 1e-10 );
    }
}

TEST(PML, GradingStrongerAtOuterEdge) {
    Simulation_Config cfg{};
    PML pml{ cfg };

    // Index 0 = outermost (depth_norm ~ 1.0), index thickness-1 = innermost (interface).
    // Kappa should be monotonically graded: highest at outer edge, approaching 1.0 at interface.
    // Note: b coefficients are NOT necessarily monotonic because the alpha counter-term
    // (alpha = alpha_max * (1 - depth)) adds damping near the interface that can
    // offset the sigma reduction. This is expected CPML behavior.

    std::size_t last = pml.thickness() - 1;

    // Kappa at outermost layer should be highest:
    ASSERT_GT( pml.kappa_Ex_ptr()[0], pml.kappa_Ex_ptr()[last] );

    // Innermost kappa should be close to 1.0 (minimal stretching at interface):
    ASSERT_NEAR( pml.kappa_Ex_ptr()[last], 1.0, 0.5 );

    // Kappa strictly monotonically decreasing toward interface:
    for ( std::size_t i = 1; i < pml.thickness(); ++i ) {
        ASSERT_LT( pml.kappa_Ex_ptr()[i], pml.kappa_Ex_ptr()[i - 1] + 1e-15 );
    }

    // c coefficients should become less negative (weaker correction) toward interface:
    ASSERT_LT( pml.c_Ex_ptr()[0], pml.c_Ex_ptr()[last] );
}

TEST(PML, SymmetricAcrossDirections) {
    Simulation_Config cfg{};
    PML pml{ cfg };

    // For isotropic grading (dx=dy=dz, same sigma_max), coefficients should
    // be identical across x, y, z directions:
    for ( std::size_t i = 0; i < pml.thickness(); ++i ) {
        ASSERT_NEAR( pml.b_Ex_ptr()[i], pml.b_Ey_ptr()[i], 1e-15 );
        ASSERT_NEAR( pml.b_Ex_ptr()[i], pml.b_Ez_ptr()[i], 1e-15 );
        ASSERT_NEAR( pml.c_Ex_ptr()[i], pml.c_Ey_ptr()[i], 1e-15 );
        ASSERT_NEAR( pml.c_Ex_ptr()[i], pml.c_Ez_ptr()[i], 1e-15 );
        ASSERT_NEAR( pml.kappa_Ex_ptr()[i], pml.kappa_Ey_ptr()[i], 1e-15 );
        ASSERT_NEAR( pml.kappa_Ex_ptr()[i], pml.kappa_Ez_ptr()[i], 1e-15 );
    }
}

TEST(PML, HalfIntegerShiftForBCoeffs) {
    Simulation_Config cfg{};
    PML pml{ cfg };

    // B-field coefficients use half-integer positions (i+0.5) so they
    // should differ from E-field coefficients (integer positions):
    for ( std::size_t i = 0; i < pml.thickness(); ++i ) {
        // They should NOT be identical (except possibly at edges):
        // At inner layer (i = thickness-1), depth_B can go negative, so
        // we only check layers where both depths are positive:
        if ( i < pml.thickness() - 1 ) {
            double diff_b = std::abs( pml.b_Ex_ptr()[i] - pml.b_Bx_ptr()[i] );
            double diff_k = std::abs( pml.kappa_Ex_ptr()[i] - pml.kappa_Bx_ptr()[i] );
            // At least one should differ (staggering effect):
            ASSERT_GT( diff_b + diff_k, 1e-15 );
        }
    }
}

TEST(PML, FaceSizesCorrect) {
    Simulation_Config cfg{};
    PML pml{ cfg };

    std::size_t Nx = cfg.Nx + 1;
    std::size_t Ny = cfg.Ny + 1;
    std::size_t Nz = cfg.Nz + 1;
    std::size_t t = cfg.pml_thickness;

    ASSERT_EQ( pml.psi_face_x(), t * Ny * Nz );
    ASSERT_EQ( pml.psi_face_y(), Nx * t * Nz );
    ASSERT_EQ( pml.psi_face_z(), Nx * Ny * t );
}