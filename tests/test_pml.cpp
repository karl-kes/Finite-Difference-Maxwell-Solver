#include "test_framework.hpp"
#include "../src/Classes/Config/config.hpp"
#include "../src/Classes/PML/pml.hpp"
#include <cmath>

// PML Unit Tests

TEST(PML, ActiveWhenConfigured) {
    Simulation_Config cfg{};
    PML pml{ cfg };

    ASSERT_TRUE( pml.is_active() );
    ASSERT_EQ( pml.thickness(), std::size_t{cfg.pml_thickness} );
}

TEST(PML, InactiveWhenDisabled) {
    Simulation_Config cfg{};
    cfg.use_pml = false;
    cfg.compute_derived();
    PML pml{ cfg };

    ASSERT_TRUE( !pml.is_active() );
    ASSERT_EQ( pml.thickness(), std::size_t{0} );
}

TEST(PML, CoefficientsInValidRange) {
    Simulation_Config cfg{};
    PML pml{ cfg };

    // b coefficients should be in (0, 1] — exponential decay factors:
    for ( std::size_t i = 0; i < pml.thickness(); ++i ) {
        ASSERT_GT( pml.b_Ex_ptr()[i], 0.0 );
        ASSERT_LT( pml.b_Ex_ptr()[i], 1.0 + 1e-10 );

        ASSERT_GT( pml.b_Hx_ptr()[i], 0.0 );
        ASSERT_LT( pml.b_Hx_ptr()[i], 1.0 + 1e-10 );
    }

    // c coefficients should be non-positive (sigma*(b-1)/denom, where b<1):
    for ( std::size_t i = 0; i < pml.thickness(); ++i ) {
        ASSERT_LT( pml.c_Ex_ptr()[i], 1e-15 );
        ASSERT_LT( pml.c_Hx_ptr()[i], 1e-15 );
    }
}

TEST(PML, GradingStrongerAtOuterEdge) {
    Simulation_Config cfg{};
    PML pml{ cfg };

    std::size_t last = pml.thickness() - 1;

    // b = exp(-(sigma + alpha) * dt / eps). Outer edge (i=0) has maximum
    // sigma, so b is smallest there. Near the interior interface sigma -> 0
    // but alpha -> alpha_max, so b plateaus below 1 and is not globally
    // monotonic (alpha term can slightly dominate for the last few layers).
    ASSERT_LT( pml.b_Ex_ptr()[0], pml.b_Ex_ptr()[last] );
    ASSERT_LT( pml.b_Ex_ptr()[last], 1.0 );
    ASSERT_GT( pml.b_Ex_ptr()[last], 0.0 );

    // c is more negative (stronger correction) at the outer edge:
    ASSERT_LT( pml.c_Ex_ptr()[0], pml.c_Ex_ptr()[last] );
}

TEST(PML, SymmetricAcrossDirections) {
    Simulation_Config cfg{};
    PML pml{ cfg };

    // For isotropic grading (dx=dy=dz, same sigma_max), coefficients
    // should be identical across x, y, z:
    for ( std::size_t i = 0; i < pml.thickness(); ++i ) {
        ASSERT_NEAR( pml.b_Ex_ptr()[i], pml.b_Ey_ptr()[i], 1e-15 );
        ASSERT_NEAR( pml.b_Ex_ptr()[i], pml.b_Ez_ptr()[i], 1e-15 );
        ASSERT_NEAR( pml.c_Ex_ptr()[i], pml.c_Ey_ptr()[i], 1e-15 );
        ASSERT_NEAR( pml.c_Ex_ptr()[i], pml.c_Ez_ptr()[i], 1e-15 );
    }
}

TEST(PML, HalfIntegerShiftForHCoeffs) {
    Simulation_Config cfg{};
    PML pml{ cfg };

    // H-field coefficients use half-integer positions (i+0.5) and must
    // differ from E-field coefficients (integer positions):
    for ( std::size_t i = 0; i < pml.thickness(); ++i ) {
        if ( i < pml.thickness() - 1 ) {
            double diff_b = std::abs( pml.b_Ex_ptr()[i] - pml.b_Hx_ptr()[i] );
            ASSERT_GT( diff_b, 1e-15 );
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

TEST(PML, AnisotropicCoefficientsDistinct) {
    // With dx != dy != dz, sigma_max differs per axis, so b/c coefficients
    // for x, y, z should differ. This catches regressions to scalar sigma_max.
    Simulation_Config cfg{};
    cfg.dx = 1.0; cfg.dy = 2.0; cfg.dz = 3.0;
    cfg.compute_derived();
    PML pml{ cfg };

    // Per-axis sigma_max should differ:
    double const eta{ std::sqrt( cfg.mu / cfg.eps ) };
    double const coeff{ 0.8 * ( cfg.pml_order + 1 ) };
    ASSERT_NEAR( cfg.pml_sigma_max_x, coeff / ( cfg.dx * eta ), 1e-14 );
    ASSERT_NEAR( cfg.pml_sigma_max_y, coeff / ( cfg.dy * eta ), 1e-14 );
    ASSERT_NEAR( cfg.pml_sigma_max_z, coeff / ( cfg.dz * eta ), 1e-14 );

    // Coefficients should differ across axes at the outermost layer:
    ASSERT_TRUE( std::abs( pml.b_Ex_ptr()[0] - pml.b_Ey_ptr()[0] ) > 1e-10 );
    ASSERT_TRUE( std::abs( pml.b_Ex_ptr()[0] - pml.b_Ez_ptr()[0] ) > 1e-10 );
    ASSERT_TRUE( std::abs( pml.c_Ex_ptr()[0] - pml.c_Ey_ptr()[0] ) > 1e-10 );
}

TEST(PML, ThicknessClampedForSmallGrid) {
    // Very small grid: cbrt(10*10*10)/6 = 1.67, which would give thickness=1.
    // The clamp should enforce a minimum of 4.
    Simulation_Config cfg{};
    cfg.Nx = 10; cfg.Ny = 10; cfg.Nz = 10;
    cfg.compute_derived();
    ASSERT_GT( cfg.pml_thickness, std::size_t{3} );
    ASSERT_EQ( cfg.pml_thickness, std::size_t{4} );
}