#include "test_framework.hpp"
#include "test_helpers.hpp"
#include "../src/Classes/Grid/grid.hpp"
#include "../src/Classes/Source/source.hpp"
#include <cmath>
#include <numbers>

// ============================================================================
// Grid Construction & Properties
// ============================================================================

TEST(Grid, DimensionsMatchConfig) {
    Simulation_Config cfg{};
    Grid grid{ cfg };

    // Grid adds +1 to each config dimension:
    ASSERT_EQ( grid.Nx(), cfg.Nx + 1 );
    ASSERT_EQ( grid.Ny(), cfg.Ny + 1 );
    ASSERT_EQ( grid.Nz(), cfg.Nz + 1 );
}

TEST(Grid, PaddedDimensionsAligned) {
    Simulation_Config cfg{};
    Grid grid{ cfg };

    constexpr std::size_t elems = SIMD_BYTES / sizeof(double);
    ASSERT_EQ( grid.Nx_padded() % elems, std::size_t{0} );
    ASSERT_EQ( grid.Ny_padded() % elems, std::size_t{0} );
    ASSERT_EQ( grid.Nz_padded() % elems, std::size_t{0} );
}

TEST(Grid, PhysicalConstantsConsistent) {
    Simulation_Config cfg{};
    Grid grid{ cfg };

    // c = 1/sqrt(mu*eps):
    double expected_c = 1.0 / std::sqrt( cfg.mu * cfg.eps );
    ASSERT_NEAR( grid.c(), expected_c, 1e-14 );

    // c^2:
    ASSERT_NEAR( grid.c_sq(), expected_c * expected_c, 1e-14 );

    // CFL condition: dt <= cfl_factor / (c * sqrt(1/dx^2 + 1/dy^2 + 1/dz^2)):
    double cfl_limit = cfg.cfl_factor / ( expected_c * std::sqrt(
        1.0/(cfg.dx*cfg.dx) + 1.0/(cfg.dy*cfg.dy) + 1.0/(cfg.dz*cfg.dz) ) );
    ASSERT_NEAR( grid.dt(), cfl_limit, 1e-14 );
}

TEST(Grid, FieldsInitializedToZero) {
    Simulation_Config cfg{};
    Grid grid{ cfg };

    // Spot-check several points in the interior:
    for ( std::size_t z = 40; z < 60; ++z ) {
        for ( std::size_t y = 40; y < 60; ++y ) {
            for ( std::size_t x = 40; x < 60; ++x ) {
                ASSERT_EQ( grid.field( Field::ELECTRIC, Component::X, x, y, z ), 0.0 );
                ASSERT_EQ( grid.field( Field::ELECTRIC, Component::Y, x, y, z ), 0.0 );
                ASSERT_EQ( grid.field( Field::ELECTRIC, Component::Z, x, y, z ), 0.0 );
                ASSERT_EQ( grid.field( Field::MAGNETIC, Component::X, x, y, z ), 0.0 );
                ASSERT_EQ( grid.field( Field::MAGNETIC, Component::Y, x, y, z ), 0.0 );
                ASSERT_EQ( grid.field( Field::MAGNETIC, Component::Z, x, y, z ), 0.0 );
            }
        }
    }
}

TEST(Grid, ZeroFieldZeroEnergy) {
    Simulation_Config cfg{};
    Grid grid{ cfg };

    ASSERT_EQ( grid.total_energy(), 0.0 );
    ASSERT_EQ( grid.source_power(), 0.0 );
}

TEST(Grid, IdxLinearization) {
    Simulation_Config cfg{};
    Grid grid{ cfg };

    // idx(0,0,0) should be 0:
    ASSERT_EQ( grid.idx(0, 0, 0), std::size_t{0} );

    // idx should be: x + Nx_padded * (y + Ny_padded * z)
    std::size_t x{5}, y{10}, z{15};
    std::size_t expected = x + grid.Nx_padded() * ( y + grid.Ny_padded() * z );
    ASSERT_EQ( grid.idx(x, y, z), expected );
}

TEST(Grid, FieldWriteRead) {
    Simulation_Config cfg{};
    Grid grid{ cfg };

    // Write through mutable reference, read through const accessor:
    grid.field( Field::ELECTRIC, Component::Y, 50, 50, 50 ) = 3.14;
    ASSERT_NEAR( grid.field( Field::ELECTRIC, Component::Y, 50, 50, 50 ), 3.14, 1e-15 );

    // Also through raw pointer:
    grid.Bz_ptr()[ grid.idx(25, 30, 35) ] = -2.71;
    ASSERT_NEAR( grid.field( Field::MAGNETIC, Component::Z, 25, 30, 35 ), -2.71, 1e-15 );
}

TEST(Grid, FieldMagnitudeComputation) {
    Simulation_Config cfg{};
    Grid grid{ cfg };

    grid.field( Field::ELECTRIC, Component::X, 50, 50, 50 ) = 3.0;
    grid.field( Field::ELECTRIC, Component::Y, 50, 50, 50 ) = 4.0;
    grid.field( Field::ELECTRIC, Component::Z, 50, 50, 50 ) = 0.0;

    ASSERT_NEAR( grid.field_magnitude( Field::ELECTRIC, 50, 50, 50 ), 5.0, 1e-14 );
}

TEST(Grid, EnergyPositiveDefinite) {
    Simulation_Config cfg{};
    Grid grid{ cfg };

    // Set a single nonzero E-field component:
    grid.field( Field::ELECTRIC, Component::X, 50, 50, 50 ) = 1.0;
    double energy = grid.total_energy();
    ASSERT_GT( energy, 0.0 );

    // Energy should be 0.5 * eps * Ex^2 * dV for that single cell:
    double dV = grid.dx() * grid.dy() * grid.dz();
    double expected = 0.5 * grid.eps() * 1.0 * dV;
    ASSERT_NEAR( energy, expected, 1e-10 );
}

// ============================================================================
// Source Application
// ============================================================================

TEST(Source, PointSourceInjectsAllComponents) {
    Simulation_Config cfg{};
    Grid grid{ cfg };

    std::size_t x{50}, y{50}, z{50};
    double val{7.5};

    grid.add_source( std::make_unique<Point_Source>( val, x, y, z ) );
    grid.apply_sources( 0 );

    std::size_t i = grid.idx(x, y, z);
    ASSERT_NEAR( grid.Jx_ptr()[i], val, 1e-15 );
    ASSERT_NEAR( grid.Jy_ptr()[i], val, 1e-15 );
    ASSERT_NEAR( grid.Jz_ptr()[i], val, 1e-15 );
}

TEST(Source, GaussianPulseTemporalShape) {
    Simulation_Config cfg{};
    Grid grid{ cfg };

    double amp{10.0};
    double t0{grid.dt() * 50.0};   // peak at step 50
    double width{grid.dt() * 10.0};
    std::size_t x{50}, y{50}, z{50};

    grid.add_source( std::make_unique<Gaussian_Pulse>( amp, t0, width, x, y, z ) );

    std::size_t i = grid.idx(x, y, z);

    // At step 50 (t=t0), Jz should be ~amplitude:
    grid.apply_sources( 50 );
    double jz_peak = grid.Jz_ptr()[i];
    ASSERT_NEAR( jz_peak, amp, 1e-10 );

    // Far from peak (step 0), should be much smaller:
    grid.apply_sources( 0 );
    double jz_far = grid.Jz_ptr()[i];
    ASSERT_LT( std::abs(jz_far), amp * 0.01 );
}

TEST(Source, StraightWireSinusoidal) {
    Simulation_Config cfg{};
    Grid grid{ cfg };

    double amp{5.0}, freq{0.1};
    std::size_t y{50}, z{50}, xs{40}, xe{60};

    grid.add_source( std::make_unique<Straight_Wire_X>( amp, freq, y, z, xs, xe ) );

    // At time_step=0, sin(0)=0, so Jx should be 0:
    grid.apply_sources( 0 );
    for ( std::size_t x = xs; x <= xe; ++x ) {
        ASSERT_NEAR( grid.Jx_ptr()[ grid.idx(x, y, z) ], 0.0, 1e-14 );
    }

    // At a later step, should be nonzero along the wire:
    grid.apply_sources( 100 );
    double sum{0};
    for ( std::size_t x = xs; x <= xe; ++x ) {
        sum += std::abs( grid.Jx_ptr()[ grid.idx(x, y, z) ] );
    }
    ASSERT_GT( sum, 0.0 );

    // All wire points should have the same current value:
    double first_val = grid.Jx_ptr()[ grid.idx(xs, y, z) ];
    for ( std::size_t x = xs + 1; x <= xe; ++x ) {
        ASSERT_NEAR( grid.Jx_ptr()[ grid.idx(x, y, z) ], first_val, 1e-15 );
    }
}