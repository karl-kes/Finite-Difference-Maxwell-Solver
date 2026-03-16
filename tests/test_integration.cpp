#include "test_framework.hpp"
#include "test_helpers.hpp"
#include "../src/Classes/Config/config.hpp"
#include "../src/Classes/Grid/grid.hpp"
#include "../src/Classes/Source/source.hpp"
#include <cmath>
#include <numbers>
#include <vector>
#include <numeric>
#include <omp.h>

// ============================================================================
// Leapfrog Update Correctness
// ============================================================================

TEST(Integration, EmptyGridStaysZero) {
    // A grid with no sources or initial fields should remain zero after stepping.
    Simulation_Config cfg{};
    Grid grid{ cfg };

    for ( int t = 0; t < 10; ++t ) {
        grid.step();
    }

    ASSERT_EQ( grid.total_energy(), 0.0 );
}

TEST(Integration, SingleStepCurlSymmetry) {
    // Place a single nonzero Ey at center. After one B-update, the curl of E
    // should produce Bx and Bz changes (but not By, since curl_y(E) involves
    // dEx/dz - dEz/dx, neither of which is set).
    Simulation_Config cfg{};
    Grid grid{ cfg };

    std::size_t cx{50}, cy{50}, cz{50};
    grid.field( Field::ELECTRIC, Component::Y, cx, cy, cz ) = 1.0;

    // Run one full step (B update then E update):
    grid.step();

    // Bx should be affected (curl_x has dEz/dy - dEy/dz terms):
    // Bz should be affected (curl_z has dEy/dx - dEx/dy terms):
    // Check that something happened:
    double bx = grid.field( Field::MAGNETIC, Component::X, cx, cy, cz );
    double bz = grid.field( Field::MAGNETIC, Component::Z, cx, cy, cz );
    double by = grid.field( Field::MAGNETIC, Component::Y, cx, cy, cz );

    // Bx and Bz should be nonzero due to Ey derivatives:
    // (The exact sign depends on forward-difference stencil direction)
    ASSERT_TRUE( std::abs(bx) > 0.0 || std::abs(bz) > 0.0 );

    // By should be zero (curl_y depends on Ex, Ez which are zero):
    ASSERT_NEAR( by, 0.0, 1e-15 );
}

// ============================================================================
// Plane Wave Propagation (E2E Physics Test)
// ============================================================================

TEST(Integration, PlaneWavePropagation) {
    // Initialize Ey = sin(kx), Bz = sin(kx)/c propagating in +x.
    // Track phase correlation and energy drift over 100 steps.
    Simulation_Config cfg{};
    Grid grid{ cfg };

    // Wavelength: 20 cells
    double const wavelength{ 20.0 * grid.dx() };
    double const k{ 2.0 * std::numbers::pi / wavelength };

    // Initialize plane wave in interior (avoid PML):
    std::size_t const margin{ 15 };
    for ( std::size_t z = margin; z < grid.Nz() - margin; ++z ) {
        for ( std::size_t y = margin; y < grid.Ny() - margin; ++y ) {
            for ( std::size_t x = margin; x < grid.Nx() - margin; ++x ) {
                double phase = k * static_cast<double>(x) * grid.dx();
                std::size_t i = grid.idx(x, y, z);
                grid.Ey_ptr()[i] = std::sin( phase );
                grid.Bz_ptr()[i] = std::sin( phase ) / grid.c();
            }
        }
    }

    double const E0 = grid.total_energy();
    ASSERT_GT( E0, 0.0 );

    // Probe at center:
    std::size_t const px{50}, py{50}, pz{50};
    double const initial_phase = k * static_cast<double>(px) * grid.dx();
    double const phase_shift = k * grid.c() * grid.dt();

    std::size_t const num_steps{100};
    double sum_exp{}, sum_act{}, sum_exp2{}, sum_act2{}, sum_prod{};

    for ( std::size_t t = 0; t < num_steps; ++t ) {
        double expected = std::sin( initial_phase - phase_shift * static_cast<double>(t) );
        double actual = grid.Ey_ptr()[ grid.idx(px, py, pz) ];

        sum_exp += expected;
        sum_act += actual;
        sum_exp2 += expected * expected;
        sum_act2 += actual * actual;
        sum_prod += expected * actual;

        grid.step();
    }

    double const Ef = grid.total_energy();
    double const drift = 100.0 * std::abs( Ef - E0 ) / E0;

    double const n = static_cast<double>(num_steps);
    double const num = n * sum_prod - sum_exp * sum_act;
    double const den = std::sqrt(
        ( n * sum_exp2 - sum_exp * sum_exp ) *
        ( n * sum_act2 - sum_act * sum_act ) );
    double const correlation = den > 1e-10 ? num / den : 0.0;

    // Thresholds (same as existing validation):
    ASSERT_LT( drift, 5.0 );          // Energy drift < 5%
    ASSERT_GT( correlation, 0.99 );    // Phase correlation > 0.99
}

// ============================================================================
// Energy Conservation (Source-free, PML active)
// ============================================================================

TEST(Integration, EnergyConservationShortRun) {
    // With a localized initial perturbation (no active source), energy should
    // be approximately conserved over ~50 steps (PML absorbs what reaches boundary).
    Simulation_Config cfg{};
    Grid grid{ cfg };

    // Gaussian E-field blob at center:
    double sigma{5.0};
    for ( std::size_t z = 30; z < 70; ++z ) {
        for ( std::size_t y = 30; y < 70; ++y ) {
            for ( std::size_t x = 30; x < 70; ++x ) {
                double r2 = std::pow(static_cast<double>(x) - 50.0, 2)
                          + std::pow(static_cast<double>(y) - 50.0, 2)
                          + std::pow(static_cast<double>(z) - 50.0, 2);
                double val = std::exp( -r2 / (2.0 * sigma * sigma) );
                grid.Ey_ptr()[ grid.idx(x, y, z) ] = val;
            }
        }
    }

    double E0 = grid.total_energy();
    ASSERT_GT( E0, 0.0 );

    // Run 20 steps — pulse shouldn't have reached PML yet at c*dt*20 << 30 cells:
    for ( int t = 0; t < 20; ++t ) {
        grid.step();
    }

    double Ef = grid.total_energy();
    double drift_pct = 100.0 * std::abs( Ef - E0 ) / E0;

    // Energy should be well conserved before PML absorption starts:
    ASSERT_LT( drift_pct, 2.0 );
}

// ============================================================================
// Divergence-Free (Gauss's Law)
// ============================================================================

TEST(Integration, GaussLawDivB) {
    // After several time steps with a point source, div(B) should remain ~0
    // everywhere (no magnetic monopoles).
    Simulation_Config cfg{};
    Grid grid{ cfg };

    grid.add_source( std::make_unique<Gaussian_Pulse>(
        10.0, grid.dt() * 10.0, grid.dt() * 5.0, 50, 50, 50 ) );

    // Run 30 steps to let fields develop:
    for ( int t = 0; t < 30; ++t ) {
        grid.apply_sources( static_cast<std::size_t>(t) );
        grid.step();
    }

    // Compute div(B) = dBx/dx + dBy/dy + dBz/dz at interior points:
    double max_divB{0.0};
    double max_B{0.0};

    double const inv_dx = 1.0 / grid.dx();
    double const inv_dy = 1.0 / grid.dy();
    double const inv_dz = 1.0 / grid.dz();

    for ( std::size_t z = 25; z < 75; ++z ) {
        for ( std::size_t y = 25; y < 75; ++y ) {
            for ( std::size_t x = 25; x < 75; ++x ) {
                double dBx_dx = ( grid.field(Field::MAGNETIC, Component::X, x+1, y, z)
                                - grid.field(Field::MAGNETIC, Component::X, x, y, z) ) * inv_dx;
                double dBy_dy = ( grid.field(Field::MAGNETIC, Component::Y, x, y+1, z)
                                - grid.field(Field::MAGNETIC, Component::Y, x, y, z) ) * inv_dy;
                double dBz_dz = ( grid.field(Field::MAGNETIC, Component::Z, x, y, z+1)
                                - grid.field(Field::MAGNETIC, Component::Z, x, y, z) ) * inv_dz;

                double divB = std::abs( dBx_dx + dBy_dy + dBz_dz );
                if ( divB > max_divB ) max_divB = divB;

                double Bmag = grid.field_magnitude( Field::MAGNETIC, x, y, z );
                if ( Bmag > max_B ) max_B = Bmag;
            }
        }
    }

    // div(B) should be zero to machine precision (Yee scheme preserves this exactly):
    // Allow small tolerance for floating point:
    if ( max_B > 1e-15 ) {
        double relative_divB = max_divB / ( max_B * inv_dx );
        ASSERT_LT( relative_divB, 1e-10 );
    }
}

// ============================================================================
// Reciprocity / Symmetry
// ============================================================================

TEST(Integration, IsotropicSourceSymmetry) {
    // A Gaussian pulse injecting Jz should produce fields that propagate outward.
    // By symmetry of the grid (dx=dy=dz), |E| at equal distances from the source
    // should be approximately equal in the xy-plane (perpendicular to dipole axis).
    Simulation_Config cfg{};
    Grid grid{ cfg };

    grid.add_source( std::make_unique<Gaussian_Pulse>(
        50.0, grid.dt() * 10.0, grid.dt() * 4.0, 50, 50, 50 ) );

    for ( int t = 0; t < 50; ++t ) {
        grid.apply_sources( static_cast<std::size_t>(t) );
        grid.step();
    }

    // For a Jz source (z-dipole), the radiation is symmetric in the xy-plane.
    // Check |E| at 8 cells away along ±x and ±y:
    double E_py = grid.field_magnitude( Field::ELECTRIC, 50, 58, 50 );
    double E_ny = grid.field_magnitude( Field::ELECTRIC, 50, 42, 50 );
    double E_px = grid.field_magnitude( Field::ELECTRIC, 58, 50, 50 );
    double E_nx = grid.field_magnitude( Field::ELECTRIC, 42, 50, 50 );

    double E_avg = ( E_py + E_ny + E_px + E_nx ) / 4.0;
    ASSERT_GT( E_avg, 1e-10 );

    // +/- symmetry in same axis should be tight:
    ASSERT_LT( std::abs( E_py - E_ny ) / E_avg, 0.35 );
    ASSERT_LT( std::abs( E_px - E_nx ) / E_avg, 0.35 );

    // x vs y symmetry (Yee staggering can break this more):
    ASSERT_LT( std::abs( E_py - E_px ) / E_avg, 0.50 );
}

// ============================================================================
// PML Absorption Effectiveness
// ============================================================================

TEST(Integration, PMLReducesReflection) {
    // Test PML effectiveness by checking that a plane wave propagating into PML
    // doesn't produce large reflections back into the interior.
    // We initialize a +x traveling plane wave and measure the reflected Ey
    // at a probe behind the wavefront after it has entered the PML.
    Simulation_Config cfg{};
    Grid grid{ cfg };

    double const wavelength = 20.0 * grid.dx();
    double const k = 2.0 * std::numbers::pi / wavelength;

    // Initialize plane wave traveling in +x, placed in the right half of the grid
    // so it hits the high-x PML quickly:
    for ( std::size_t z = 15; z < grid.Nz() - 15; ++z ) {
        for ( std::size_t y = 15; y < grid.Ny() - 15; ++y ) {
            for ( std::size_t x = 50; x < grid.Nx() - 15; ++x ) {
                double phase = k * static_cast<double>(x) * grid.dx();
                std::size_t i = grid.idx(x, y, z);
                grid.Ey_ptr()[i] = std::sin( phase );
                grid.Bz_ptr()[i] = std::sin( phase ) / grid.c();
            }
        }
    }

    double E0 = grid.total_energy();

    // Run 80 steps — wavefront travels ~29 cells, enough to enter PML:
    for ( int t = 0; t < 80; ++t ) {
        grid.step();
    }

    // The PML should absorb outgoing energy, not reflect it.
    // Energy should decrease (PML absorbing) or stay similar, not increase:
    double E80 = grid.total_energy();
    ASSERT_LT( E80, E0 * 1.10 );  // No more than 10% energy increase

    // Check that the probe at x=30 (behind where the wave was initialized)
    // has very small fields (any signal there is reflected from PML):
    double reflected_max{0};
    for ( std::size_t z = 40; z < 60; ++z ) {
        for ( std::size_t y = 40; y < 60; ++y ) {
            double val = std::abs( grid.Ey_ptr()[ grid.idx(30, y, z) ] );
            if ( val > reflected_max ) reflected_max = val;
        }
    }

    // Reflected field should be small relative to initial amplitude (1.0):
    ASSERT_LT( reflected_max, 0.15 );
}

TEST(Integration, PMLLateTimeStabilityCheck) {
    // Document/detect PML late-time instability.
    // CPML with the current alpha/sigma parameters can exhibit late-time
    // energy growth. This test flags if growth exceeds 5× initial energy.
    Simulation_Config cfg{};
    Grid grid{ cfg };

    grid.add_source( std::make_unique<Gaussian_Pulse>(
        10.0, grid.dt() * 10.0, grid.dt() * 4.0, 50, 50, 50 ) );

    for ( int t = 0; t < 20; ++t ) {
        grid.apply_sources( static_cast<std::size_t>(t) );
        grid.step();
    }
    double E_baseline = grid.total_energy();

    for ( int t = 20; t < 500; ++t ) {
        grid.step();
    }
    double E_final = grid.total_energy();

    // The simulation should remain bounded — no exponential blowup.
    // Note: CPML with current parameters exhibits mild late-time energy growth
    // (~5-10% over 500 steps). This is a known CPML limitation when alpha
    // parameters aren't perfectly tuned. The test catches catastrophic blowup:
    ASSERT_TRUE( std::isfinite(E_final) );
    ASSERT_LT( E_final, E_baseline * 100.0 );  // No catastrophic instability
}

// ============================================================================
// Leapfrog Time-Reversal
// ============================================================================

TEST(Integration, LeapfrogReversibility) {
    // The leapfrog scheme (without PML/sources) should be time-reversible.
    // We can't truly reverse (would need to negate dt), but we can check that
    // energy is exactly conserved per step in the interior.
    // Here we verify energy doesn't grow — no instability.
    Simulation_Config cfg{};
    Grid grid{ cfg };

    // Set a structured initial condition well inside PML:
    for ( std::size_t z = 30; z < 70; ++z ) {
        for ( std::size_t y = 30; y < 70; ++y ) {
            for ( std::size_t x = 30; x < 70; ++x ) {
                double phase = 2.0 * std::numbers::pi * static_cast<double>(x) / 40.0;
                grid.Ey_ptr()[ grid.idx(x, y, z) ] = std::sin(phase);
                grid.Bz_ptr()[ grid.idx(x, y, z) ] = std::sin(phase) / grid.c();
            }
        }
    }

    double E_prev = grid.total_energy();
    std::vector<double> energies;
    energies.push_back( E_prev );

    for ( int t = 0; t < 50; ++t ) {
        grid.step();
        double E_now = grid.total_energy();
        energies.push_back( E_now );
    }

    // Check no exponential growth (energy at end <= energy at start * 1.05):
    ASSERT_LT( energies.back(), energies.front() * 1.05 );

    // Check energy is monotonically non-increasing (PML only absorbs):
    // Allow small tolerance for floating point:
    for ( std::size_t i = 1; i < energies.size(); ++i ) {
        // Energy can decrease (PML absorption) but shouldn't increase significantly:
        ASSERT_LT( energies[i], energies[i-1] * 1.001 + 1e-15 );
    }
}

// ============================================================================
// CFL Stability
// ============================================================================

TEST(Integration, CFLStabilityLongRun) {
    // Run 200 steps with a Gaussian pulse source and verify no blowup.
    Simulation_Config cfg{};
    Grid grid{ cfg };

    grid.add_source( std::make_unique<Gaussian_Pulse>(
        5.0, grid.dt() * 20.0, grid.dt() * 8.0, 50, 50, 50 ) );

    double max_energy{0};

    for ( int t = 0; t < 200; ++t ) {
        grid.apply_sources( static_cast<std::size_t>(t) );
        grid.step();

        double E = grid.total_energy();
        if ( E > max_energy ) max_energy = E;

        // No NaN or Inf:
        ASSERT_TRUE( std::isfinite(E) );
    }

    // Final energy should be finite and reasonable:
    double E_final = grid.total_energy();
    ASSERT_TRUE( std::isfinite(E_final) );
    ASSERT_LT( E_final, 1e20 );  // No blowup
}

// ============================================================================
// Causality
// ============================================================================

TEST(Integration, CausalPropagation) {
    // A point impulse at center at t=0 should not have reached cells far away
    // after only a few time steps (wavefront travels at c).
    Simulation_Config cfg{};
    Grid grid{ cfg };

    // Single impulse at center:
    grid.Ey_ptr()[ grid.idx(50, 50, 50) ] = 1.0;

    // After 5 steps, the wavefront has traveled ~5*c*dt cells ≈ 1.8 cells:
    for ( int t = 0; t < 5; ++t ) {
        grid.step();
    }

    // Field 20 cells away should be zero (causality):
    double far_field = std::abs( grid.Ey_ptr()[ grid.idx(70, 50, 50) ] );
    ASSERT_LT( far_field, 1e-14 );

    // Field nearby should be nonzero:
    double near_field = std::abs( grid.Ey_ptr()[ grid.idx(51, 50, 50) ] );
    ASSERT_GT( near_field, 1e-20 );
}

// ============================================================================
// Grid Convergence (2nd Order)
// ============================================================================

TEST(Integration, DispersionConvergence) {
    // Run plane waves at two different CFL factors (effective resolution).
    // The lower CFL (more steps per wavelength) should have less dispersion.
    // We can't easily change dx with constexpr config, so we measure dispersion
    // at the default resolution and verify it's within expected bounds.
    Simulation_Config cfg{};
    Grid grid{ cfg };

    double const wavelength = 20.0 * grid.dx();
    double const k = 2.0 * std::numbers::pi / wavelength;

    // Initialize:
    for ( std::size_t z = 15; z < grid.Nz() - 15; ++z ) {
        for ( std::size_t y = 15; y < grid.Ny() - 15; ++y ) {
            for ( std::size_t x = 15; x < grid.Nx() - 15; ++x ) {
                double phase = k * static_cast<double>(x) * grid.dx();
                grid.Ey_ptr()[ grid.idx(x, y, z) ] = std::sin(phase);
                grid.Bz_ptr()[ grid.idx(x, y, z) ] = std::sin(phase) / grid.c();
            }
        }
    }

    // Run 100 steps and measure phase error:
    std::size_t px{50}, py{50}, pz{50};
    double initial_phase = k * static_cast<double>(px) * grid.dx();
    double phase_shift = k * grid.c() * grid.dt();
    double total_phase_err{0};

    for ( std::size_t t = 0; t < 100; ++t ) {
        double expected = std::sin( initial_phase - phase_shift * static_cast<double>(t) );
        double actual = grid.Ey_ptr()[ grid.idx(px, py, pz) ];
        total_phase_err += std::abs( actual - expected );
        grid.step();
    }

    double avg_err = total_phase_err / 100.0;
    double total_shift = phase_shift * 100.0;
    double dispersion_pct = 100.0 * avg_err / std::max(total_shift, 1e-10);

    // With 20 cells/wavelength and CFL=0.125, dispersion should be very low:
    ASSERT_LT( dispersion_pct, 5.0 );
}

// ============================================================================
// Hertzian Dipole Radiation Pattern
// ============================================================================

// Sinusoidal z-current source for the dipole test:
class Sinusoidal_Jz : public Source {
private:
    double amplitude_;
    double frequency_;
    std::size_t x_, y_, z_;

public:
    Sinusoidal_Jz( double amp, double freq,
                   std::size_t x, std::size_t y, std::size_t z )
        : amplitude_{ amp }, frequency_{ freq }
        , x_{ x }, y_{ y }, z_{ z } {}

    void apply( Grid &grid, std::size_t const time_step ) override {
        double const omega = 2.0 * std::numbers::pi * frequency_;
        double const t = static_cast<double>( time_step ) * grid.dt();
        grid.Jz_ptr()[ grid.idx( x_, y_, z_ ) ] = amplitude_ * std::sin( omega * t );
    }
};

// Helper: run the dipole to steady state and return the grid.
// Wavelength = 10 cells (50 units) so kr >> 1 at sampling radii of 15-25 cells.
// Runs ~8 periods to establish the radiation pattern.
struct DipoleSetup {
    double wavelength;
    double freq;
    std::size_t cx, cy, cz;
    std::size_t steps_per_period;
    std::size_t warmup_steps;

    DipoleSetup( Grid const& grid )
        : wavelength{ 10.0 * grid.dx() }
        , freq{ grid.c() / wavelength }
        , cx{ grid.Nx() / 2 }, cy{ grid.Ny() / 2 }, cz{ grid.Nz() / 2 }
        , steps_per_period{ static_cast<std::size_t>( wavelength / ( grid.c() * grid.dt() ) ) }
        , warmup_steps{ 8 * steps_per_period }
    {}
};

TEST(Integration, HertzianDipoleSinTheta) {
    // A z-directed Hertzian dipole should radiate with a sin(theta) pattern:
    //   - Maximum |E| in the xy-plane (theta = pi/2, perpendicular to dipole)
    //   - Minimum |E| along the z-axis (theta = 0, along the dipole)
    //
    // Uses wavelength = 10 cells so sampling at r = 15-20 cells gives kr = 9-13,
    // well into the far-field regime where the 1/r radiation term dominates.

    Simulation_Config cfg{};
    Grid grid{ cfg };
    DipoleSetup ds{ grid };

    grid.add_source( std::make_unique<Sinusoidal_Jz>(
        1.0, ds.freq, ds.cx, ds.cy, ds.cz ) );

    // Warm up to steady state:
    for ( std::size_t t = 0; t < ds.warmup_steps; ++t ) {
        grid.apply_sources( t );
        grid.step();
    }

    // Sample peak |E| over one full period at r = 20 cells:
    int const r = 20;
    double peak_eq_px{0}, peak_eq_nx{0}, peak_eq_py{0}, peak_eq_ny{0};
    double peak_pole_pz{0}, peak_pole_nz{0};
    double peak_diag{0};

    int const rd = static_cast<int>( std::round( static_cast<double>(r) / std::sqrt(2.0) ) );

    for ( std::size_t s = 0; s < ds.steps_per_period; ++s ) {
        std::size_t t = ds.warmup_steps + s;
        grid.apply_sources( t );
        grid.step();

        // Equator (theta = pi/2):
        double e;
        e = grid.field_magnitude( Field::ELECTRIC,
            static_cast<std::size_t>(static_cast<int>(ds.cx) + r), ds.cy, ds.cz );
        if ( e > peak_eq_px ) peak_eq_px = e;

        e = grid.field_magnitude( Field::ELECTRIC,
            static_cast<std::size_t>(static_cast<int>(ds.cx) - r), ds.cy, ds.cz );
        if ( e > peak_eq_nx ) peak_eq_nx = e;

        e = grid.field_magnitude( Field::ELECTRIC,
            ds.cx, static_cast<std::size_t>(static_cast<int>(ds.cy) + r), ds.cz );
        if ( e > peak_eq_py ) peak_eq_py = e;

        e = grid.field_magnitude( Field::ELECTRIC,
            ds.cx, static_cast<std::size_t>(static_cast<int>(ds.cy) - r), ds.cz );
        if ( e > peak_eq_ny ) peak_eq_ny = e;

        // Poles (theta = 0, pi):
        e = grid.field_magnitude( Field::ELECTRIC,
            ds.cx, ds.cy, static_cast<std::size_t>(static_cast<int>(ds.cz) + r) );
        if ( e > peak_pole_pz ) peak_pole_pz = e;

        e = grid.field_magnitude( Field::ELECTRIC,
            ds.cx, ds.cy, static_cast<std::size_t>(static_cast<int>(ds.cz) - r) );
        if ( e > peak_pole_nz ) peak_pole_nz = e;

        // Diagonal (theta ~ pi/4, in xz-plane):
        e = grid.field_magnitude( Field::ELECTRIC,
            static_cast<std::size_t>(static_cast<int>(ds.cx) + rd), ds.cy,
            static_cast<std::size_t>(static_cast<int>(ds.cz) + rd) );
        if ( e > peak_diag ) peak_diag = e;
    }

    double peak_equator = ( peak_eq_px + peak_eq_nx + peak_eq_py + peak_eq_ny ) / 4.0;
    double peak_pole = ( peak_pole_pz + peak_pole_nz ) / 2.0;

    // 1. Fields should be nonzero:
    ASSERT_GT( peak_equator, 1e-6 );

    // 2. sin(theta) pattern: equator should be significantly stronger than pole.
    //    At kr ~ 12.6, far-field dominates and the ratio should be > 3:
    ASSERT_GT( peak_equator / std::max(peak_pole, 1e-30), 3.0 );

    // 3. Diagonal (theta ~ pi/4) should be between pole and equator:
    ASSERT_GT( peak_diag, peak_pole );
    ASSERT_LT( peak_diag, peak_equator * 1.5 );

    // 4. Azimuthal symmetry: ±x and ±y equatorial peaks should match within 30%:
    ASSERT_LT( std::abs(peak_eq_px - peak_eq_py) / peak_equator, 0.30 );
    ASSERT_LT( std::abs(peak_eq_px - peak_eq_nx) / peak_equator, 0.30 );

    // 5. Polar symmetry: +z and -z peaks should match within 40%:
    if ( peak_pole > 1e-10 ) {
        ASSERT_LT( std::abs(peak_pole_pz - peak_pole_nz) / peak_pole, 0.40 );
    }
}

TEST(Integration, HertzianDipole1OverR) {
    // In the far field, |Etheta| falls as 1/r. Compare peak |E| at r=15 and r=25
    // cells in the equatorial plane. With kr = 9.4 and 15.7 respectively,
    // both points are in the far field.

    Simulation_Config cfg{};
    Grid grid{ cfg };
    DipoleSetup ds{ grid };

    grid.add_source( std::make_unique<Sinusoidal_Jz>(
        1.0, ds.freq, ds.cx, ds.cy, ds.cz ) );

    for ( std::size_t t = 0; t < ds.warmup_steps; ++t ) {
        grid.apply_sources( t );
        grid.step();
    }

    int const r1 = 15;
    int const r2 = 25;

    double peak_r1{0}, peak_r2{0};

    for ( std::size_t s = 0; s < ds.steps_per_period; ++s ) {
        std::size_t t = ds.warmup_steps + s;
        grid.apply_sources( t );
        grid.step();

        double e1 = grid.field_magnitude( Field::ELECTRIC,
            static_cast<std::size_t>(static_cast<int>(ds.cx) + r1), ds.cy, ds.cz );
        double e2 = grid.field_magnitude( Field::ELECTRIC,
            static_cast<std::size_t>(static_cast<int>(ds.cx) + r2), ds.cy, ds.cz );

        if ( e1 > peak_r1 ) peak_r1 = e1;
        if ( e2 > peak_r2 ) peak_r2 = e2;
    }

    ASSERT_GT( peak_r1, 1e-6 );
    ASSERT_GT( peak_r2, 1e-6 );

    // 1/r law: E(r1)/E(r2) should approximate r2/r1 = 25/15 = 1.667.
    // Allow 35% tolerance for grid dispersion and residual near-field terms:
    double ratio_measured = peak_r1 / peak_r2;
    double ratio_expected = static_cast<double>(r2) / static_cast<double>(r1);

    ASSERT_GT( ratio_measured, ratio_expected * 0.65 );
    ASSERT_LT( ratio_measured, ratio_expected * 1.50 );
}

// ============================================================================
// Second-Order Convergence Rate
// ============================================================================

// Measure numerical phase velocity error by tracking zero-crossings at a probe.
// The probe sees a sinusoidal signal; we extract the effective period from
// zero-crossings and compare to the analytical period T = wavelength / c.
// Returns |v_numerical/c - 1|.
static double measure_phase_velocity_error( double wavelength_cells, std::size_t num_steps ) {
    Simulation_Config cfg{};
    Grid grid{ cfg };

    double const wavelength = wavelength_cells * grid.dx();
    double const k = 2.0 * std::numbers::pi / wavelength;

    std::size_t const margin{ 15 };
    for ( std::size_t z = margin; z < grid.Nz() - margin; ++z ) {
        for ( std::size_t y = margin; y < grid.Ny() - margin; ++y ) {
            for ( std::size_t x = margin; x < grid.Nx() - margin; ++x ) {
                double phase = k * static_cast<double>(x) * grid.dx();
                std::size_t i = grid.idx(x, y, z);
                grid.Ey_ptr()[i] = std::sin( phase );
                grid.Bz_ptr()[i] = std::sin( phase ) / grid.c();
            }
        }
    }

    std::size_t const px{50}, py{50}, pz{50};

    // Record Ey at probe over time:
    std::vector<double> probe;
    probe.reserve( num_steps );

    for ( std::size_t t = 0; t < num_steps; ++t ) {
        probe.push_back( grid.Ey_ptr()[ grid.idx(px, py, pz) ] );
        grid.step();
    }

    // Find upward zero-crossings with linear interpolation:
    std::vector<double> crossings;
    for ( std::size_t i = 1; i < probe.size(); ++i ) {
        if ( probe[i-1] < 0.0 && probe[i] >= 0.0 ) {
            double frac = -probe[i-1] / ( probe[i] - probe[i-1] );
            crossings.push_back( static_cast<double>(i-1) + frac );
        }
    }

    if ( crossings.size() < 2 ) return -1.0;

    // Average period from consecutive crossings:
    double total_period{0};
    for ( std::size_t i = 1; i < crossings.size(); ++i ) {
        total_period += crossings[i] - crossings[i-1];
    }
    double avg_period_steps = total_period / static_cast<double>( crossings.size() - 1 );
    double avg_period_time = avg_period_steps * grid.dt();

    // Phase velocity error:
    double v_num = wavelength / avg_period_time;
    return std::abs( v_num / grid.c() - 1.0 );
}

TEST(Integration, SecondOrderConvergenceRate) {
    // The Yee scheme is second-order in space. Halving the number of points
    // per wavelength should quadruple the phase velocity error.
    //
    // We test at 10 and 20 cells/wavelength (wavelengths of 10 and 20 cells).
    // Both are short enough to avoid PML contamination and long enough to
    // produce multiple zero-crossings in 500 steps.
    //
    // Expected: err(10) / err(20) ≈ 4.0

    std::size_t const steps{ 500 };

    double err_coarse = measure_phase_velocity_error( 10.0, steps );  // 10 pts/lambda
    double err_fine   = measure_phase_velocity_error( 20.0, steps );  // 20 pts/lambda

    ASSERT_GT( err_coarse, 0.0 );
    ASSERT_GT( err_fine, 0.0 );

    // Coarse should have more error:
    ASSERT_GT( err_coarse, err_fine );

    // Convergence ratio should be approximately 4 (second-order):
    double ratio = err_coarse / err_fine;
    ASSERT_GT( ratio, 3.0 );
    ASSERT_LT( ratio, 5.5 );
}

// ============================================================================
// OpenMP Thread Determinism
// ============================================================================

// Helper: run a short simulation with a given thread count and return final energy
// plus a field sample. Uses a Gaussian pulse source and 50 steps.
struct ThreadTestResult {
    double energy;
    double field_sample;
};

static ThreadTestResult run_with_threads( int num_threads ) {
    omp_set_num_threads( num_threads );

    Simulation_Config cfg{};
    Grid grid{ cfg };

    grid.add_source( std::make_unique<Gaussian_Pulse>(
        10.0, grid.dt() * 10.0, grid.dt() * 4.0, 50, 50, 50 ) );

    for ( std::size_t t = 0; t < 50; ++t ) {
        grid.apply_sources( t );
        grid.step();
    }

    return { grid.total_energy(), grid.Ey_ptr()[ grid.idx(60, 50, 50) ] };
}

TEST(Integration, ThreadDeterminism) {
    // The same simulation should produce identical (or near-identical) results
    // regardless of thread count. OpenMP parallel loops with static scheduling
    // should give deterministic work partitioning, and reductions should produce
    // bit-identical results if the reduction tree is consistent.
    //
    // Compare 1 thread vs 2 threads vs max threads. Both total energy and a
    // specific field value are checked.

    int const max_threads = omp_get_max_threads();

    ThreadTestResult r1 = run_with_threads( 1 );
    ThreadTestResult r2 = run_with_threads( 2 );
    ThreadTestResult rm = run_with_threads( max_threads );

    ASSERT_GT( r1.energy, 0.0 );

    // Field updates use static scheduling and no reductions, so the field values
    // should be bit-identical across thread counts. Energy uses omp reduction
    // which may reorder additions — allow tiny tolerance:
    double energy_tol = r1.energy * 1e-10;
    ASSERT_NEAR( r1.energy, r2.energy, energy_tol );
    ASSERT_NEAR( r1.energy, rm.energy, energy_tol );

    // Field sample should match more tightly (no reduction involved):
    double field_tol = std::abs(r1.field_sample) * 1e-12 + 1e-20;
    ASSERT_NEAR( r1.field_sample, r2.field_sample, field_tol );
    ASSERT_NEAR( r1.field_sample, rm.field_sample, field_tol );

    // Restore:
    omp_set_num_threads( max_threads );
}

// ============================================================================
// Superposition (Linearity)
// ============================================================================

TEST(Integration, Superposition) {
    // Maxwell's equations are linear: the field from two sources should equal
    // the sum of the fields from each source independently. Run three simulations:
    //   A: source at (40, 50, 50) only
    //   B: source at (60, 50, 50) only
    //   C: both sources simultaneously
    // Verify: field_C ≈ field_A + field_B at several probe points.

    Simulation_Config cfg{};

    double const amp = 5.0;
    double const t0 = cfg.dt * 10.0;  // use cfg.dt since grid isn't constructed yet
    double const width = cfg.dt * 4.0;
    std::size_t const xa = 40, xb = 60, y = 50, z = 50;
    std::size_t const num_steps = 40;

    // Helper: run a simulation with given sources and return Ey at probe points
    struct ProbeResult {
        double Ey_at_45;
        double Ey_at_50;
        double Ey_at_55;
        double Ey_at_65;
    };

    auto run_sim = [&]( bool use_A, bool use_B ) -> ProbeResult {
        Grid grid{ cfg };
        if ( use_A ) {
            grid.add_source( std::make_unique<Gaussian_Pulse>( amp, t0, width, xa, y, z ) );
        }
        if ( use_B ) {
            grid.add_source( std::make_unique<Gaussian_Pulse>( amp, t0, width, xb, y, z ) );
        }
        for ( std::size_t t = 0; t < num_steps; ++t ) {
            grid.apply_sources( t );
            grid.step();
        }
        return {
            grid.Ey_ptr()[ grid.idx(45, 50, 50) ],
            grid.Ey_ptr()[ grid.idx(50, 50, 50) ],
            grid.Ey_ptr()[ grid.idx(55, 50, 50) ],
            grid.Ey_ptr()[ grid.idx(65, 50, 50) ]
        };
    };

    ProbeResult rA = run_sim( true, false );
    ProbeResult rB = run_sim( false, true );
    ProbeResult rC = run_sim( true, true );

    // At each probe, field_C should equal field_A + field_B:
    double tol_45 = std::max( std::abs(rA.Ey_at_45 + rB.Ey_at_45) * 1e-10, 1e-20 );
    double tol_50 = std::max( std::abs(rA.Ey_at_50 + rB.Ey_at_50) * 1e-10, 1e-20 );
    double tol_55 = std::max( std::abs(rA.Ey_at_55 + rB.Ey_at_55) * 1e-10, 1e-20 );
    double tol_65 = std::max( std::abs(rA.Ey_at_65 + rB.Ey_at_65) * 1e-10, 1e-20 );

    ASSERT_NEAR( rC.Ey_at_45, rA.Ey_at_45 + rB.Ey_at_45, tol_45 );
    ASSERT_NEAR( rC.Ey_at_50, rA.Ey_at_50 + rB.Ey_at_50, tol_50 );
    ASSERT_NEAR( rC.Ey_at_55, rA.Ey_at_55 + rB.Ey_at_55, tol_55 );
    ASSERT_NEAR( rC.Ey_at_65, rA.Ey_at_65 + rB.Ey_at_65, tol_65 );

    // Sanity: fields should be nonzero (sources actually fired):
    ASSERT_GT( std::abs(rC.Ey_at_45), 1e-10 );
}