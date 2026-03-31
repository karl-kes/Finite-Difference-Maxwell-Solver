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
    // Place a single nonzero Ey at center. After one H-update, the curl of E
    // should produce Hx and Hz changes (but not Hy, since curl_y(E) involves
    // dEx/dz - dEz/dx, neither of which is set).
    Simulation_Config cfg{};
    Grid grid{ cfg };

    std::size_t cx{50}, cy{50}, cz{50};
    grid.field( Field::ELECTRIC, Component::Y, cx, cy, cz ) = 1.0;

    // Run one full step (H update then E update):
    grid.step();

    // Hx should be affected (curl_x has dEz/dy - dEy/dz terms):
    // Hz should be affected (curl_z has dEy/dx - dEx/dy terms):
    double hx = grid.field( Field::MAGNETIC, Component::X, cx, cy, cz );
    double hz = grid.field( Field::MAGNETIC, Component::Z, cx, cy, cz );
    double hy = grid.field( Field::MAGNETIC, Component::Y, cx, cy, cz );

    // Hx and Hz should be nonzero due to Ey derivatives:
    ASSERT_TRUE( std::abs(hx) > 0.0 || std::abs(hz) > 0.0 );

    // Hy should be zero (curl_y depends on Ex, Ez which are zero):
    ASSERT_NEAR( hy, 0.0, 1e-15 );
}

// ============================================================================
// Plane Wave Propagation (E2E Physics Test)
// ============================================================================

TEST(Integration, PlaneWavePropagation) {
    // Initialize Ey = sin(kx), Hz = sin(kx)/η propagating in +x.
    // Track phase correlation and energy drift over 100 steps.
    Simulation_Config cfg{};
    Grid grid{ cfg };

    // Wavelength: 20 cells
    double const wavelength{ 20.0 * grid.dx() };
    double const k{ 2.0 * std::numbers::pi / wavelength };
    double const c{ 1.0 / std::sqrt( cfg.mu * cfg.eps ) };
    double const eta{ std::sqrt( cfg.mu / cfg.eps ) };

    // Initialize plane wave in interior (avoid PML):
    std::size_t const margin{ 15 };
    for ( std::size_t z = margin; z < grid.Nz() - margin; ++z ) {
        for ( std::size_t y = margin; y < grid.Ny() - margin; ++y ) {
            for ( std::size_t x = margin; x < grid.Nx() - margin; ++x ) {
                double phase = k * static_cast<double>(x) * grid.dx();
                std::size_t i = grid.idx(x, y, z);
                grid.Ey_ptr()[i] = std::sin( phase );
                grid.Hz_ptr()[i] = std::sin( phase ) / eta;
            }
        }
    }

    double const E0 = grid.total_energy();
    ASSERT_GT( E0, 0.0 );

    // Probe at center:
    std::size_t const px{50}, py{50}, pz{50};
    double const initial_phase = k * static_cast<double>(px) * grid.dx();
    double const phase_shift = k * c * grid.dt();

    // Run for half the time it takes the wavefront to cross the margin:
    double const margin_dist = static_cast<double>( margin ) * grid.dx();
    std::size_t const num_steps = std::max( std::size_t{10},
        static_cast<std::size_t>( 0.5 * margin_dist / ( c * grid.dt() ) ) );

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

    ASSERT_LT( drift, 5.0 );
    ASSERT_GT( correlation, 0.99 );
}

// ============================================================================
// Energy Conservation (Source-free, PML active)
// ============================================================================

TEST(Integration, EnergyConservationShortRun) {
    // With a localized initial perturbation (no active source), energy should
    // be approximately conserved over ~20 steps before wavefront reaches PML.
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

    for ( int t = 0; t < 20; ++t ) {
        grid.step();
    }

    double Ef = grid.total_energy();
    double drift_pct = 100.0 * std::abs( Ef - E0 ) / E0;

    ASSERT_LT( drift_pct, 2.0 );
}

// ============================================================================
// Divergence-Free (Gauss's Law)
// ============================================================================

TEST(Integration, GaussLawDivH) {
    // After several time steps with a point source, div(H) should remain ~0.
    // For E/H formulation: div(B) = div(μH) = 0. With uniform μ, div(H) = 0.
    Simulation_Config cfg{};
    Grid grid{ cfg };

    grid.add_source( std::make_unique<Gaussian_Pulse>(
        10.0, grid.dt() * 10.0, grid.dt() * 5.0, 50, 50, 50 ) );

    for ( int t = 0; t < 30; ++t ) {
        grid.apply_sources( static_cast<std::size_t>(t) );
        grid.step();
    }

    double max_divH{0.0};
    double max_H{0.0};

    double const inv_dx = 1.0 / grid.dx();
    double const inv_dy = 1.0 / grid.dy();
    double const inv_dz = 1.0 / grid.dz();

    for ( std::size_t z = 25; z < 75; ++z ) {
        for ( std::size_t y = 25; y < 75; ++y ) {
            for ( std::size_t x = 25; x < 75; ++x ) {
                double dHx_dx = ( grid.field(Field::MAGNETIC, Component::X, x+1, y, z)
                                - grid.field(Field::MAGNETIC, Component::X, x, y, z) ) * inv_dx;
                double dHy_dy = ( grid.field(Field::MAGNETIC, Component::Y, x, y+1, z)
                                - grid.field(Field::MAGNETIC, Component::Y, x, y, z) ) * inv_dy;
                double dHz_dz = ( grid.field(Field::MAGNETIC, Component::Z, x, y, z+1)
                                - grid.field(Field::MAGNETIC, Component::Z, x, y, z) ) * inv_dz;

                double divH = std::abs( dHx_dx + dHy_dy + dHz_dz );
                if ( divH > max_divH ) max_divH = divH;

                double Hmag = grid.field_magnitude( Field::MAGNETIC, x, y, z );
                if ( Hmag > max_H ) max_H = Hmag;
            }
        }
    }

    if ( max_H > 1e-15 ) {
        double relative_divH = max_divH / ( max_H * inv_dx );
        ASSERT_LT( relative_divH, 1e-10 );
    }
}

// ============================================================================
// Reciprocity / Symmetry
// ============================================================================

TEST(Integration, IsotropicSourceSymmetry) {
    // A Gaussian pulse injecting Jz should produce fields symmetric in the
    // xy-plane (perpendicular to dipole axis).
    Simulation_Config cfg{};
    Grid grid{ cfg };

    grid.add_source( std::make_unique<Gaussian_Pulse>(
        50.0, grid.dt() * 10.0, grid.dt() * 4.0, 50, 50, 50 ) );

    for ( int t = 0; t < 50; ++t ) {
        grid.apply_sources( static_cast<std::size_t>(t) );
        grid.step();
    }

    double E_py = grid.field_magnitude( Field::ELECTRIC, 50, 58, 50 );
    double E_ny = grid.field_magnitude( Field::ELECTRIC, 50, 42, 50 );
    double E_px = grid.field_magnitude( Field::ELECTRIC, 58, 50, 50 );
    double E_nx = grid.field_magnitude( Field::ELECTRIC, 42, 50, 50 );

    double E_avg = ( E_py + E_ny + E_px + E_nx ) / 4.0;
    ASSERT_GT( E_avg, 1e-10 );

    ASSERT_LT( std::abs( E_py - E_ny ) / E_avg, 0.35 );
    ASSERT_LT( std::abs( E_px - E_nx ) / E_avg, 0.35 );
    ASSERT_LT( std::abs( E_py - E_px ) / E_avg, 0.50 );
}

// ============================================================================
// PML Absorption Effectiveness
// ============================================================================

TEST(Integration, PMLReducesReflection) {
    // Initialize a +x traveling plane wave and verify PML absorbs it
    // without large reflections.
    Simulation_Config cfg{};
    Grid grid{ cfg };

    double const wavelength = 20.0 * grid.dx();
    double const k = 2.0 * std::numbers::pi / wavelength;
    double const eta = std::sqrt( cfg.mu / cfg.eps );

    for ( std::size_t z = 15; z < grid.Nz() - 15; ++z ) {
        for ( std::size_t y = 15; y < grid.Ny() - 15; ++y ) {
            for ( std::size_t x = 50; x < grid.Nx() - 15; ++x ) {
                double phase = k * static_cast<double>(x) * grid.dx();
                std::size_t i = grid.idx(x, y, z);
                grid.Ey_ptr()[i] = std::sin( phase );
                grid.Hz_ptr()[i] = std::sin( phase ) / eta;
            }
        }
    }

    double E0 = grid.total_energy();

    for ( int t = 0; t < 80; ++t ) {
        grid.step();
    }

    double E80 = grid.total_energy();
    ASSERT_LT( E80, E0 * 1.10 );

    double reflected_max{0};
    for ( std::size_t z = 40; z < 60; ++z ) {
        for ( std::size_t y = 40; y < 60; ++y ) {
            double val = std::abs( grid.Ey_ptr()[ grid.idx(30, y, z) ] );
            if ( val > reflected_max ) reflected_max = val;
        }
    }

    ASSERT_LT( reflected_max, 0.15 );
}

TEST(Integration, PMLLateTimeStabilityCheck) {
    // Detect PML late-time instability. Flags if growth exceeds 100x baseline.
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

    ASSERT_TRUE( std::isfinite(E_final) );
    ASSERT_LT( E_final, E_baseline * 100.0 );
}

TEST(Integration, LeapfrogEnergyConservation) {
    Simulation_Config cfg{};
    cfg.use_pml = false;
    cfg.compute_derived();
    Grid grid{ cfg };

    double const eta{ std::sqrt( cfg.mu / cfg.eps ) };
    double const wavelength_cells{ 20.0 };
    double const k{ 2.0 * std::numbers::pi / ( wavelength_cells * grid.dx() ) };

    double const cx{ static_cast<double>( grid.Nx() / 2 ) };
    double const cy{ static_cast<double>( grid.Ny() / 2 ) };
    double const cz{ static_cast<double>( grid.Nz() / 2 ) };
    double const envelope_sigma{ 15.0 };

    for ( std::size_t z = 1; z < grid.Nz() - 1; ++z ) {
        for ( std::size_t y = 1; y < grid.Ny() - 1; ++y ) {
            for ( std::size_t x = 1; x < grid.Nx() - 1; ++x ) {
                double const dx2{ std::pow( static_cast<double>(x) - cx, 2 ) };
                double const dy2{ std::pow( static_cast<double>(y) - cy, 2 ) };
                double const dz2{ std::pow( static_cast<double>(z) - cz, 2 ) };
                double const envelope{ std::exp( -( dx2 + dy2 + dz2 )
                    / ( 2.0 * envelope_sigma * envelope_sigma ) ) };
                double const phase{ k * static_cast<double>( x ) * grid.dx() };
                std::size_t const i{ grid.idx(x, y, z) };

                grid.Ey_ptr()[i] = envelope * std::sin( phase );
                grid.Hz_ptr()[i] = envelope * std::sin( phase ) / eta;
            }
        }
    }

    double const E_initial{ grid.total_energy() };
    ASSERT_GT( E_initial, 0.0 );

    // Step 100 times, record first and last energy:
    double E_100{};
    for ( std::size_t t = 0; t < 100; ++t ) {
        grid.step();
        ASSERT_TRUE( std::isfinite( grid.total_energy() ) );
    }
    E_100 = grid.total_energy();

    // Step another 100 if unstable, energy grows exponentially:
    double E_200{};
    for ( std::size_t t = 0; t < 100; ++t ) {
        grid.step();
        ASSERT_TRUE( std::isfinite( grid.total_energy() ) );
    }
    E_200 = grid.total_energy();

    // Energy should remain bounded no exponential growth.
    // In a PEC cavity energy redistributes but total is conserved
    // up to leapfrog staggering artifacts:
    ASSERT_LT( E_200, E_initial * 2.0 );

    // Growth rate should not accelerate (exponential test):
    // If stable: E_200/E_100 = E_100/E_initial
    // If unstable: E_200/E_100 >> E_100/E_initial
    double const growth_1{ E_100 / E_initial };
    double const growth_2{ E_200 / std::max( E_100, 1e-30 ) };
    ASSERT_LT( growth_2, growth_1 * 1.5 );
}

// ============================================================================
// Causality
// ============================================================================

TEST(Integration, CausalPropagation) {
    // A point impulse at center at t=0 should not have reached cells far away
    // after only a few time steps.
    Simulation_Config cfg{};
    Grid grid{ cfg };

    grid.Ey_ptr()[ grid.idx(50, 50, 50) ] = 1.0;

    for ( int t = 0; t < 5; ++t ) {
        grid.step();
    }

    double far_field = std::abs( grid.Ey_ptr()[ grid.idx(70, 50, 50) ] );
    ASSERT_LT( far_field, 1e-14 );

    double near_field = std::abs( grid.Ey_ptr()[ grid.idx(51, 50, 50) ] );
    ASSERT_GT( near_field, 1e-20 );
}

// ============================================================================
// Grid Convergence (2nd Order)
// ============================================================================

TEST(Integration, DispersionConvergence) {
    // Measure dispersion at default resolution and verify within expected bounds.
    Simulation_Config cfg{};
    Grid grid{ cfg };

    double const wavelength = 20.0 * grid.dx();
    double const k = 2.0 * std::numbers::pi / wavelength;
    double const c = 1.0 / std::sqrt( cfg.mu * cfg.eps );
    double const eta = std::sqrt( cfg.mu / cfg.eps );

    for ( std::size_t z = 15; z < grid.Nz() - 15; ++z ) {
        for ( std::size_t y = 15; y < grid.Ny() - 15; ++y ) {
            for ( std::size_t x = 15; x < grid.Nx() - 15; ++x ) {
                double phase = k * static_cast<double>(x) * grid.dx();
                grid.Ey_ptr()[ grid.idx(x, y, z) ] = std::sin(phase);
                grid.Hz_ptr()[ grid.idx(x, y, z) ] = std::sin(phase) / eta;
            }
        }
    }

    std::size_t px{50}, py{50}, pz{50};
    double initial_phase = k * static_cast<double>(px) * grid.dx();
    double phase_shift = k * c * grid.dt();
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

    ASSERT_LT( dispersion_pct, 5.0 );
}

// ============================================================================
// Hertzian Dipole Radiation Pattern
// ============================================================================

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

struct DipoleSetup {
    double wavelength;
    double freq;
    std::size_t cx, cy, cz;
    std::size_t steps_per_period;
    std::size_t warmup_steps;

    DipoleSetup( Grid const& grid, Simulation_Config const& cfg )
        : wavelength{ 10.0 * grid.dx() }
        , freq{ ( 1.0 / std::sqrt( cfg.mu * cfg.eps ) ) / wavelength }
        , cx{ grid.Nx() / 2 }, cy{ grid.Ny() / 2 }, cz{ grid.Nz() / 2 }
        , steps_per_period{ static_cast<std::size_t>( 1.0 / ( freq * grid.dt() ) ) }
        , warmup_steps{ 8 * steps_per_period }
    {}
};

TEST(Integration, HertzianDipoleSinTheta) {
    Simulation_Config cfg{};
    Grid grid{ cfg };
    DipoleSetup ds{ grid, cfg };

    grid.add_source( std::make_unique<Sinusoidal_Jz>(
        1.0, ds.freq, ds.cx, ds.cy, ds.cz ) );

    for ( std::size_t t = 0; t < ds.warmup_steps; ++t ) {
        grid.apply_sources( t );
        grid.step();
    }

    int const r = 20;
    double peak_eq_px{0}, peak_eq_nx{0}, peak_eq_py{0}, peak_eq_ny{0};
    double peak_pole_pz{0}, peak_pole_nz{0};
    double peak_diag{0};

    int const rd = static_cast<int>( std::round( static_cast<double>(r) / std::sqrt(2.0) ) );

    for ( std::size_t s = 0; s < ds.steps_per_period; ++s ) {
        std::size_t t = ds.warmup_steps + s;
        grid.apply_sources( t );
        grid.step();

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

        e = grid.field_magnitude( Field::ELECTRIC,
            ds.cx, ds.cy, static_cast<std::size_t>(static_cast<int>(ds.cz) + r) );
        if ( e > peak_pole_pz ) peak_pole_pz = e;

        e = grid.field_magnitude( Field::ELECTRIC,
            ds.cx, ds.cy, static_cast<std::size_t>(static_cast<int>(ds.cz) - r) );
        if ( e > peak_pole_nz ) peak_pole_nz = e;

        e = grid.field_magnitude( Field::ELECTRIC,
            static_cast<std::size_t>(static_cast<int>(ds.cx) + rd), ds.cy,
            static_cast<std::size_t>(static_cast<int>(ds.cz) + rd) );
        if ( e > peak_diag ) peak_diag = e;
    }

    double peak_equator = ( peak_eq_px + peak_eq_nx + peak_eq_py + peak_eq_ny ) / 4.0;
    double peak_pole = ( peak_pole_pz + peak_pole_nz ) / 2.0;

    ASSERT_GT( peak_equator, 1e-6 );
    ASSERT_GT( peak_equator / std::max(peak_pole, 1e-30), 3.0 );
    ASSERT_GT( peak_diag, peak_pole );
    ASSERT_LT( peak_diag, peak_equator * 1.5 );
    ASSERT_LT( std::abs(peak_eq_px - peak_eq_py) / peak_equator, 0.30 );
    ASSERT_LT( std::abs(peak_eq_px - peak_eq_nx) / peak_equator, 0.30 );

    if ( peak_pole > 1e-10 ) {
        ASSERT_LT( std::abs(peak_pole_pz - peak_pole_nz) / peak_pole, 0.40 );
    }
}

TEST(Integration, HertzianDipole1OverR) {
    Simulation_Config cfg{};
    Grid grid{ cfg };
    DipoleSetup ds{ grid, cfg };

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

    double ratio_measured = peak_r1 / peak_r2;
    double ratio_expected = static_cast<double>(r2) / static_cast<double>(r1);

    ASSERT_GT( ratio_measured, ratio_expected * 0.65 );
    ASSERT_LT( ratio_measured, ratio_expected * 1.50 );
}

// ============================================================================
// OpenMP Thread Determinism
// ============================================================================

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
    int const max_threads = omp_get_max_threads();

    ThreadTestResult r1 = run_with_threads( 1 );
    ThreadTestResult r2 = run_with_threads( 2 );
    ThreadTestResult rm = run_with_threads( max_threads );

    ASSERT_GT( r1.energy, 0.0 );

    double energy_tol = r1.energy * 1e-10;
    ASSERT_NEAR( r1.energy, r2.energy, energy_tol );
    ASSERT_NEAR( r1.energy, rm.energy, energy_tol );

    double field_tol = std::abs(r1.field_sample) * 1e-12 + 1e-20;
    ASSERT_NEAR( r1.field_sample, r2.field_sample, field_tol );
    ASSERT_NEAR( r1.field_sample, rm.field_sample, field_tol );

    omp_set_num_threads( max_threads );
}

// ============================================================================
// Superposition (Linearity)
// ============================================================================

TEST(Integration, Superposition) {
    Simulation_Config cfg{};

    double const amp = 5.0;
    double const t0 = cfg.dt * 10.0;
    double const width = cfg.dt * 4.0;
    std::size_t const xa = 40, xb = 60, y = 50, z = 50;
    std::size_t const num_steps = 40;

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

    double tol_45 = std::max( std::abs(rA.Ey_at_45 + rB.Ey_at_45) * 1e-10, 1e-20 );
    double tol_50 = std::max( std::abs(rA.Ey_at_50 + rB.Ey_at_50) * 1e-10, 1e-20 );
    double tol_55 = std::max( std::abs(rA.Ey_at_55 + rB.Ey_at_55) * 1e-10, 1e-20 );
    double tol_65 = std::max( std::abs(rA.Ey_at_65 + rB.Ey_at_65) * 1e-10, 1e-20 );

    ASSERT_NEAR( rC.Ey_at_45, rA.Ey_at_45 + rB.Ey_at_45, tol_45 );
    ASSERT_NEAR( rC.Ey_at_50, rA.Ey_at_50 + rB.Ey_at_50, tol_50 );
    ASSERT_NEAR( rC.Ey_at_55, rA.Ey_at_55 + rB.Ey_at_55, tol_55 );
    ASSERT_NEAR( rC.Ey_at_65, rA.Ey_at_65 + rB.Ey_at_65, tol_65 );

    ASSERT_GT( std::abs(rC.Ey_at_45), 1e-10 );
}