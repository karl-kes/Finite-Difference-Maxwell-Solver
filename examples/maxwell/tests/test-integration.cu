#include "test-framework.hpp"
#include "test-helpers.hpp"

#include "../src/classes/config/config.hpp"
#include "../src/classes/grid/grid.hpp"
#include "../src/classes/source/source.hpp"

#include <cmath>
#include <numbers>
#include <vector>
#include <numeric>
#include <omp.h>

// Leapfrog Update Correctness

TEST(Integration, SingleStepCurlAnalytical) {
    // Place Ey=1 at center, step once, verify H at 4 grid points against
    // exact analytical single-step values: H -= (dt/mu) * curl(E).
    Simulation_Config cfg{};
    Grid grid{ cfg };
    std::size_t cx{50}, cy{50}, cz{50};
    grid.field( Field::ELECTRIC, Component::Y, cx, cy, cz ) = 1.0;
    double const Db{ cfg.dt / cfg.mu };
    double const inv_dx{ 1.0 / cfg.dx };
    double const inv_dz{ 1.0 / cfg.dz };

    grid.step();

    // Hx(cx,cy,cz): curl_x = -(0-1)/dz = +1/dz  =>  Hx = -Db/dz
    ASSERT_NEAR( grid.field(Field::MAGNETIC, Component::X, cx, cy, cz),
                 -Db * inv_dz, 1e-14 );
    // Hx(cx,cy,cz-1): curl_x = -(1-0)/dz = -1/dz  =>  Hx = +Db/dz
    ASSERT_NEAR( grid.field(Field::MAGNETIC, Component::X, cx, cy, cz-1),
                 Db * inv_dz, 1e-14 );
    // Hz(cx,cy,cz): curl_z = (0-1)/dx = -1/dx  =>  Hz = +Db/dx
    ASSERT_NEAR( grid.field(Field::MAGNETIC, Component::Z, cx, cy, cz),
                 Db * inv_dx, 1e-14 );
    // Hz(cx-1,cy,cz): curl_z = (1-0)/dx = +1/dx  =>  Hz = -Db/dx
    ASSERT_NEAR( grid.field(Field::MAGNETIC, Component::Z, cx-1, cy, cz),
                 -Db * inv_dx, 1e-14 );
    // Hy: curl_y depends on Ex, Ez which are zero:
    ASSERT_NEAR( grid.field(Field::MAGNETIC, Component::Y, cx, cy, cz), 0.0, 1e-15 );
}

// Plane Wave Propagation

TEST(Integration, PlaneWavePropagation) {
    Simulation_Config cfg{};
    Grid grid{ cfg };
    double const wavelength{ 20.0 * grid.dx() };
    double const k{ 2.0 * std::numbers::pi / wavelength };
    double const c{ 1.0 / std::sqrt( cfg.mu * cfg.eps ) };
    double const eta{ std::sqrt( cfg.mu / cfg.eps ) };
    std::size_t const margin{ 15 };
    for ( std::size_t z = margin; z < grid.Nz() - margin; ++z )
        for ( std::size_t y = margin; y < grid.Ny() - margin; ++y )
            for ( std::size_t x = margin; x < grid.Nx() - margin; ++x ) {
                double phase = k * static_cast<double>(x) * grid.dx();
                std::size_t i = grid.idx(x, y, z);
                grid.Ey_ptr()[i] = std::sin( phase );
                grid.Hz_ptr()[i] = std::sin( phase ) / eta;
            }

    double const E0 = grid.total_energy();
    ASSERT_GT( E0, 0.0 );
    std::size_t const px{50}, py{50}, pz{50};
    double const initial_phase = k * static_cast<double>(px) * grid.dx();
    double const phase_shift = k * c * grid.dt();
    double const margin_dist = static_cast<double>(margin) * grid.dx();
    std::size_t const num_steps = std::max( std::size_t{10},
        static_cast<std::size_t>( 0.5 * margin_dist / ( c * grid.dt() ) ) );
    double sum_exp{}, sum_act{}, sum_exp2{}, sum_act2{}, sum_prod{};
    for ( std::size_t t = 0; t < num_steps; ++t ) {
        double expected = std::sin( initial_phase - phase_shift * static_cast<double>(t) );
        double actual = grid.Ey_ptr()[ grid.idx(px, py, pz) ];
        sum_exp += expected; sum_act += actual;
        sum_exp2 += expected*expected; sum_act2 += actual*actual;
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

// Energy Conservation

TEST(Integration, EnergyConservationNoPML) {
    // PML OFF (PEC cavity). Gaussian E-field blob. No dissipation mechanism.
    // The Gaussian blob isn't a matched mode — it radiates, hits PEC walls,
    // and reflects. The leapfrog staggering artifact (E and H measured at
    // different half-steps) causes measured energy to oscillate. We verify
    // the oscillation stays bounded and doesn't blow up.
    Simulation_Config cfg{};
    cfg.use_pml = false;
    cfg.compute_derived();
    Grid grid{ cfg };
    double sigma{5.0};
    for ( std::size_t z = 30; z < 70; ++z )
        for ( std::size_t y = 30; y < 70; ++y )
            for ( std::size_t x = 30; x < 70; ++x ) {
                double r2 = std::pow(static_cast<double>(x) - 50.0, 2)
                          + std::pow(static_cast<double>(y) - 50.0, 2)
                          + std::pow(static_cast<double>(z) - 50.0, 2);
                grid.Ey_ptr()[ grid.idx(x, y, z) ] = std::exp( -r2 / (2.0 * sigma * sigma) );
            }
    double E0 = grid.total_energy();
    ASSERT_GT( E0, 0.0 );

    // Track max drift over 200 steps:
    double max_drift{ 0.0 };
    for ( int t = 0; t < 200; ++t ) {
        grid.step();
        double drift = std::abs( grid.total_energy() - E0 ) / E0;
        if ( drift > max_drift ) max_drift = drift;
    }

    // Leapfrog staggering artifact grows as the field becomes more spatially
    // complex (boundary reflections), but should stay bounded. Actual max
    // drift ~3% over 200 steps. True instability would give 100%+.
    ASSERT_LT( max_drift, 0.10 );
    ASSERT_TRUE( std::isfinite( grid.total_energy() ) );
}

// Divergence-Free (Gauss's Law)

TEST(Integration, GaussLawDivH) {
    Simulation_Config cfg{};
    Grid grid{ cfg };
    grid.add_source( std::make_unique<Gaussian_Pulse>(
        10.0, grid.dt() * 10.0, grid.dt() * 5.0, 50, 50, 50 ) );
    for ( int t = 0; t < 30; ++t ) {
        grid.apply_sources( static_cast<std::size_t>(t) );
        grid.step();
    }
    double max_divH{0.0}, max_H{0.0};
    double const inv_dx = 1.0 / grid.dx();
    double const inv_dy = 1.0 / grid.dy();
    double const inv_dz = 1.0 / grid.dz();
    for ( std::size_t z = 25; z < 75; ++z )
        for ( std::size_t y = 25; y < 75; ++y )
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
    if ( max_H > 1e-15 ) {
        double relative_divH = max_divH / ( max_H * inv_dx );
        ASSERT_LT( relative_divH, 1e-10 );
    }
}

// Symmetry

TEST(Integration, IsotropicSourceSymmetry) {
    // Gaussian Jz pulse: xy-plane symmetry. Actual asymmetries ~2.3%:
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
    ASSERT_LT( std::abs( E_py - E_ny ) / E_avg, 0.10 );
    ASSERT_LT( std::abs( E_px - E_nx ) / E_avg, 0.10 );
    ASSERT_LT( std::abs( E_py - E_px ) / E_avg, 0.10 );
}

// PML

TEST(Integration, PMLAbsorbsVsPEC) {
    // Compare total energy after a Gaussian E-field blob has expanded and
    // hit boundaries: PML absorbs outgoing energy, PEC reflects it.
    auto run = []( bool use_pml ) -> double {
        Simulation_Config cfg{};
        cfg.use_pml = use_pml;
        cfg.compute_derived();
        Grid grid{ cfg };
        // Gaussian blob centered in grid:
        double const sigma{ 8.0 };
        for ( std::size_t z = 20; z < grid.Nz()-20; ++z )
            for ( std::size_t y = 20; y < grid.Ny()-20; ++y )
                for ( std::size_t x = 20; x < grid.Nx()-20; ++x ) {
                    double r2 = std::pow(static_cast<double>(x)-50.0, 2)
                              + std::pow(static_cast<double>(y)-50.0, 2)
                              + std::pow(static_cast<double>(z)-50.0, 2);
                    grid.Ey_ptr()[ grid.idx(x,y,z) ] = std::exp(-r2/(2.0*sigma*sigma));
                }
        for ( int t = 0; t < 120; ++t ) grid.step();
        return grid.total_energy();
    };

    double E_pml = run( true );
    double E_pec = run( false );

    ASSERT_GT( E_pml, 0.0 );
    ASSERT_GT( E_pec, 0.0 );
    // PEC conserves; PML absorbs. PML energy should be well below PEC:
    ASSERT_LT( E_pml, E_pec * 0.50 );
}

TEST(Integration, PMLLateTimeStability) {
    // Detect exponential CPML instability. Use a 40^3 grid to afford 5000
    // steps. Compare growth rates in two halves — exponential instability
    // accelerates, while normal residual energy does not.
    Simulation_Config cfg{};
    cfg.Nx = 40; cfg.Ny = 40; cfg.Nz = 40;
    cfg.compute_derived();
    Grid grid{ cfg };
    std::size_t cx = grid.Nx()/2, cy = grid.Ny()/2, cz = grid.Nz()/2;
    grid.add_source( std::make_unique<Gaussian_Pulse>(
        10.0, grid.dt() * 10.0, grid.dt() * 4.0, cx, cy, cz ) );

    for ( std::size_t t = 0; t < 20; ++t ) {
        grid.apply_sources( t );
        grid.step();
    }
    double E_baseline = grid.total_energy();

    for ( std::size_t t = 20; t < 2520; ++t ) grid.step();
    double E_mid = grid.total_energy();

    for ( std::size_t t = 2520; t < 5020; ++t ) grid.step();
    double E_final = grid.total_energy();

    ASSERT_TRUE( std::isfinite(E_mid) );
    ASSERT_TRUE( std::isfinite(E_final) );

    // Exponential growth detection: if unstable, growth accelerates
    // (growth_2/growth_1 > 1 and increasing). Stable: decelerates.
    // Measured: g1~440, g2~3.8, ratio~0.009. Unstable would give ratio>1.
    // Threshold of 2.0 tolerates linear residual growth while catching
    // exponential divergence:
    double const growth_1{ E_mid / std::max( E_baseline, 1e-30 ) };
    double const growth_2{ E_final / std::max( E_mid, 1e-30 ) };
    ASSERT_LT( growth_2, growth_1 * 2.0 );
}

TEST(Integration, LeapfrogEnergyConservation) {
    // PML off. Gaussian-enveloped plane wave in PEC cavity.
    // No exponential growth allowed. Actual E_200/E_initial ~1.018.
    Simulation_Config cfg{};
    cfg.use_pml = false;
    cfg.compute_derived();
    Grid grid{ cfg };
    double const eta{ std::sqrt( cfg.mu / cfg.eps ) };
    double const k{ 2.0 * std::numbers::pi / ( 20.0 * grid.dx() ) };
    double const cx{ static_cast<double>( grid.Nx() / 2 ) };
    double const cy{ static_cast<double>( grid.Ny() / 2 ) };
    double const cz{ static_cast<double>( grid.Nz() / 2 ) };
    double const envelope_sigma{ 15.0 };
    for ( std::size_t z = 1; z < grid.Nz() - 1; ++z )
        for ( std::size_t y = 1; y < grid.Ny() - 1; ++y )
            for ( std::size_t x = 1; x < grid.Nx() - 1; ++x ) {
                double dx2 = std::pow( static_cast<double>(x) - cx, 2 );
                double dy2 = std::pow( static_cast<double>(y) - cy, 2 );
                double dz2 = std::pow( static_cast<double>(z) - cz, 2 );
                double env = std::exp( -(dx2+dy2+dz2) / (2.0*envelope_sigma*envelope_sigma) );
                double phase = k * static_cast<double>(x) * grid.dx();
                std::size_t i = grid.idx(x, y, z);
                grid.Ey_ptr()[i] = env * std::sin(phase);
                grid.Hz_ptr()[i] = env * std::sin(phase) / eta;
            }
    double const E_initial{ grid.total_energy() };
    ASSERT_GT( E_initial, 0.0 );
    for ( std::size_t t = 0; t < 100; ++t ) {
        grid.step();
        ASSERT_TRUE( std::isfinite( grid.total_energy() ) );
    }
    double E_100 = grid.total_energy();
    for ( std::size_t t = 0; t < 100; ++t ) {
        grid.step();
        ASSERT_TRUE( std::isfinite( grid.total_energy() ) );
    }
    double E_200 = grid.total_energy();
    ASSERT_LT( E_200, E_initial * 1.10 );
    double const growth_1{ E_100 / E_initial };
    double const growth_2{ E_200 / std::max( E_100, 1e-30 ) };
    ASSERT_LT( growth_2, growth_1 * 1.5 );
}

// Causality

TEST(Integration, CausalPropagation) {
    Simulation_Config cfg{};
    Grid grid{ cfg };
    grid.Ey_ptr()[ grid.idx(50, 50, 50) ] = 1.0;
    for ( int t = 0; t < 5; ++t ) grid.step();
    ASSERT_LT( std::abs( grid.Ey_ptr()[ grid.idx(70, 50, 50) ] ), 1e-14 );
    ASSERT_GT( std::abs( grid.Ey_ptr()[ grid.idx(51, 50, 50) ] ), 1e-20 );
}

// Dispersion

TEST(Integration, DispersionConvergence) {
    Simulation_Config cfg{};
    Grid grid{ cfg };
    double const wavelength = 20.0 * grid.dx();
    double const k = 2.0 * std::numbers::pi / wavelength;
    double const c = 1.0 / std::sqrt( cfg.mu * cfg.eps );
    double const eta = std::sqrt( cfg.mu / cfg.eps );
    for ( std::size_t z = 15; z < grid.Nz() - 15; ++z )
        for ( std::size_t y = 15; y < grid.Ny() - 15; ++y )
            for ( std::size_t x = 15; x < grid.Nx() - 15; ++x ) {
                double phase = k * static_cast<double>(x) * grid.dx();
                grid.Ey_ptr()[ grid.idx(x, y, z) ] = std::sin(phase);
                grid.Hz_ptr()[ grid.idx(x, y, z) ] = std::sin(phase) / eta;
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
    double dispersion_pct = 100.0 * (total_phase_err / 100.0) / std::max(phase_shift * 100.0, 1e-10);
    ASSERT_LT( dispersion_pct, 5.0 );
}

// Hertzian Dipole

class Sinusoidal_Jz : public Source {
    double amplitude_, frequency_;
    std::size_t x_, y_, z_;
public:
    Sinusoidal_Jz( double amp, double freq, std::size_t x, std::size_t y, std::size_t z )
        : amplitude_{amp}, frequency_{freq}, x_{x}, y_{y}, z_{z} {}
    void apply( Grid &grid, std::size_t const time_step ) override {
        double const omega = 2.0 * std::numbers::pi * frequency_;
        double const t = static_cast<double>(time_step) * grid.dt();
        grid.Jz_ptr()[ grid.idx(x_, y_, z_) ] = amplitude_ * std::sin(omega * t);
    }
};

struct DipoleSetup {
    double wavelength, freq;
    std::size_t cx, cy, cz, steps_per_period, warmup_steps;
    DipoleSetup( Grid const& grid, Simulation_Config const& cfg )
        : wavelength{10.0 * grid.dx()}
        , freq{ (1.0/std::sqrt(cfg.mu*cfg.eps)) / wavelength }
        , cx{grid.Nx()/2}, cy{grid.Ny()/2}, cz{grid.Nz()/2}
        , steps_per_period{ static_cast<std::size_t>(1.0/(freq*grid.dt())) }
        , warmup_steps{8*steps_per_period} {}
};

TEST(Integration, HertzianDipoleSinTheta) {
    Simulation_Config cfg{};
    Grid grid{ cfg };
    DipoleSetup ds{ grid, cfg };
    grid.add_source( std::make_unique<Sinusoidal_Jz>(1.0, ds.freq, ds.cx, ds.cy, ds.cz) );
    for ( std::size_t t = 0; t < ds.warmup_steps; ++t ) { grid.apply_sources(t); grid.step(); }

    int const r = 20;
    double peak_eq_px{0}, peak_eq_nx{0}, peak_eq_py{0}, peak_eq_ny{0};
    double peak_pole_pz{0}, peak_pole_nz{0}, peak_diag{0};
    int const rd = static_cast<int>( std::round( static_cast<double>(r)/std::sqrt(2.0) ) );

    for ( std::size_t s = 0; s < ds.steps_per_period; ++s ) {
        std::size_t t = ds.warmup_steps + s;
        grid.apply_sources(t); grid.step();
        auto mag = [&](std::size_t x, std::size_t y, std::size_t z) {
            return grid.field_magnitude(Field::ELECTRIC, x, y, z); };
        double e;
        e = mag(static_cast<std::size_t>(static_cast<int>(ds.cx)+r), ds.cy, ds.cz); if (e>peak_eq_px) peak_eq_px=e;
        e = mag(static_cast<std::size_t>(static_cast<int>(ds.cx)-r), ds.cy, ds.cz); if (e>peak_eq_nx) peak_eq_nx=e;
        e = mag(ds.cx, static_cast<std::size_t>(static_cast<int>(ds.cy)+r), ds.cz); if (e>peak_eq_py) peak_eq_py=e;
        e = mag(ds.cx, static_cast<std::size_t>(static_cast<int>(ds.cy)-r), ds.cz); if (e>peak_eq_ny) peak_eq_ny=e;
        e = mag(ds.cx, ds.cy, static_cast<std::size_t>(static_cast<int>(ds.cz)+r)); if (e>peak_pole_pz) peak_pole_pz=e;
        e = mag(ds.cx, ds.cy, static_cast<std::size_t>(static_cast<int>(ds.cz)-r)); if (e>peak_pole_nz) peak_pole_nz=e;
        e = mag(static_cast<std::size_t>(static_cast<int>(ds.cx)+rd), ds.cy,
                static_cast<std::size_t>(static_cast<int>(ds.cz)+rd)); if (e>peak_diag) peak_diag=e;
    }
    double peak_equator = (peak_eq_px+peak_eq_nx+peak_eq_py+peak_eq_ny)/4.0;
    double peak_pole = (peak_pole_pz+peak_pole_nz)/2.0;

    ASSERT_GT( peak_equator, 1e-6 );
    // Equator >> pole (sin(theta) pattern). Actual ratio ~5.5:
    // Equator >> pole (sin(theta) pattern). Actual ratio ~5.5:
    ASSERT_GT( peak_equator / std::max(peak_pole, 1e-30), 4.0 );
    ASSERT_GT( peak_diag, peak_pole );
    ASSERT_LT( peak_diag, peak_equator * 1.2 );
    // Equatorial symmetry — actual ~0.16%. Threshold 5%:
    ASSERT_LT( std::abs(peak_eq_px - peak_eq_py) / peak_equator, 0.05 );
    ASSERT_LT( std::abs(peak_eq_px - peak_eq_nx) / peak_equator, 0.05 );
    if ( peak_pole > 1e-10 )
        ASSERT_LT( std::abs(peak_pole_pz - peak_pole_nz) / peak_pole, 0.10 );
}

TEST(Integration, HertzianDipole1OverR) {
    Simulation_Config cfg{};
    Grid grid{ cfg };
    DipoleSetup ds{ grid, cfg };
    grid.add_source( std::make_unique<Sinusoidal_Jz>(1.0, ds.freq, ds.cx, ds.cy, ds.cz) );
    for ( std::size_t t = 0; t < ds.warmup_steps; ++t ) { grid.apply_sources(t); grid.step(); }
    int const r1 = 15, r2 = 25;
    double peak_r1{0}, peak_r2{0};
    for ( std::size_t s = 0; s < ds.steps_per_period; ++s ) {
        std::size_t t = ds.warmup_steps + s;
        grid.apply_sources(t); grid.step();
        double e1 = grid.field_magnitude(Field::ELECTRIC,
            static_cast<std::size_t>(static_cast<int>(ds.cx)+r1), ds.cy, ds.cz);
        double e2 = grid.field_magnitude(Field::ELECTRIC,
            static_cast<std::size_t>(static_cast<int>(ds.cx)+r2), ds.cy, ds.cz);
        if (e1>peak_r1) peak_r1=e1;
        if (e2>peak_r2) peak_r2=e2;
    }
    ASSERT_GT( peak_r1, 1e-6 );
    ASSERT_GT( peak_r2, 1e-6 );
    double ratio_measured = peak_r1 / peak_r2;
    double ratio_expected = static_cast<double>(r2) / static_cast<double>(r1);
    // Actual ratio within 0.01% of expected. Threshold ±5%:
    ASSERT_GT( ratio_measured, ratio_expected * 0.95 );
    ASSERT_LT( ratio_measured, ratio_expected * 1.05 );
}

// Thread Determinism

struct ThreadTestResult { double energy; double field_sample; };

static ThreadTestResult run_with_threads( int num_threads ) {
    omp_set_num_threads( num_threads );
    Simulation_Config cfg{};
    Grid grid{ cfg };
    grid.add_source( std::make_unique<Gaussian_Pulse>(
        10.0, grid.dt() * 10.0, grid.dt() * 4.0, 50, 50, 50 ) );
    for ( std::size_t t = 0; t < 50; ++t ) { grid.apply_sources(t); grid.step(); }
    return { grid.total_energy(), grid.Ey_ptr()[ grid.idx(60, 50, 50) ] };
}

TEST(Integration, ThreadDeterminism) {
    int const max_threads = omp_get_max_threads();
    ThreadTestResult r1 = run_with_threads(1);
    ThreadTestResult r2 = run_with_threads(2);
    ThreadTestResult rm = run_with_threads(max_threads);
    ASSERT_GT( r1.energy, 0.0 );
    double energy_tol = r1.energy * 1e-10;
    ASSERT_NEAR( r1.energy, r2.energy, energy_tol );
    ASSERT_NEAR( r1.energy, rm.energy, energy_tol );
    double field_tol = std::abs(r1.field_sample) * 1e-12 + 1e-20;
    ASSERT_NEAR( r1.field_sample, r2.field_sample, field_tol );
    ASSERT_NEAR( r1.field_sample, rm.field_sample, field_tol );
    omp_set_num_threads( max_threads );
}

// Superposition

TEST(Integration, Superposition) {
    Simulation_Config cfg{};
    double const amp = 5.0, t0 = cfg.dt * 10.0, width = cfg.dt * 4.0;
    std::size_t const xa = 40, xb = 60, y = 50, z = 50, num_steps = 40;
    struct P { double a, b, c, d; };
    auto run_sim = [&]( bool A, bool B ) -> P {
        Grid grid{ cfg };
        if (A) grid.add_source(std::make_unique<Gaussian_Pulse>(amp, t0, width, xa, y, z));
        if (B) grid.add_source(std::make_unique<Gaussian_Pulse>(amp, t0, width, xb, y, z));
        for ( std::size_t t = 0; t < num_steps; ++t ) { grid.apply_sources(t); grid.step(); }
        return { grid.Ey_ptr()[grid.idx(45,50,50)], grid.Ey_ptr()[grid.idx(50,50,50)],
                 grid.Ey_ptr()[grid.idx(55,50,50)], grid.Ey_ptr()[grid.idx(65,50,50)] };
    };
    P rA = run_sim(true, false);
    P rB = run_sim(false, true);
    P rC = run_sim(true, true);
    auto tol = [](double v){ return std::max(std::abs(v)*1e-10, 1e-20); };
    ASSERT_NEAR( rC.a, rA.a + rB.a, tol(rA.a+rB.a) );
    ASSERT_NEAR( rC.b, rA.b + rB.b, tol(rA.b+rB.b) );
    ASSERT_NEAR( rC.c, rA.c + rB.c, tol(rA.c+rB.c) );
    ASSERT_NEAR( rC.d, rA.d + rB.d, tol(rA.d+rB.d) );
    ASSERT_GT( std::abs(rC.a), 1e-10 );
}