#include "test_framework.hpp"
#include "test_helpers.hpp"
#include "../src/Classes/Config/config.hpp"
#include "../src/Classes/Grid/grid.hpp"
#include "../src/Classes/Source/source.hpp"
#include "../src/Classes/Write_Output/output.hpp"
#include "../src/Classes/Validation/validation.hpp"
#include <cmath>
#include <numbers>
#include <fstream>
#include <cstdint>
#include <vector>
#include <filesystem>

// ============================================================================
// Output Binary Round-Trip
// ============================================================================

TEST(Output, BinaryRoundTrip) {
    Simulation_Config cfg{};
    Grid grid{ cfg };

    for ( std::size_t z = 0; z < grid.Nz(); ++z ) {
        for ( std::size_t y = 0; y < grid.Ny(); ++y ) {
            for ( std::size_t x = 0; x < grid.Nx(); ++x ) {
                grid.Ey_ptr()[ grid.idx(x, y, z) ] = static_cast<double>(x);
            }
        }
    }

    std::string test_dir{ "test_output_roundtrip" };
    Output output{ test_dir };
    output.initialize( grid );
    output.write_field( grid );
    output.finalize();

    std::string path{ test_dir + "/E.bin" };
    std::ifstream file( path, std::ios::binary );
    ASSERT_TRUE( file.is_open() );

    uint64_t dims[3];
    file.read( reinterpret_cast<char*>(dims), sizeof(dims) );
    std::size_t nx = static_cast<std::size_t>(dims[0]);
    std::size_t ny = static_cast<std::size_t>(dims[1]);
    std::size_t nz = static_cast<std::size_t>(dims[2]);

    ASSERT_EQ( nx, grid.Nx() - 1 );
    ASSERT_EQ( ny, grid.Ny() - 1 );
    ASSERT_EQ( nz, grid.Nz() - 1 );

    std::size_t total_floats = nx * ny * nz * 3;
    std::vector<double> data( total_floats );
    file.read( reinterpret_cast<char*>(data.data()),
               static_cast<std::streamsize>(total_floats * sizeof(double)) );
    ASSERT_TRUE( file.good() );
    file.close();

    for ( std::size_t z = 10; z < 20; ++z ) {
        for ( std::size_t y = 10; y < 20; ++y ) {
            for ( std::size_t x = 10; x < 20; ++x ) {
                std::size_t offset = ( z * ny * nx + y * nx + x ) * 3;
                double Ex_avg = data[offset + 0];
                double Ey_avg = data[offset + 1];
                double Ez_avg = data[offset + 2];

                ASSERT_NEAR( Ex_avg, 0.0, 1e-14 );
                ASSERT_NEAR( Ez_avg, 0.0, 1e-14 );

                double expected_Ey = static_cast<double>(x);
                ASSERT_NEAR( Ey_avg, expected_Ey, 1e-12 );
            }
        }
    }

    ASSERT_TRUE( std::filesystem::exists( test_dir + "/H.bin" ) );
    std::filesystem::remove_all( test_dir );
}

TEST(Output, HeaderDimensionsConsistent) {
    Simulation_Config cfg{};
    Grid grid{ cfg };

    std::string test_dir{ "test_output_header" };
    Output output{ test_dir };
    output.initialize( grid );

    grid.Ey_ptr()[ grid.idx(50, 50, 50) ] = 1.0;
    output.write_field( grid );
    output.finalize();

    std::size_t const frame_bytes{ (grid.Nx() - 1) * (grid.Ny() - 1) * (grid.Nz() - 1) * 3 * sizeof(double) };

    for ( char const* name : { "/E.bin", "/H.bin" } ) {
        std::string path{ test_dir + name };
        std::ifstream file( path, std::ios::binary );
        ASSERT_TRUE( file.is_open() );

        uint64_t dims[3];
        file.read( reinterpret_cast<char*>(dims), sizeof(dims) );

        ASSERT_EQ( static_cast<std::size_t>(dims[0]), grid.Nx() - 1 );
        ASSERT_EQ( static_cast<std::size_t>(dims[1]), grid.Ny() - 1 );
        ASSERT_EQ( static_cast<std::size_t>(dims[2]), grid.Nz() - 1 );

        file.seekg( 0, std::ios::end );
        auto file_size = file.tellg();
        auto expected_size = static_cast<std::streampos>( 24 + frame_bytes );
        ASSERT_EQ( file_size, expected_size );

        file.close();
    }

    std::filesystem::remove_all( test_dir );
}

// ============================================================================
// Validation Class
// ============================================================================

TEST(ValidationClass, MetricsPhysicallyReasonable) {
    Simulation_Config cfg{};
    Plane_Wave_Test test{ cfg };
    Validation_Result result = test.run();

    ASSERT_TRUE( result.passed );

    // Energy drift present but small (actual ~3.8%):
    ASSERT_GT( result.energy_drift_percent, 0.0 );
    ASSERT_LT( result.energy_drift_percent, 5.0 );

    // Phase correlation very high (actual ~0.999998):
    ASSERT_GT( result.phase_correlation, 0.999 );

    // Dispersion small (actual ~0.33%):
    ASSERT_GT( result.dispersion_percent, 0.0 );
    ASSERT_LT( result.dispersion_percent, 2.0 );
}

// ============================================================================
// Quantitative CPML Reflection Coefficient
// ============================================================================

TEST(Integration, CPMLReflectionCoefficient) {
    Simulation_Config cfg{};
    Grid grid{ cfg };

    std::size_t const pml_t{ cfg.pml_thickness };
    std::size_t const nx{ grid.Nx() };
    std::size_t const ny{ grid.Ny() };
    std::size_t const nz{ grid.Nz() };

    double const wavelength = 20.0 * grid.dx();
    double const k = 2.0 * std::numbers::pi / wavelength;
    double const c = 1.0 / std::sqrt( cfg.mu * cfg.eps );
    double const eta = std::sqrt( cfg.mu / cfg.eps );
    double const incident_amplitude = 1.0;

    std::size_t const yz_margin{ pml_t + 10 };

    std::size_t const pml_inner_hi{ nx - pml_t };
    std::size_t const pml_inner_lo{ pml_t };

    std::size_t const interior_lo{ pml_inner_lo + 2 };
    std::size_t const interior_hi{ pml_inner_hi - 1 };
    std::size_t const interior_width{ interior_hi - interior_lo };

    std::size_t const wave_start{ interior_lo + interior_width * 2 / 5 };
    std::size_t const wave_end{ interior_hi };

    std::size_t const taper_cells{ std::min(
        static_cast<std::size_t>( 2.0 * wavelength / grid.dx() ),
        ( wave_end - wave_start ) / 4 ) };

    std::size_t const probe_x{ pml_inner_lo + ( wave_start - pml_inner_lo ) / 2 };

    std::size_t const probe_y_lo{ ny * 2 / 5 };
    std::size_t const probe_y_hi{ ny * 3 / 5 };
    std::size_t const probe_z_lo{ nz * 2 / 5 };
    std::size_t const probe_z_hi{ nz * 3 / 5 };

    for ( std::size_t z = yz_margin; z < nz - yz_margin; ++z ) {
        for ( std::size_t y = yz_margin; y < ny - yz_margin; ++y ) {
            for ( std::size_t x = wave_start; x < wave_end; ++x ) {
                double envelope = 1.0;
                if ( x < wave_start + taper_cells ) {
                    double t = static_cast<double>( x - wave_start ) / static_cast<double>( taper_cells );
                    envelope = 0.5 * ( 1.0 - std::cos( std::numbers::pi * t ) );
                } else if ( x >= wave_end - taper_cells ) {
                    double t = static_cast<double>( wave_end - 1 - x ) / static_cast<double>( taper_cells );
                    envelope = 0.5 * ( 1.0 - std::cos( std::numbers::pi * t ) );
                }

                double phase = k * static_cast<double>( x ) * grid.dx();
                std::size_t i = grid.idx( x, y, z );
                grid.Ey_ptr()[i] = incident_amplitude * envelope * std::sin( phase );
                grid.Hz_ptr()[i] = incident_amplitude * envelope * std::sin( phase ) / eta;
            }
        }
    }

    double const dist_to_pml{ static_cast<double>( wave_end - wave_start ) * grid.dx() };
    double const dist_return{ static_cast<double>( wave_start - probe_x ) * grid.dx() };
    double const travel_time{ ( dist_to_pml + dist_return ) / c };
    std::size_t const num_steps{ static_cast<std::size_t>(
        std::ceil( 1.5 * travel_time / grid.dt() ) ) };

    for ( std::size_t t = 0; t < num_steps; ++t ) {
        grid.step();
    }

    double peak_reflected{ 0 };
    for ( std::size_t z = probe_z_lo; z < probe_z_hi; ++z ) {
        for ( std::size_t y = probe_y_lo; y < probe_y_hi; ++y ) {
            double val = std::abs( grid.Ey_ptr()[ grid.idx( probe_x, y, z ) ] );
            if ( val > peak_reflected ) peak_reflected = val;
        }
    }

    double R_dB = 20.0 * std::log10( std::max( peak_reflected, 1e-30 ) / incident_amplitude );

    double const threshold_dB{ std::min( -12.0,
        -0.8 * static_cast<double>( pml_t ) - 8.0 ) };
    ASSERT_LT( R_dB, threshold_dB );

    ASSERT_GT( peak_reflected, 1e-20 );
}

// ============================================================================
// Multi-Axis CPML Reflection
// ============================================================================

static double cpml_reflection_on_axis( int axis ) {
    // Launch a tapered plane wave along the given axis (0=x, 1=y, 2=z),
    // wait for the round-trip to a probe behind the wave start, and
    // measure the peak reflected amplitude in dB.
    Simulation_Config cfg{};
    cfg.Nx = 60; cfg.Ny = 60; cfg.Nz = 60;
    cfg.compute_derived();
    Grid grid{ cfg };

    std::size_t const pml_t{ cfg.pml_thickness };
    std::size_t const n[3]{ grid.Nx(), grid.Ny(), grid.Nz() };
    double const d[3]{ cfg.dx, cfg.dy, cfg.dz };
    double const c{ 1.0 / std::sqrt( cfg.mu * cfg.eps ) };
    double const eta{ std::sqrt( cfg.mu / cfg.eps ) };

    double const wavelength{ 15.0 * d[axis] };
    double const k{ 2.0 * std::numbers::pi / wavelength };

    std::size_t const wave_start{ pml_t + 5 };
    std::size_t const wave_end{ n[axis] - pml_t - 1 };
    std::size_t const taper{ std::min( std::size_t{10}, ( wave_end - wave_start ) / 4 ) };
    std::size_t const probe_pos{ pml_t + 2 };

    for ( std::size_t z = 1; z < n[2] - 1; ++z )
        for ( std::size_t y = 1; y < n[1] - 1; ++y )
            for ( std::size_t x = 1; x < n[0] - 1; ++x ) {
                std::size_t const coords[3]{ x, y, z };
                std::size_t const ap{ coords[axis] };
                if ( ap < wave_start || ap >= wave_end ) continue;

                double envelope{ 1.0 };
                if ( ap < wave_start + taper ) {
                    double t = static_cast<double>( ap - wave_start ) / static_cast<double>( taper );
                    envelope = 0.5 * ( 1.0 - std::cos( std::numbers::pi * t ) );
                } else if ( ap >= wave_end - taper ) {
                    double t = static_cast<double>( wave_end - 1 - ap ) / static_cast<double>( taper );
                    envelope = 0.5 * ( 1.0 - std::cos( std::numbers::pi * t ) );
                }

                double const phase{ k * static_cast<double>( ap ) * d[axis] };
                std::size_t const i{ grid.idx( x, y, z ) };
                double const val{ envelope * std::sin( phase ) };

                // Plane wave polarization: axis 0→(Ey,Hz), 1→(Ez,Hx), 2→(Ex,Hy)
                if ( axis == 0 ) { grid.Ey_ptr()[i] = val; grid.Hz_ptr()[i] = val / eta; }
                else if ( axis == 1 ) { grid.Ez_ptr()[i] = val; grid.Hx_ptr()[i] = val / eta; }
                else { grid.Ex_ptr()[i] = val; grid.Hy_ptr()[i] = val / eta; }
            }

    double const dist_fwd{ static_cast<double>( wave_end - wave_start ) * d[axis] };
    double const dist_ret{ static_cast<double>( wave_start - probe_pos ) * d[axis] };
    std::size_t const steps{ static_cast<std::size_t>(
        std::ceil( 1.5 * ( dist_fwd + dist_ret ) / ( c * cfg.dt ) ) ) };

    for ( std::size_t t = 0; t < steps; ++t ) grid.step();

    double peak{ 0 };
    for ( std::size_t a = n[axis] / 3; a < 2 * n[axis] / 3; ++a )
        for ( std::size_t b = n[axis] / 3; b < 2 * n[axis] / 3; ++b ) {
            std::size_t x, y, z;
            if ( axis == 0 ) { x = probe_pos; y = a; z = b; }
            else if ( axis == 1 ) { x = a; y = probe_pos; z = b; }
            else { x = a; y = b; z = probe_pos; }

            double val;
            if ( axis == 0 ) val = std::abs( grid.Ey_ptr()[ grid.idx(x,y,z) ] );
            else if ( axis == 1 ) val = std::abs( grid.Ez_ptr()[ grid.idx(x,y,z) ] );
            else val = std::abs( grid.Ex_ptr()[ grid.idx(x,y,z) ] );
            if ( val > peak ) peak = val;
        }

    return 20.0 * std::log10( std::max( peak, 1e-30 ) );
}

TEST(Integration, CPMLReflectionAllAxes) {
    // The CPML has independently coded face loops for x, y, and z with
    // different ψ indexing. A sign or index error in one axis would cause
    // that face to amplify instead of absorb. Test all three.
    // Measured: all axes give -36 to -37 dB on a 60^3 isotropic grid.
    double R[3];
    for ( int axis = 0; axis < 3; ++axis ) {
        R[axis] = cpml_reflection_on_axis( axis );
        ASSERT_LT( R[axis], -30.0 );
        ASSERT_GT( R[axis], -80.0 );  // sanity: signal should exist
    }

    // All axes should agree within 5 dB (isotropic grid):
    ASSERT_LT( std::abs( R[0] - R[1] ), 5.0 );
    ASSERT_LT( std::abs( R[0] - R[2] ), 5.0 );
}

TEST(Integration, CPMLWorksWithNonUnityMaterials) {
    // The PML ψ corrections must use Cb/Db (material-dependent) not bare dt.
    // With the correct formulation, all material combinations give < -25 dB.
    auto measure = []( double eps_val, double mu_val ) -> double {
        Simulation_Config cfg{};
        cfg.Nx = 60; cfg.Ny = 60; cfg.Nz = 60;
        cfg.eps = eps_val; cfg.mu = mu_val;
        cfg.compute_derived();
        Grid grid{ cfg };
        double const c{ 1.0 / std::sqrt( cfg.mu * cfg.eps ) };
        double const eta{ std::sqrt( cfg.mu / cfg.eps ) };
        double const wavelength{ 15.0 * cfg.dx };
        double const k{ 2.0 * std::numbers::pi / wavelength };
        std::size_t const pml_t{ cfg.pml_thickness };
        std::size_t const wave_start{ pml_t + 5 };
        std::size_t const wave_end{ grid.Nx() - pml_t - 1 };
        std::size_t const taper{ std::min( std::size_t{10}, (wave_end-wave_start)/4 ) };
        std::size_t const probe_pos{ pml_t + 2 };

        for ( std::size_t z = 1; z < grid.Nz()-1; ++z )
            for ( std::size_t y = 1; y < grid.Ny()-1; ++y )
                for ( std::size_t x = 1; x < grid.Nx()-1; ++x ) {
                    if ( x < wave_start || x >= wave_end ) continue;
                    double envelope{ 1.0 };
                    if ( x < wave_start + taper ) {
                        double t = static_cast<double>(x-wave_start)/static_cast<double>(taper);
                        envelope = 0.5*(1.0-std::cos(std::numbers::pi*t));
                    } else if ( x >= wave_end - taper ) {
                        double t = static_cast<double>(wave_end-1-x)/static_cast<double>(taper);
                        envelope = 0.5*(1.0-std::cos(std::numbers::pi*t));
                    }
                    double phase = k * static_cast<double>(x) * cfg.dx;
                    std::size_t i = grid.idx(x,y,z);
                    grid.Ey_ptr()[i] = envelope * std::sin(phase);
                    grid.Hz_ptr()[i] = envelope * std::sin(phase) / eta;
                }

        double dist = static_cast<double>(wave_end-wave_start+wave_start-probe_pos)*cfg.dx;
        std::size_t steps = static_cast<std::size_t>(std::ceil(1.5*dist/(c*cfg.dt)));
        for ( std::size_t t = 0; t < steps; ++t ) grid.step();

        double peak{ 0 };
        for ( std::size_t z = grid.Nz()/3; z < 2*grid.Nz()/3; ++z )
            for ( std::size_t y = grid.Ny()/3; y < 2*grid.Ny()/3; ++y ) {
                double val = std::abs(grid.Ey_ptr()[grid.idx(probe_pos,y,z)]);
                if ( val > peak ) peak = val;
            }
        return 20.0 * std::log10( std::max(peak, 1e-30) );
    };

    ASSERT_LT( measure(4.0, 1.0), -25.0 );
    ASSERT_LT( measure(1.0, 4.0), -25.0 );
    ASSERT_LT( measure(9.0, 1.0), -25.0 );
    ASSERT_LT( measure(2.0, 2.0), -25.0 );
}