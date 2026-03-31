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
    // Write a known field state to binary, read it back using the same format
    // as render.py, and verify the Yee-averaged values match.
    Simulation_Config cfg{};
    Grid grid{ cfg };

    for ( std::size_t z = 0; z < grid.Nz(); ++z ) {
        for ( std::size_t y = 0; y < grid.Ny(); ++y ) {
            for ( std::size_t x = 0; x < grid.Nx(); ++x ) {
                grid.Ey_ptr()[ grid.idx(x, y, z) ] = static_cast<double>(x);
            }
        }
    }

    // Write:
    std::string test_dir{ "test_output_roundtrip" };
    Output output{ test_dir };
    output.initialize( grid );
    output.write_field( grid );
    output.finalize();

    // Read back the E-field file:
    std::string path{ test_dir + "/E.bin" };
    std::ifstream file( path, std::ios::binary );
    ASSERT_TRUE( file.is_open() );

    // Read header:
    uint64_t dims[3];
    file.read( reinterpret_cast<char*>(dims), sizeof(dims) );
    std::size_t nx = static_cast<std::size_t>(dims[0]);
    std::size_t ny = static_cast<std::size_t>(dims[1]);
    std::size_t nz = static_cast<std::size_t>(dims[2]);

    ASSERT_EQ( nx, grid.Nx() - 1 );
    ASSERT_EQ( ny, grid.Ny() - 1 );
    ASSERT_EQ( nz, grid.Nz() - 1 );

    // Read all data:
    std::size_t total_floats = nx * ny * nz * 3;
    std::vector<double> data( total_floats );
    file.read( reinterpret_cast<char*>(data.data()),
               static_cast<std::streamsize>(total_floats * sizeof(double)) );
    ASSERT_TRUE( file.good() );
    file.close();

    // Verify Yee-averaged Ey values at sample points.
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

    // Also verify the H-field file was created:
    std::string path_H{ test_dir + "/H.bin" };
    ASSERT_TRUE( std::filesystem::exists( path_H ) );

    // Cleanup:
    std::filesystem::remove_all( test_dir );
}

TEST(Output, HeaderDimensionsConsistent) {
    // The binary header should match the grid dimensions for both E and H files.
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
// Validation Class Coverage
// ============================================================================

TEST(ValidationClass, PlaneWaveTestPasses) {
    Simulation_Config cfg{};
    Plane_Wave_Test test{ cfg };
    Validation_Result result = test.run();

    ASSERT_TRUE( result.passed );
    ASSERT_LT( result.energy_drift_percent, 5.0 );
    ASSERT_GT( result.phase_correlation, 0.97 );
    ASSERT_LT( result.dispersion_percent, 10.0 );
}

TEST(ValidationClass, MetricsPhysicallyReasonable) {
    Simulation_Config cfg{};
    Plane_Wave_Test test{ cfg };
    Validation_Result result = test.run();

    ASSERT_GT( result.energy_drift_percent, 0.0 );
    ASSERT_GT( result.phase_correlation, 0.999 );
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