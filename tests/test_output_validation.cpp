#include "test_framework.hpp"
#include "test_helpers.hpp"
#include "../src/Classes/Config/config.hpp"
#include "../src/Classes/Grid/grid.hpp"
#include "../src/Classes/Source/source.hpp"
#include "../src/Classes/Write_Output/output.hpp"
#include "../src/Classes/Validation/Validation.hpp"
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

    // Set a known pattern: Ey = x-index at every cell, everything else zero.
    // This makes the Yee-averaged Ey at cell (x,y,z) = 0.5*(x + (x+1)) = x + 0.5
    // (since Ey is averaged along y, but we set uniform in y, so it's just Ey[i]).
    // Actually, Ey is averaged along y: 0.5*(Ey[x,y,z] + Ey[x,y+1,z]).
    // If Ey = x everywhere, then avg = x. Let's use that.
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
    output.initialize();
    output.write_field( grid, 1 );

    // Read back the E-field file:
    std::string path = output.file_name( Field::ELECTRIC, 1 );
    std::ifstream file( path, std::ios::binary );
    ASSERT_TRUE( file.is_open() );

    // Read header:
    uint64_t dims[3];
    file.read( reinterpret_cast<char*>(dims), sizeof(dims) );
    std::size_t nx = static_cast<std::size_t>(dims[0]);
    std::size_t ny = static_cast<std::size_t>(dims[1]);
    std::size_t nz = static_cast<std::size_t>(dims[2]);

    // Output dimensions should be Nx-1, Ny-1, Nz-1:
    ASSERT_EQ( nx, grid.Nx() - 1 );
    ASSERT_EQ( ny, grid.Ny() - 1 );
    ASSERT_EQ( nz, grid.Nz() - 1 );

    // Read all data:
    std::size_t total_floats = nx * ny * nz * 4;
    std::vector<double> data( total_floats );
    file.read( reinterpret_cast<char*>(data.data()),
               static_cast<std::streamsize>(total_floats * sizeof(double)) );
    ASSERT_TRUE( file.good() );
    file.close();

    // Verify Yee-averaged Ey values at a few sample points.
    // Data layout per cell: [Ex_avg, Ey_avg, Ez_avg, |E|] — 4 doubles.
    // Storage order: x fastest, then y, then z (matching the write loops).
    for ( std::size_t z = 10; z < 20; ++z ) {
        for ( std::size_t y = 10; y < 20; ++y ) {
            for ( std::size_t x = 10; x < 20; ++x ) {
                std::size_t offset = ( z * ny * nx + y * nx + x ) * 4;
                double Ex_avg = data[offset + 0];
                double Ey_avg = data[offset + 1];
                double Ez_avg = data[offset + 2];
                double E_mag  = data[offset + 3];

                // Ex and Ez are zero everywhere:
                ASSERT_NEAR( Ex_avg, 0.0, 1e-14 );
                ASSERT_NEAR( Ez_avg, 0.0, 1e-14 );

                // Ey Yee-average: 0.5*(Ey[x,y,z] + Ey[x,y+1,z]) = 0.5*(x + x) = x
                double expected_Ey = static_cast<double>(x);
                ASSERT_NEAR( Ey_avg, expected_Ey, 1e-12 );

                // Magnitude should equal |Ey_avg| since Ex=Ez=0:
                ASSERT_NEAR( E_mag, std::abs(expected_Ey), 1e-12 );
            }
        }
    }

    // Also verify the B-field file was created:
    std::string path_B = output.file_name( Field::MAGNETIC, 1 );
    ASSERT_TRUE( std::filesystem::exists( path_B ) );

    // Cleanup:
    std::filesystem::remove_all( test_dir );
}

TEST(Output, HeaderDimensionsConsistent) {
    // The binary header should match the grid dimensions for both E and B files.
    Simulation_Config cfg{};
    Grid grid{ cfg };

    std::string test_dir{ "test_output_header" };
    Output output{ test_dir };
    output.initialize();

    // Inject a small nonzero field so the file isn't trivially zeros:
    grid.Ey_ptr()[ grid.idx(50, 50, 50) ] = 1.0;
    output.write_field( grid, 42 );

    for ( Field f : { Field::ELECTRIC, Field::MAGNETIC } ) {
        std::string path = output.file_name( f, 42 );
        std::ifstream file( path, std::ios::binary );
        ASSERT_TRUE( file.is_open() );

        uint64_t dims[3];
        file.read( reinterpret_cast<char*>(dims), sizeof(dims) );

        ASSERT_EQ( static_cast<std::size_t>(dims[0]), grid.Nx() - 1 );
        ASSERT_EQ( static_cast<std::size_t>(dims[1]), grid.Ny() - 1 );
        ASSERT_EQ( static_cast<std::size_t>(dims[2]), grid.Nz() - 1 );

        // File size should be header (24 bytes) + nx*ny*nz*4*8 bytes:
        file.seekg( 0, std::ios::end );
        auto file_size = file.tellg();
        auto expected_size = static_cast<std::streampos>(
            24 + dims[0] * dims[1] * dims[2] * 4 * sizeof(double) );
        ASSERT_EQ( file_size, expected_size );

        file.close();
    }

    std::filesystem::remove_all( test_dir );
}

// ============================================================================
// Validation Class Coverage
// ============================================================================

TEST(ValidationClass, PlaneWaveTestPasses) {
    // Exercise the actual Plane_Wave_Test class and verify it reports passing.
    Simulation_Config cfg{};
    Plane_Wave_Test test{ cfg };
    Validation_Result result = test.run( 100 );

    ASSERT_TRUE( result.passed );
    ASSERT_LT( result.energy_drift_percent, 5.0 );
    ASSERT_GT( result.phase_correlation, 0.99 );
    ASSERT_LT( result.dispersion_percent, 10.0 );
}

TEST(ValidationClass, MetricsPhysicallyReasonable) {
    // The validation metrics should be in physically sensible ranges.
    Simulation_Config cfg{};
    Plane_Wave_Test test{ cfg };
    Validation_Result result = test.run( 100 );

    // Energy drift should be positive (it's defined as |final-initial|/initial):
    ASSERT_GT( result.energy_drift_percent, 0.0 );

    // Phase correlation should be very close to 1 for a well-resolved wave:
    ASSERT_GT( result.phase_correlation, 0.999 );

    // Dispersion should be small but nonzero (numerical dispersion always exists):
    ASSERT_GT( result.dispersion_percent, 0.0 );
    ASSERT_LT( result.dispersion_percent, 2.0 );
}

// ============================================================================
// Quantitative CPML Reflection Coefficient
// ============================================================================

TEST(Integration, CPMLReflectionCoefficient) {
    // Measure CPML reflection in dB by launching a narrow-band plane wave
    // toward the high-x PML boundary and measuring the reflected amplitude
    // at a probe behind the wavefront.
    //
    // Setup: Initialize a +x propagating plane wave (Ey, Bz) in the right
    // half of the grid. After enough steps for the wavefront to enter and
    // be absorbed by PML, measure the peak |Ey| at a probe in the left half
    // (where no initial field was set — any signal there is reflection).

    Simulation_Config cfg{};
    Grid grid{ cfg };

    double const wavelength = 20.0 * grid.dx();
    double const k = 2.0 * std::numbers::pi / wavelength;

    // Initialize plane wave in x ∈ [40, 85] (right half, away from low-x PML,
    // heading into high-x PML at x ≈ 92):
    double incident_amplitude = 1.0;
    std::size_t const margin = 15;

    for ( std::size_t z = margin; z < grid.Nz() - margin; ++z ) {
        for ( std::size_t y = margin; y < grid.Ny() - margin; ++y ) {
            for ( std::size_t x = 40; x < grid.Nx() - margin; ++x ) {
                double phase = k * static_cast<double>(x) * grid.dx();
                std::size_t i = grid.idx(x, y, z);
                grid.Ey_ptr()[i] = incident_amplitude * std::sin( phase );
                grid.Bz_ptr()[i] = incident_amplitude * std::sin( phase ) / grid.c();
            }
        }
    }

    // Run long enough for the wave to travel from x=85 into the PML and back.
    // Distance: ~50 cells round trip. Steps: 50/(c*dt) ≈ 139. Use 200 for margin.
    for ( int t = 0; t < 200; ++t ) {
        grid.step();
    }

    // Measure peak reflected |Ey| at x=25 (well behind the initial wave region).
    // Any signal here must be reflected from the high-x PML.
    double peak_reflected{0};
    for ( std::size_t z = 40; z < 60; ++z ) {
        for ( std::size_t y = 40; y < 60; ++y ) {
            double val = std::abs( grid.Ey_ptr()[ grid.idx(25, y, z) ] );
            if ( val > peak_reflected ) peak_reflected = val;
        }
    }

    // Reflection coefficient in dB:
    double R_dB = 20.0 * std::log10( std::max(peak_reflected, 1e-30) / incident_amplitude );

    // With 8 PML layers, polynomial order 3, well-tuned sigma and alpha,
    // we expect reflection below -20 dB (1% reflected amplitude).
    // This is a realistic target for 8-layer CPML.
    ASSERT_LT( R_dB, -20.0 );

    // Also verify it's not unrealistically low (sanity check — there should be
    // *some* measurable reflection from a finite-thickness PML):
    ASSERT_GT( peak_reflected, 1e-20 );
}
