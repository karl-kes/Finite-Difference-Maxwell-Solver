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
    Validation_Result result = test.run();

    ASSERT_TRUE( result.passed );
    ASSERT_LT( result.energy_drift_percent, 5.0 );
    ASSERT_GT( result.phase_correlation, 0.99 );
    ASSERT_LT( result.dispersion_percent, 10.0 );
}

TEST(ValidationClass, MetricsPhysicallyReasonable) {
    // The validation metrics should be in physically sensible ranges.
    Simulation_Config cfg{};
    Plane_Wave_Test test{ cfg };
    Validation_Result result = test.run();

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
    // Key design decisions (all derived from config):
    //   1. Wave is tapered with a raised-cosine window to eliminate truncation
    //      artifacts that would produce spurious leftward signals.
    //   2. y/z margins keep the wave well inside the transverse PML.
    //   3. The probe waits long enough for the PML reflection to arrive but
    //      not so long that multi-bounce artifacts accumulate.
    //   4. The dB threshold scales with PML thickness.

    Simulation_Config cfg{};
    Grid grid{ cfg };

    std::size_t const pml_t{ cfg.pml_thickness };
    std::size_t const nx{ grid.Nx() };
    std::size_t const ny{ grid.Ny() };
    std::size_t const nz{ grid.Nz() };

    double const wavelength = 20.0 * grid.dx();
    double const k = 2.0 * std::numbers::pi / wavelength;
    double const incident_amplitude = 1.0;

    // y/z margin: well inside the transverse PML faces:
    std::size_t const yz_margin{ pml_t + 10 };

    // x layout:
    //   PML inner face (high side): nx - pml_t
    //   Wave region: [wave_start, wave_end) with taper zones on both edges
    //   Probe: in the left half, well separated from wave_start
    std::size_t const pml_inner_hi{ nx - pml_t };
    std::size_t const pml_inner_lo{ pml_t };

    // Place the wave in the right ~60% of the interior. Leave a gap between
    // the probe and the wave so even the tapered edge has negligible amplitude:
    std::size_t const interior_lo{ pml_inner_lo + 2 };
    std::size_t const interior_hi{ pml_inner_hi - 1 };
    std::size_t const interior_width{ interior_hi - interior_lo };

    std::size_t const wave_start{ interior_lo + interior_width * 2 / 5 };
    std::size_t const wave_end{ interior_hi };

    // Taper width: 2 wavelengths (in cells) on each edge, clamped to fit:
    std::size_t const taper_cells{ std::min(
        static_cast<std::size_t>( 2.0 * wavelength / grid.dx() ),
        ( wave_end - wave_start ) / 4 ) };

    // Probe x: midway between low-x PML inner face and wave_start:
    std::size_t const probe_x{ pml_inner_lo + ( wave_start - pml_inner_lo ) / 2 };

    // y/z probe region — central strip:
    std::size_t const probe_y_lo{ ny * 2 / 5 };
    std::size_t const probe_y_hi{ ny * 3 / 5 };
    std::size_t const probe_z_lo{ nz * 2 / 5 };
    std::size_t const probe_z_hi{ nz * 3 / 5 };

    // Initialize +x plane wave with raised-cosine taper on both x edges:
    for ( std::size_t z = yz_margin; z < nz - yz_margin; ++z ) {
        for ( std::size_t y = yz_margin; y < ny - yz_margin; ++y ) {
            for ( std::size_t x = wave_start; x < wave_end; ++x ) {
                // Taper envelope: 1.0 in the bulk, raised-cosine at edges:
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
                grid.Bz_ptr()[i] = incident_amplitude * envelope * std::sin( phase ) / grid.c();
            }
        }
    }

    // Step count: wave must reach high-x PML, be absorbed/reflected, and the
    // reflection must travel back to the probe. Then add margin for the
    // reflection to fully arrive:
    double const dist_to_pml{ static_cast<double>( wave_end - wave_start ) * grid.dx() };
    double const dist_return{ static_cast<double>( wave_start - probe_x ) * grid.dx() };
    double const travel_time{ ( dist_to_pml + dist_return ) / grid.c() };
    std::size_t const num_steps{ static_cast<std::size_t>(
        std::ceil( 1.5 * travel_time / grid.dt() ) ) };

    for ( std::size_t t = 0; t < num_steps; ++t ) {
        grid.step();
    }

    // Measure peak reflected |Ey| at the probe:
    double peak_reflected{ 0 };
    for ( std::size_t z = probe_z_lo; z < probe_z_hi; ++z ) {
        for ( std::size_t y = probe_y_lo; y < probe_y_hi; ++y ) {
            double val = std::abs( grid.Ey_ptr()[ grid.idx( probe_x, y, z ) ] );
            if ( val > peak_reflected ) peak_reflected = val;
        }
    }

    // Reflection coefficient in dB:
    double R_dB = 20.0 * std::log10( std::max( peak_reflected, 1e-30 ) / incident_amplitude );

    // Threshold scales with PML thickness. Empirical baseline:
    //   8 layers  → ~-15 dB      (conservative for high CFL)
    //   12 layers → ~-18 dB
    //   16 layers → ~-22 dB
    // Use: threshold = -0.8 * thickness - 8, clamped to at most -12 dB:
    double const threshold_dB{ std::min( -12.0,
        -0.8 * static_cast<double>( pml_t ) - 8.0 ) };
    ASSERT_LT( R_dB, threshold_dB );

    // Sanity: there should be some measurable reflection:
    ASSERT_GT( peak_reflected, 1e-20 );
}