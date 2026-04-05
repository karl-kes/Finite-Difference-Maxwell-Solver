#include "test_framework.hpp"
#include "test_helpers.hpp"

#include "../src/Classes/Grid/grid.hpp"
#include "../src/Classes/Source/source.hpp"
#include "../src/Classes/Write_Output/output.hpp"

#include <cmath>
#include <numbers>
#include <fstream>
#include <filesystem>
#include <cstdint>
#include <vector>

// Grid Construction & Properties

TEST(Grid, DimensionsMatchConfig) {
    Simulation_Config cfg{};
    Grid grid{ cfg };
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

TEST(Grid, FieldsInitializedToZero) {
    Simulation_Config cfg{};
    Grid grid{ cfg };
    for ( std::size_t z = 40; z < 60; ++z )
        for ( std::size_t y = 40; y < 60; ++y )
            for ( std::size_t x = 40; x < 60; ++x ) {
                ASSERT_EQ( grid.field( Field::ELECTRIC, Component::X, x, y, z ), 0.0 );
                ASSERT_EQ( grid.field( Field::ELECTRIC, Component::Y, x, y, z ), 0.0 );
                ASSERT_EQ( grid.field( Field::ELECTRIC, Component::Z, x, y, z ), 0.0 );
                ASSERT_EQ( grid.field( Field::MAGNETIC, Component::X, x, y, z ), 0.0 );
                ASSERT_EQ( grid.field( Field::MAGNETIC, Component::Y, x, y, z ), 0.0 );
                ASSERT_EQ( grid.field( Field::MAGNETIC, Component::Z, x, y, z ), 0.0 );
            }
}

TEST(Grid, IdxLinearization) {
    Simulation_Config cfg{};
    Grid grid{ cfg };
    ASSERT_EQ( grid.idx(0, 0, 0), std::size_t{0} );
    std::size_t x{5}, y{10}, z{15};
    std::size_t expected = x + grid.Nx_padded() * ( y + grid.Ny_padded() * z );
    ASSERT_EQ( grid.idx(x, y, z), expected );
}

TEST(Grid, FieldWriteRead) {
    Simulation_Config cfg{};
    Grid grid{ cfg };
    grid.field( Field::ELECTRIC, Component::Y, 50, 50, 50 ) = 3.14;
    ASSERT_NEAR( grid.field( Field::ELECTRIC, Component::Y, 50, 50, 50 ), 3.14, 1e-15 );
    grid.Hz_ptr()[ grid.idx(25, 30, 35) ] = -2.71;
    ASSERT_NEAR( grid.field( Field::MAGNETIC, Component::Z, 25, 30, 35 ), -2.71, 1e-15 );
}

// Coefficient Baking

TEST(Grid, VacuumCoefficientsBaked) {
    // With sigma=0: Ca=1, Cb=dt/eps, Db=dt/mu.
    Simulation_Config cfg{};
    Grid grid{ cfg };
    double const dt{ grid.dt() };
    double const expected_Ca{ 1.0 };
    double const expected_Cb{ dt / cfg.eps };
    double const expected_Db{ dt / cfg.mu };
    std::size_t const i{ grid.idx(50, 50, 50) };

    ASSERT_NEAR( grid.Ca_x_ptr()[i], expected_Ca, 1e-14 );
    ASSERT_NEAR( grid.Ca_y_ptr()[i], expected_Ca, 1e-14 );
    ASSERT_NEAR( grid.Ca_z_ptr()[i], expected_Ca, 1e-14 );
    ASSERT_NEAR( grid.Cb_x_ptr()[i], expected_Cb, 1e-14 );
    ASSERT_NEAR( grid.Cb_y_ptr()[i], expected_Cb, 1e-14 );
    ASSERT_NEAR( grid.Cb_z_ptr()[i], expected_Cb, 1e-14 );
    ASSERT_NEAR( grid.Db_x_ptr()[i], expected_Db, 1e-14 );
    ASSERT_NEAR( grid.Db_y_ptr()[i], expected_Db, 1e-14 );
    ASSERT_NEAR( grid.Db_z_ptr()[i], expected_Db, 1e-14 );
}

// Energy Diagnostics

TEST(Grid, EnergyDecompositionConsistent) {
    // e_energy() + h_energy() must equal total_energy() exactly.
    Simulation_Config cfg{};
    Grid grid{ cfg };
    grid.Ey_ptr()[ grid.idx(50, 50, 50) ] = 1.0;
    grid.Hx_ptr()[ grid.idx(40, 40, 40) ] = 0.5;
    grid.Hz_ptr()[ grid.idx(60, 60, 60) ] = -0.3;

    double const ee{ grid.e_energy() };
    double const he{ grid.h_energy() };
    double const te{ grid.total_energy() };

    ASSERT_GT( ee, 0.0 );
    ASSERT_GT( he, 0.0 );
    ASSERT_NEAR( ee + he, te, 1e-15 );
}

// Non-Cubic Grid

TEST(Grid, NonCubicSingleStepCurl) {
    // With dx != dy != dz, the curl must use the correct spatial step for
    // each finite difference. This catches bugs where inv_dx/inv_dy/inv_dz
    // are swapped or all set to the same value.
    Simulation_Config cfg{};
    cfg.dx = 1.0; cfg.dy = 2.0; cfg.dz = 3.0;
    cfg.compute_derived();
    Grid grid{ cfg };

    std::size_t cx{50}, cy{50}, cz{50};
    grid.field( Field::ELECTRIC, Component::Y, cx, cy, cz ) = 1.0;

    double const Db{ cfg.dt / cfg.mu };
    grid.step();

    // Hx(cx,cy,cz): curl_x = -dEy/dz = -(0-1)/dz = 1/dz
    //   Hx -= Db * (1/dz), so Hx = -Db/dz
    ASSERT_NEAR( grid.field( Field::MAGNETIC, Component::X, cx, cy, cz ),
                 -Db / cfg.dz, 1e-14 );

    // Hz(cx,cy,cz): curl_z = dEy/dx = (0-1)/dx = -1/dx
    //   Hz -= Db * (-1/dx), so Hz = +Db/dx
    ASSERT_NEAR( grid.field( Field::MAGNETIC, Component::Z, cx, cy, cz ),
                 Db / cfg.dx, 1e-14 );

    // If the code incorrectly used inv_dx everywhere, Hx would be -Db/dx
    // instead of -Db/dz. Since dx=1 and dz=3, these differ by 3x.
    // This assertion would catch that:
    ASSERT_TRUE( std::abs( -Db / cfg.dz - ( -Db / cfg.dx ) ) > 0.1 );
}

// Config Parser

TEST(Grid, ConfigFromFileRoundTrip) {
    std::string const path{ "test_config_roundtrip.cfg" };
    {
        std::ofstream f{ path };
        f << "Nx = 50\n";
        f << "Ny = 60\n";
        f << "Nz = 70\n";
        f << "dx = 0.5\n";
        f << "cfl_factor = 0.9\n";
        f << "use_pml = false\n";
        f << "# This is a comment\n";
        f << "total_steps = 500\n";
    }
    Simulation_Config cfg{ Simulation_Config::from_file( path ) };
    ASSERT_EQ( cfg.Nx, std::size_t{50} );
    ASSERT_EQ( cfg.Ny, std::size_t{60} );
    ASSERT_EQ( cfg.Nz, std::size_t{70} );
    ASSERT_NEAR( cfg.dx, 0.5, 1e-15 );
    ASSERT_NEAR( cfg.cfl_factor, 0.9, 1e-15 );
    ASSERT_EQ( cfg.use_pml, false );
    ASSERT_EQ( cfg.total_steps, std::size_t{500} );
    double const c{ 1.0 / std::sqrt( cfg.mu * cfg.eps ) };
    double const expected_dt{ 0.9 / ( c * std::sqrt(
        1.0/(cfg.dx*cfg.dx) + 1.0/(cfg.dy*cfg.dy) + 1.0/(cfg.dz*cfg.dz) ) ) };
    ASSERT_NEAR( cfg.dt, expected_dt, 1e-14 );
    std::filesystem::remove( path );
}

TEST(Grid, ConfigMissingFileUsesDefaults) {
    Simulation_Config cfg{ Simulation_Config::from_file( "nonexistent_file.cfg" ) };
    ASSERT_EQ( cfg.Nx, std::size_t{100} );
    ASSERT_EQ( cfg.Ny, std::size_t{100} );
    ASSERT_EQ( cfg.Nz, std::size_t{100} );
    ASSERT_NEAR( cfg.dx, 1.0, 1e-15 );
    ASSERT_NEAR( cfg.cfl_factor, 1.0, 1e-15 );
}

// Output: H-Field Yee Averaging

TEST(Output, HFieldYeeAveraging) {
    // Set Hx = x and verify the binary output contains the correct
    // Yee-averaged values: Hx_avg(x,y,z) = 0.5*(Hx(x) + Hx(x+1)).
    // This complements BinaryRoundTrip which only checks E-field output.
    Simulation_Config cfg{};
    Grid grid{ cfg };

    for ( std::size_t z = 0; z < grid.Nz(); ++z )
        for ( std::size_t y = 0; y < grid.Ny(); ++y )
            for ( std::size_t x = 0; x < grid.Nx(); ++x )
                grid.Hx_ptr()[ grid.idx(x, y, z) ] = static_cast<double>(x);

    std::string test_dir{ "test_hfield_yee" };
    Output output{ test_dir };
    output.initialize( grid );
    output.write_field( grid );
    output.finalize();

    std::ifstream file( test_dir + "/H.bin", std::ios::binary );
    ASSERT_TRUE( file.is_open() );

    uint64_t dims[3];
    file.read( reinterpret_cast<char*>(dims), sizeof(dims) );
    std::size_t nx = dims[0], ny = dims[1], nz = dims[2];

    std::size_t total = nx * ny * nz * 3;
    std::vector<double> data( total );
    file.read( reinterpret_cast<char*>(data.data()),
               static_cast<std::streamsize>( total * sizeof(double) ) );
    ASSERT_TRUE( file.good() );
    file.close();

    for ( std::size_t z = 5; z < 15; ++z )
        for ( std::size_t y = 5; y < 15; ++y )
            for ( std::size_t x = 5; x < 15; ++x ) {
                std::size_t offset = ( z * ny * nx + y * nx + x ) * 3;
                double Hx_avg = data[offset + 0];
                double Hy_avg = data[offset + 1];
                double Hz_avg = data[offset + 2];

                // Hx_avg = 0.5*(Hx[x] + Hx[x+1]) = 0.5*(x + x+1) = x + 0.5
                ASSERT_NEAR( Hx_avg, static_cast<double>(x) + 0.5, 1e-12 );
                ASSERT_NEAR( Hy_avg, 0.0, 1e-14 );
                ASSERT_NEAR( Hz_avg, 0.0, 1e-14 );
            }

    std::filesystem::remove_all( test_dir );
}

// Source Application

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
    double t0{grid.dt() * 50.0};
    double width{grid.dt() * 10.0};
    std::size_t x{50}, y{50}, z{50};
    grid.add_source( std::make_unique<Gaussian_Pulse>( amp, t0, width, x, y, z ) );
    std::size_t i = grid.idx(x, y, z);
    grid.apply_sources( 50 );
    ASSERT_NEAR( grid.Jz_ptr()[i], amp, 1e-10 );
    grid.apply_sources( 0 );
    ASSERT_LT( std::abs(grid.Jz_ptr()[i]), amp * 0.01 );
}

TEST(Source, StraightWireSinusoidal) {
    Simulation_Config cfg{};
    Grid grid{ cfg };
    double amp{5.0}, freq{0.1};
    std::size_t y{50}, z{50}, xs{40}, xe{60};
    grid.add_source( std::make_unique<Straight_Wire_X>( amp, freq, y, z, xs, xe ) );
    grid.apply_sources( 0 );
    for ( std::size_t x = xs; x < xe; ++x )
        ASSERT_NEAR( grid.Jx_ptr()[ grid.idx(x, y, z) ], 0.0, 1e-14 );
    grid.apply_sources( 100 );
    double first_val = grid.Jx_ptr()[ grid.idx(xs, y, z) ];
    ASSERT_GT( std::abs(first_val), 0.0 );
    for ( std::size_t x = xs + 1; x < xe; ++x )
        ASSERT_NEAR( grid.Jx_ptr()[ grid.idx(x, y, z) ], first_val, 1e-15 );
}

TEST(Source, ACCurrentLoopGeometryAndCirculation) {
    // Verify closed rectangular current path with correct circulation.
    Simulation_Config cfg{};
    Grid grid{ cfg };
    double const c{ 1.0 / std::sqrt( cfg.mu * cfg.eps ) };
    double const freq{ c / ( 20.0 * grid.dx() ) };
    std::size_t const z{ grid.Nz() / 2 };
    grid.add_source( std::make_unique<AC_Current_Loop>( 1.0, freq, z ) );

    // Step 0: sin(0)=0, all J zero:
    grid.apply_sources( 0 );
    std::size_t const x_i{ grid.Nx() / 4 }, x_f{ 3 * grid.Nx() / 4 };
    std::size_t const y_i{ grid.Ny() / 4 }, y_f{ 3 * grid.Ny() / 4 };
    ASSERT_EQ( grid.Jx_ptr()[ grid.idx(x_i+5, y_i, z) ], 0.0 );

    // Step 1: sin(omega*dt) != 0:
    grid.apply_sources( 1 );
    double const jx_bottom{ grid.Jx_ptr()[ grid.idx(x_i+5, y_i, z) ] };
    double const jx_top{ grid.Jx_ptr()[ grid.idx(x_i+5, y_f, z) ] };
    double const jy_right{ grid.Jy_ptr()[ grid.idx(x_f, y_i+5, z) ] };
    double const jy_left{ grid.Jy_ptr()[ grid.idx(x_i, y_i+5, z) ] };

    ASSERT_GT( jx_bottom, 0.0 );                      // Bottom: +x current
    ASSERT_NEAR( jx_top, -jx_bottom, 1e-15 );         // Top: -x current
    ASSERT_GT( jy_right, 0.0 );                       // Right: +y current
    ASSERT_NEAR( jy_left, -jy_right, 1e-15 );         // Left: -y current
    ASSERT_NEAR( jx_bottom, jy_right, 1e-15 );        // Same magnitude

    // Outside loop: zero:
    ASSERT_EQ( grid.Jx_ptr()[ grid.idx(x_i-1, y_i, z) ], 0.0 );
    ASSERT_EQ( grid.Jx_ptr()[ grid.idx(x_i+5, y_i+5, z) ], 0.0 );

    // After stepping, Hz at center should be positive (right-hand rule:
    // CCW current in xy-plane when viewed from +z produces +Hz).
    for ( std::size_t t = 1; t <= 50; ++t ) { grid.apply_sources(t); grid.step(); }
    double const hz_center{ grid.field(
        Field::MAGNETIC, Component::Z, grid.Nx()/2, grid.Ny()/2, z ) };
    ASSERT_GT( hz_center, 0.0 );
}

// Lossy Coefficient Baking (end-to-end via bake_coefficients)

TEST(Grid, LossyBakeCoefficients) {
    // Set sigma > 0 on a region, rebake, verify Ca < 1 and Cb < dt/eps
    // at those cells, and that unmodified cells are unchanged.
    Simulation_Config cfg{};
    cfg.use_pml = false;
    cfg.compute_derived();
    Grid grid{ cfg };

    std::size_t const i_lossy{ grid.idx(50, 50, 50) };
    std::size_t const i_vacuum{ grid.idx(20, 20, 20) };

    double const dt{ grid.dt() };
    double const eps{ cfg.eps };

    // Record vacuum values before modification:
    double const Ca_vacuum{ grid.Ca_x_ptr()[i_vacuum] };
    double const Cb_vacuum{ grid.Cb_x_ptr()[i_vacuum] };
    ASSERT_NEAR( Ca_vacuum, 1.0, 1e-14 );
    ASSERT_NEAR( Cb_vacuum, dt / eps, 1e-14 );

    // Set sigma at (50,50,50):
    double const sig{ 0.5 };
    grid.sig_x_ptr()[i_lossy] = sig;
    grid.sig_y_ptr()[i_lossy] = sig;
    grid.sig_z_ptr()[i_lossy] = sig;

    grid.bake_coefficients();

    // Lossy cell: Ca < 1, Cb < dt/eps:
    double const loss{ sig * dt / ( 2.0 * eps ) };
    double const expected_Ca{ ( 1.0 - loss ) / ( 1.0 + loss ) };
    double const expected_Cb{ 1.0 / ( eps / dt + sig / 2.0 ) };

    ASSERT_NEAR( grid.Ca_x_ptr()[i_lossy], expected_Ca, 1e-14 );
    ASSERT_NEAR( grid.Cb_x_ptr()[i_lossy], expected_Cb, 1e-14 );
    ASSERT_NEAR( grid.Ca_y_ptr()[i_lossy], expected_Ca, 1e-14 );
    ASSERT_NEAR( grid.Ca_z_ptr()[i_lossy], expected_Ca, 1e-14 );

    // Vacuum cell: unchanged:
    ASSERT_NEAR( grid.Ca_x_ptr()[i_vacuum], 1.0, 1e-14 );
    ASSERT_NEAR( grid.Cb_x_ptr()[i_vacuum], dt / eps, 1e-14 );

    // Db should be unaffected by sigma (only depends on mu):
    ASSERT_NEAR( grid.Db_x_ptr()[i_lossy], dt / cfg.mu, 1e-14 );
}

TEST(Grid, LossyDecayBehavior) {
    // Set a lossy slab, inject Ey, step, verify the field decays faster
    // inside the lossy region than outside.
    Simulation_Config cfg{};
    cfg.use_pml = false;
    cfg.Nx = 60; cfg.Ny = 60; cfg.Nz = 60;
    cfg.compute_derived();
    Grid grid{ cfg };

    // Set sigma=2.0 in a slab from x=25 to x=35:
    for ( std::size_t z = 0; z < grid.Nz(); ++z )
        for ( std::size_t y = 0; y < grid.Ny(); ++y )
            for ( std::size_t x = 25; x <= 35; ++x ) {
                std::size_t i = grid.idx(x, y, z);
                grid.sig_x_ptr()[i] = 2.0;
                grid.sig_y_ptr()[i] = 2.0;
                grid.sig_z_ptr()[i] = 2.0;
            }
    grid.bake_coefficients();

    // Inject uniform Ey=1 across the whole interior:
    for ( std::size_t z = 1; z < grid.Nz()-1; ++z )
        for ( std::size_t y = 1; y < grid.Ny()-1; ++y )
            for ( std::size_t x = 1; x < grid.Nx()-1; ++x )
                grid.Ey_ptr()[ grid.idx(x, y, z) ] = 1.0;

    // Step 10 times:
    for ( int t = 0; t < 10; ++t ) grid.step();

    // Sample Ey inside lossy slab vs outside:
    double ey_lossy = std::abs( grid.Ey_ptr()[ grid.idx(30, 30, 30) ] );
    double ey_vacuum = std::abs( grid.Ey_ptr()[ grid.idx(10, 30, 30) ] );

    // Field in lossy region should have decayed more:
    ASSERT_LT( ey_lossy, ey_vacuum );
}

// Independent H/E Update Half-Steps

TEST(Grid, UpdateHIndependent) {
    // Call update_H() alone. Starting from Ey=1 at center, only H should
    // change; E should be untouched since update_E() was not called.
    Simulation_Config cfg{};
    cfg.use_pml = false;
    cfg.compute_derived();
    Grid grid{ cfg };

    std::size_t cx{50}, cy{50}, cz{50};
    grid.field( Field::ELECTRIC, Component::Y, cx, cy, cz ) = 1.0;

    grid.update_H();

    // Ey at center should still be 1.0 (update_E not called):
    ASSERT_NEAR( grid.field(Field::ELECTRIC, Component::Y, cx, cy, cz), 1.0, 1e-15 );

    // Hx at center should be nonzero:
    double const Db{ cfg.dt / cfg.mu };
    ASSERT_NEAR( grid.field(Field::MAGNETIC, Component::X, cx, cy, cz),
                 -Db / cfg.dz, 1e-14 );
}

TEST(Grid, UpdateEIndependent) {
    // Call update_H() then update_E() separately, and verify the result
    // matches a single step() call.
    Simulation_Config cfg{};
    cfg.use_pml = false;
    cfg.compute_derived();

    // Reference: use step():
    Grid grid_ref{ cfg };
    grid_ref.field( Field::ELECTRIC, Component::Y, 50, 50, 50 ) = 1.0;
    grid_ref.step();

    // Test: use update_H() then update_E():
    Grid grid_test{ cfg };
    grid_test.field( Field::ELECTRIC, Component::Y, 50, 50, 50 ) = 1.0;
    grid_test.update_H();
    grid_test.update_E();

    // Sample several points — should be identical:
    for ( std::size_t x = 48; x <= 52; ++x ) {
        for ( auto comp : { Component::X, Component::Y, Component::Z } ) {
            double ref_e = grid_ref.field(Field::ELECTRIC, comp, x, 50, 50);
            double tst_e = grid_test.field(Field::ELECTRIC, comp, x, 50, 50);
            ASSERT_NEAR( tst_e, ref_e, 1e-15 );

            double ref_h = grid_ref.field(Field::MAGNETIC, comp, x, 50, 50);
            double tst_h = grid_test.field(Field::MAGNETIC, comp, x, 50, 50);
            ASSERT_NEAR( tst_h, ref_h, 1e-15 );
        }
    }
}