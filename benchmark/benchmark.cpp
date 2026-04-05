#include "../src/Classes/Config/config.hpp"
#include "../src/Classes/Grid/grid.hpp"
#include "../src/Classes/Source/source.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <vector>
#include <string>
#include <algorithm>
#include <memory>
#include <stdexcept>

#include <omp.h>

// Helpers

static Grid make_grid( std::size_t const N, bool const use_pml ) {
    Simulation_Config cfg{};
    cfg.Nx = N;
    cfg.Ny = N;
    cfg.Nz = N;
    cfg.dx = 1.0;
    cfg.dy = 1.0;
    cfg.dz = 1.0;
    cfg.eps = 1.0;
    cfg.mu = 1.0;
    cfg.cfl_factor = 1.0;
    cfg.total_steps = 1;
    cfg.use_pml = use_pml;
    cfg.pml_order = 3;
    cfg.pml_kappa_max = 15.0;
    cfg.pml_alpha_max = 0.3;
    cfg.run_validation = false;
    cfg.compute_derived();

    Grid grid{ cfg };

    double const c{ 1.0 / std::sqrt( cfg.eps * cfg.mu ) };
    double const cells_per_wavelength{ static_cast<double>( grid.Nx() ) / 2.0 };
    double const freq{ c / ( cells_per_wavelength * grid.dx() ) };
    double const amp{ 1.0 };

    grid.add_source( std::make_unique<AC_Current_Loop>(
        amp,
        freq,
        grid.Nz() / 2
    ) );

    return grid;
}

struct BenchResult {
    std::size_t N;
    std::size_t total_cells;
    std::size_t steps;
    bool        pml;
    double      serial_ms;
    double      omp_ms;
    double      serial_mcells;
    double      omp_mcells;
    double      speedup;
    int         omp_threads;
};

static void advance_grid_with_sources( Grid &grid, std::size_t const steps ) {
    for ( std::size_t s{}; s < steps; ++s ) {
        grid.apply_sources( s );
        grid.step();
    }
}

static double run_timed(
    std::size_t const N,
    std::size_t const steps,
    bool const use_pml,
    int const num_threads
) {
    int const prev_threads{ omp_get_max_threads() };
    omp_set_num_threads( num_threads );

    Grid grid{ make_grid( N, use_pml ) };

    // Warmup
    advance_grid_with_sources( grid, 3 );

    auto const start{ std::chrono::high_resolution_clock::now() };
    advance_grid_with_sources( grid, steps );
    auto const end{ std::chrono::high_resolution_clock::now() };

    omp_set_num_threads( prev_threads );

    return std::chrono::duration<double, std::milli>( end - start ).count();
}

static double median_ms(
    std::size_t const N,
    std::size_t const steps,
    bool const use_pml,
    int const num_threads,
    std::size_t const trials
) {
    std::vector<double> samples{};
    samples.reserve( trials );

    for ( std::size_t t{}; t < trials; ++t ) {
        samples.push_back( run_timed( N, steps, use_pml, num_threads ) );
    }

    std::sort( samples.begin(), samples.end() );

    if ( trials % 2 == 1 ) {
        return samples[trials / 2];
    }

    std::size_t const hi{ trials / 2 };
    std::size_t const lo{ hi - 1 };
    return 0.5 * ( samples[lo] + samples[hi] );
}

static double mcells_per_sec(
    std::size_t const total_cells,
    std::size_t const steps,
    double const ms
) {
    return static_cast<double>( total_cells ) * static_cast<double>( steps )
         / ( ms * 1e3 );
}

int main( int argc, char* argv[] ) {
    std::size_t steps{ 100 };
    std::size_t max_n{ 400 };
    std::size_t trials{ 5 };
    bool pml_on{ true };

    int const max_threads{ omp_get_max_threads() };
    int omp_threads{ max_threads };

    for ( int i{ 1 }; i < argc; ++i ) {
        std::string const arg{ argv[i] };
        if      ( arg == "--steps"   && i + 1 < argc ) { steps = std::stoull( argv[++i] ); }
        else if ( arg == "--max-n"   && i + 1 < argc ) { max_n = std::stoull( argv[++i] ); }
        else if ( arg == "--threads" && i + 1 < argc ) { omp_threads = std::stoi( argv[++i] ); }
        else if ( arg == "--trials"  && i + 1 < argc ) { trials = std::stoull( argv[++i] ); }
        else if ( arg == "--no-pml" )                  { pml_on = false; }
        else if ( arg == "-h" || arg == "--help" ) {
            std::cout << "Usage: benchmark [--steps N] [--max-n N] [--threads N] [--trials N] [--no-pml]\n";
            std::cout << "  --threads N   OpenMP thread count for parallel run (default: max available)\n";
            std::cout << "  --trials N    Timed trials per configuration, report median (default: 3)\n";
            return 0;
        }
    }

    if ( omp_threads < 1 ) {
        throw std::invalid_argument{ "--threads must be >= 1" };
    }
    if ( trials < 1 ) {
        throw std::invalid_argument{ "--trials must be >= 1" };
    }
    if ( omp_threads > max_threads ) {
        std::cerr << "Warning: requested " << omp_threads
                  << " threads, but omp_get_max_threads() = " << max_threads
                  << ". Using " << max_threads << " instead.\n";
        omp_threads = max_threads;
    }

    std::cout << "\n=== FDTD Yee-Kernel Scaling Benchmark ===\n";
    std::cout << "  Scheme:      Yee leapfrog (update_H + update_E)\n";
    std::cout << "  Source:      AC current loop applied every step\n";
    std::cout << "  PML:         " << ( pml_on ? "CPML (Roden-Gedney)" : "OFF" ) << "\n";
    std::cout << "  Steps/trial: " << steps << "\n";
    std::cout << "  Trials/conf: " << trials << " (median)\n";
    std::cout << "  OMP threads: " << omp_threads << "\n\n";

    std::vector<std::size_t> N_values{};
    for ( std::size_t const n : { 32, 50, 64, 80, 100, 128, 150, 200, 256, 300, 400 } ) {
        if ( n <= max_n ) {
            N_values.push_back( n );
        }
    }

    std::vector<BenchResult> results{};

    std::cout << std::left
              << std::setw(6)  << "N"
              << std::setw(12) << "Cells"
              << std::setw(7)  << "Steps"
              << std::setw(14) << "Serial (ms)"
              << std::setw(14) << "OMP (ms)"
              << std::setw(10) << "Speedup"
              << std::setw(18) << "Serial Mcells/s"
              << std::setw(18) << "OMP Mcells/s"
              << "\n";
    std::cout << std::string( 99, '=' ) << "\n";

    for ( std::size_t const N : N_values ) {
        std::size_t const total_cells{ N * N * N };

        std::size_t trial_steps{ steps };
        if ( N >= 150 ) { trial_steps = std::max<std::size_t>( steps / 4, 10 ); }
        if ( N >= 256 ) { trial_steps = std::max<std::size_t>( steps / 16, 5 ); }
        if ( N >= 400 ) { trial_steps = std::max<std::size_t>( steps / 64, 3 ); }

        double const serial_ms{ median_ms( N, trial_steps, pml_on, 1, trials ) };
        double const omp_ms{ median_ms( N, trial_steps, pml_on, omp_threads, trials ) };

        double const serial_mc{ mcells_per_sec( total_cells, trial_steps, serial_ms ) };
        double const omp_mc{ mcells_per_sec( total_cells, trial_steps, omp_ms ) };
        double const speedup{ serial_ms / omp_ms };

        results.push_back( BenchResult{
            N,
            total_cells,
            trial_steps,
            pml_on,
            serial_ms,
            omp_ms,
            serial_mc,
            omp_mc,
            speedup,
            omp_threads
        } );

        std::cout << std::left << std::fixed
                  << std::setw(6)  << N
                  << std::setw(12) << total_cells
                  << std::setw(7)  << trial_steps
                  << std::setw(14) << std::setprecision( 1 ) << serial_ms
                  << std::setw(14) << std::setprecision( 1 ) << omp_ms
                  << std::setw(10) << std::setprecision( 2 ) << speedup
                  << std::setw(18) << std::setprecision( 2 ) << serial_mc
                  << std::setw(18) << std::setprecision( 2 ) << omp_mc
                  << "\n";
    }

    std::cout << std::string( 99, '=' ) << "\n\n";

    std::cout << "Notes:\n";
    std::cout << "  - 'Serial' uses 1 OpenMP thread; 'OMP' uses " << omp_threads << " threads.\n";
    std::cout << "  - Mcells/s = (N^3 * steps) / wall_time. Standard FDTD throughput metric.\n";
    std::cout << "  - Each timed result is the median of " << trials << " runs.\n";
    std::cout << "  - Each trial includes 3 warmup steps (excluded from timing).\n";
    std::cout << "  - Steps scaled down for large grids to keep total runtime manageable.\n";
    std::cout << "  - Source ON: includes AC current loop application every step.\n";
    std::cout << "  - PML " << ( pml_on ? "ON" : "OFF" ) << ": "
              << ( pml_on ? "includes CPML ψ-update overhead on all 6 faces."
                          : "pure Yee interior kernel only." )
              << "\n\n";

    std::cout << "CSV (for plotting):\n";
    std::cout << "N,cells,steps,pml,serial_ms,omp_ms,speedup,serial_mcells_s,omp_mcells_s,threads\n";
    for ( auto const &r : results ) {
        std::cout << r.N << ","
                  << r.total_cells << ","
                  << r.steps << ","
                  << ( r.pml ? 1 : 0 ) << ","
                  << std::fixed << std::setprecision( 2 )
                  << r.serial_ms << ","
                  << r.omp_ms << ","
                  << std::setprecision( 3 ) << r.speedup << ","
                  << std::setprecision( 3 ) << r.serial_mcells << ","
                  << r.omp_mcells << ","
                  << r.omp_threads << "\n";
    }

    return 0;
}