#include "../src/classes/config/config.hpp"
#include "../src/classes/grid/grid.hpp"
#include "../src/classes/source/source.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <chrono>
#include <vector>
#include <string>
#include <algorithm>
#include <memory>
#include <stdexcept>
#include <thread>
#include <cstdlib>

#ifdef _OPENMP
#include <omp.h>
#endif

// -----------------------------------------------------------------------------
// Per-cell cost model (doubles, 8 bytes each). Derived from the kernels in
// Grid::update_H and Grid::update_E in the Maxwell source tree.
//
// Amortized bytes assume perfect neighbor reuse within a thread's tile
// (each field value loaded once and reused by adjacent stencil points).
// Naive bytes count every stencil access as a distinct memory transaction,
// giving an upper bound on memory traffic. Achieved bandwidth lies between
// the two, so we report both.
//
// Flops exclude the three divisions used to form inv_dx, inv_dy, inv_dz
// (hoisted out of the loop) and the neighbor-index integer math.
// -----------------------------------------------------------------------------
struct CostModel {
    // Per cell, per full step (update_H + update_E).
    static constexpr double bytes_amortized{ 240.0 }; // 30 doubles * 8
    static constexpr double bytes_naive{ 336.0 };     // 42 doubles * 8
    static constexpr double flops{ 48.0 };            // 21 (H) + 27 (E)
};

// -----------------------------------------------------------------------------
// Grid construction. Kept as a kernel microbenchmark: cfl_factor = 1.0 and
// unit material parameters mean this is a throughput test, not a validated
// physics run. Correctness is covered by the test suite.
// -----------------------------------------------------------------------------
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
    double      serial_gbs_amort;
    double      omp_gbs_amort;
    double      serial_gbs_naive;
    double      omp_gbs_naive;
    double      serial_gflops;
    double      omp_gflops;
    int         omp_threads;
};

static void advance_grid_with_sources( Grid &grid, std::size_t const steps ) {
    for ( std::size_t s{}; s < steps; ++s ) {
        grid.apply_sources( s );
        grid.step();
    }
}

// time_one_trial: runs `steps` timed kernel iterations on an already-warm
// grid. Caller owns the grid; the grid is reused across trials to keep the
// page cache and translation cache warm, which drops variance significantly
// at large N.
static double time_one_trial(
    Grid &grid,
    std::size_t const steps,
    std::size_t const source_step_offset
) {
    auto const start{ std::chrono::high_resolution_clock::now() };
    for ( std::size_t s{}; s < steps; ++s ) {
        grid.apply_sources( source_step_offset + s );
        grid.step();
    }
    auto const end{ std::chrono::high_resolution_clock::now() };
    return std::chrono::duration<double, std::milli>( end - start ).count();
}

static double median_ms(
    std::size_t const N,
    std::size_t const steps,
    bool const use_pml,
    int const num_threads,
    std::size_t const trials,
    std::size_t const warmup_steps
) {
#ifdef _OPENMP
    int const prev_threads{ omp_get_max_threads() };
    omp_set_num_threads( num_threads );
#else
    (void)num_threads;
#endif

    Grid grid{ make_grid( N, use_pml ) };

    // Warmup: first-touch page faults on every array, primes caches,
    // settles OMP runtime state, and brings frequency scaling to steady state.
    advance_grid_with_sources( grid, warmup_steps );

    std::vector<double> samples{};
    samples.reserve( trials );

    // Reuse the grid across trials: caches stay warm, so we measure steady-state
    // kernel throughput rather than cold-start cost.
    std::size_t source_offset{ warmup_steps };
    for ( std::size_t t{}; t < trials; ++t ) {
        samples.push_back( time_one_trial( grid, steps, source_offset ) );
        source_offset += steps;
    }

#ifdef _OPENMP
    omp_set_num_threads( prev_threads );
#endif

    std::sort( samples.begin(), samples.end() );

    if ( trials % 2 == 1 ) {
        return samples[trials / 2];
    }

    std::size_t const hi{ trials / 2 };
    std::size_t const lo{ hi - 1 };
    return 0.5 * ( samples[lo] + samples[hi] );
}

// When compiled WITH OpenMP, the single-thread baseline is OMP pinned to 1
// thread. This still pays parallel-region fork/barrier overhead (typically a
// few percent). For a true serial compile, build this TU with -fno-openmp; the
// _OPENMP macro then disappears, the #pragma omp directives are ignored, and
// median_ms runs the kernels straight. The `kestrel-maxwell-benchmark-serial`
// below exposes that.
static double serial_baseline_ms(
    std::size_t const N,
    std::size_t const steps,
    bool const use_pml,
    std::size_t const trials,
    std::size_t const warmup_steps
) {
    return median_ms( N, steps, use_pml, 1, trials, warmup_steps );
}

static double mcells_per_sec(
    std::size_t const total_cells,
    std::size_t const steps,
    double const ms
) {
    return static_cast<double>( total_cells ) * static_cast<double>( steps )
         / ( ms * 1e3 );
}

static double gbytes_per_sec(
    std::size_t const total_cells,
    std::size_t const steps,
    double const ms,
    double const bytes_per_cell
) {
    double const total_bytes{
        static_cast<double>( total_cells ) * static_cast<double>( steps ) * bytes_per_cell
    };
    return total_bytes / ( ms * 1e6 ); // ms -> s, bytes -> GB (1e9)
}

static double gflops_per_sec(
    std::size_t const total_cells,
    std::size_t const steps,
    double const ms
) {
    double const total_flops{
        static_cast<double>( total_cells ) * static_cast<double>( steps ) * CostModel::flops
    };
    return total_flops / ( ms * 1e6 );
}

// Machine / build identification. Reviewers expect this in any serious
// benchmark output so results are reproducible and comparable across hosts.
static std::string compiler_id() {
    return std::string{ "gcc " } + __VERSION__;
}

static std::string build_type() {
#ifdef NDEBUG
    return "Release";
#else
    return "Debug";
#endif
}

static void print_machine_info( std::ostream &os, int const omp_threads ) {
    os << "  Host hw threads: " << std::thread::hardware_concurrency() << "\n";
#ifdef _OPENMP
    os << "  OpenMP:          yes (_OPENMP=" << _OPENMP << ")\n";
    os << "  OMP max threads: " << omp_get_max_threads() << "\n";
#else
    os << "  OpenMP:          no (serial build)\n";
    (void)omp_threads;
#endif
    char const* places{ std::getenv( "OMP_PLACES" ) };
    char const* bind{ std::getenv( "OMP_PROC_BIND" ) };
    os << "  OMP_PLACES:      " << ( places ? places : "<unset>" ) << "\n";
    os << "  OMP_PROC_BIND:   " << ( bind ? bind : "<unset>" ) << "\n";
    os << "  Compiler:        " << compiler_id() << "\n";
    os << "  Build type:      " << build_type() << "\n";
    os << "  sizeof(double):  " << sizeof( double ) << "\n";
}

static void warn_if_unpinned() {
#ifdef _OPENMP
    char const* places{ std::getenv( "OMP_PLACES" ) };
    char const* bind{ std::getenv( "OMP_PROC_BIND" ) };
    if ( !places || !bind ) {
        std::cerr << "\n[warning] OMP_PLACES and/or OMP_PROC_BIND not set. "
                     "Thread migration will add variance to OMP timings. "
                     "Recommended: OMP_PLACES=cores OMP_PROC_BIND=close.\n\n";
    }
#endif
}

int main( int argc, char* argv[] ) {
    std::size_t steps{ 100 };
    std::size_t max_n{ 400 };
    std::size_t trials{ 5 };
    std::size_t warmup_steps{ 20 };
    bool pml_on{ true };
    std::string csv_path{};

#ifdef _OPENMP
    int const max_threads{ omp_get_max_threads() };
#else
    int const max_threads{ 1 };
#endif
    int omp_threads{ max_threads };

    for ( int i{ 1 }; i < argc; ++i ) {
        std::string const arg{ argv[i] };
        if      ( arg == "--steps"   && i + 1 < argc ) { steps = std::stoull( argv[++i] ); }
        else if ( arg == "--max-n"   && i + 1 < argc ) { max_n = std::stoull( argv[++i] ); }
        else if ( arg == "--threads" && i + 1 < argc ) { omp_threads = std::stoi( argv[++i] ); }
        else if ( arg == "--trials"  && i + 1 < argc ) { trials = std::stoull( argv[++i] ); }
        else if ( arg == "--warmup"  && i + 1 < argc ) { warmup_steps = std::stoull( argv[++i] ); }
        else if ( arg == "--output"  && i + 1 < argc ) { csv_path = argv[++i]; }
        else if ( arg == "--no-pml" )                  { pml_on = false; }
        else if ( arg == "-h" || arg == "--help" ) {
            std::cout << "Usage: benchmark [options]\n"
                      << "  --steps N      Timed steps per trial at base N (scaled at large N). Default 100.\n"
                      << "  --max-n N      Largest grid edge to sweep. Default 400.\n"
                      << "  --threads N    OpenMP threads for the parallel run. Default: max available.\n"
                      << "  --trials N     Timed trials per configuration; median reported. Default 5.\n"
                      << "  --warmup N     Warmup steps before timing (first-touch + cache prime). Default 20.\n"
                      << "  --no-pml       Disable CPML (pure interior kernel).\n"
                      << "  --output PATH  Write CSV to PATH (in addition to stdout).\n";
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
                  << " threads, but max available = " << max_threads
                  << ". Using " << max_threads << " instead.\n";
        omp_threads = max_threads;
    }

    warn_if_unpinned();

    std::cout << "\n=== FDTD Yee-Kernel Scaling Benchmark ===\n";
    std::cout << "  Scheme:      Yee leapfrog (update_H + update_E)\n";
    std::cout << "  Source:      AC current loop applied every step\n";
    std::cout << "  PML:         " << ( pml_on ? "CPML (Roden-Gedney)" : "OFF" ) << "\n";
    std::cout << "  Steps/trial: " << steps << " (scaled at large N)\n";
    std::cout << "  Warmup:      " << warmup_steps << " steps (excluded from timing)\n";
    std::cout << "  Trials/conf: " << trials << " (median reported)\n";
    std::cout << "  OMP threads: " << omp_threads << "\n";
    std::cout << "  Cost model:  " << CostModel::bytes_amortized << " B/cell (amort), "
              << CostModel::bytes_naive << " B/cell (naive), "
              << CostModel::flops << " flops/cell per full step\n";
    print_machine_info( std::cout, omp_threads );
    std::cout << "\n";

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
              << std::setw(13) << "1thr (ms)"
              << std::setw(13) << "OMP (ms)"
              << std::setw(9)  << "Speedup"
              << std::setw(12) << "1thr Mc/s"
              << std::setw(12) << "OMP Mc/s"
              << std::setw(12) << "OMP GB/s*"
              << std::setw(12) << "OMP GF/s"
              << "\n";
    std::cout << std::string( 108, '=' ) << "\n";

    for ( std::size_t const N : N_values ) {
        std::size_t const total_cells{ N * N * N };

        // Step scaling: large grids get fewer steps to bound wall time, with
        // a raised floor so the timed region dominates chrono jitter and
        // parallel-region entry overhead.
        std::size_t trial_steps{ steps };
        if ( N >= 150 ) { trial_steps = std::max<std::size_t>( steps / 4,  20 ); }
        if ( N >= 256 ) { trial_steps = std::max<std::size_t>( steps / 8,  15 ); }
        if ( N >= 400 ) { trial_steps = std::max<std::size_t>( steps / 16, 10 ); }

        double const serial_ms{
            serial_baseline_ms( N, trial_steps, pml_on, trials, warmup_steps )
        };
        double const omp_ms{
            median_ms( N, trial_steps, pml_on, omp_threads, trials, warmup_steps )
        };

        double const serial_mc{ mcells_per_sec( total_cells, trial_steps, serial_ms ) };
        double const omp_mc{ mcells_per_sec( total_cells, trial_steps, omp_ms ) };
        double const speedup{ serial_ms / omp_ms };

        double const serial_gbs_amort{ gbytes_per_sec( total_cells, trial_steps, serial_ms, CostModel::bytes_amortized ) };
        double const omp_gbs_amort{    gbytes_per_sec( total_cells, trial_steps, omp_ms,    CostModel::bytes_amortized ) };
        double const serial_gbs_naive{ gbytes_per_sec( total_cells, trial_steps, serial_ms, CostModel::bytes_naive ) };
        double const omp_gbs_naive{    gbytes_per_sec( total_cells, trial_steps, omp_ms,    CostModel::bytes_naive ) };

        double const serial_gflops{ gflops_per_sec( total_cells, trial_steps, serial_ms ) };
        double const omp_gflops{    gflops_per_sec( total_cells, trial_steps, omp_ms ) };

        results.push_back( BenchResult{
            N, total_cells, trial_steps, pml_on,
            serial_ms, omp_ms,
            serial_mc, omp_mc, speedup,
            serial_gbs_amort, omp_gbs_amort,
            serial_gbs_naive, omp_gbs_naive,
            serial_gflops, omp_gflops,
            omp_threads
        } );

        std::cout << std::left << std::fixed
                  << std::setw(6)  << N
                  << std::setw(12) << total_cells
                  << std::setw(7)  << trial_steps
                  << std::setw(13) << std::setprecision( 1 ) << serial_ms
                  << std::setw(13) << std::setprecision( 1 ) << omp_ms
                  << std::setw(9)  << std::setprecision( 2 ) << speedup
                  << std::setw(12) << std::setprecision( 2 ) << serial_mc
                  << std::setw(12) << std::setprecision( 2 ) << omp_mc
                  << std::setw(12) << std::setprecision( 2 ) << omp_gbs_amort
                  << std::setw(12) << std::setprecision( 2 ) << omp_gflops
                  << "\n";
    }

    std::cout << std::string( 108, '=' ) << "\n\n";

    std::cout << "Notes:\n";
    std::cout << "  - '1thr' = single-thread run. With _OPENMP defined, this is OMP pinned to 1 thread;\n";
    std::cout << "    build `kestrel-maxwell-benchmark-serial` (no -fopenmp) for a true serial baseline.\n";
    std::cout << "  - Mcells/s = (N^3 * steps) / wall_time.\n";
    std::cout << "  - GB/s* uses the AMORTIZED cost model (" << CostModel::bytes_amortized
              << " B/cell). Naive-model GB/s also in the CSV.\n";
    std::cout << "  - GF/s uses " << CostModel::flops
              << " flops/cell/step (curls + field updates; divides hoisted out).\n";
    std::cout << "  - Each timed result is the median of " << trials
              << " back-to-back trials on a warm grid (caches stay hot).\n";
    std::cout << "  - Warmup: " << warmup_steps << " steps, excluded from timing.\n";
    std::cout << "  - PML " << ( pml_on ? "ON" : "OFF" ) << ": "
              << ( pml_on ? "includes CPML psi-update overhead on all 6 faces."
                          : "pure Yee interior kernel only." )
              << "\n";
    std::cout << "  - Kernel microbenchmark: unit dx/eps/mu, cfl_factor=1.0. Not a validated physics run.\n\n";

    // --- CSV output ---------------------------------------------------------
    auto write_csv = [&]( std::ostream &os ) {
        os << "N,cells,steps,pml,threads,"
           << "serial_ms,omp_ms,speedup,"
           << "serial_mcells_s,omp_mcells_s,"
           << "serial_gbs_amort,omp_gbs_amort,"
           << "serial_gbs_naive,omp_gbs_naive,"
           << "serial_gflops,omp_gflops\n";
        for ( auto const &r : results ) {
            os << r.N << ","
               << r.total_cells << ","
               << r.steps << ","
               << ( r.pml ? 1 : 0 ) << ","
               << r.omp_threads << ","
               << std::fixed << std::setprecision( 3 )
               << r.serial_ms << ","
               << r.omp_ms << ","
               << r.speedup << ","
               << r.serial_mcells << ","
               << r.omp_mcells << ","
               << r.serial_gbs_amort << ","
               << r.omp_gbs_amort << ","
               << r.serial_gbs_naive << ","
               << r.omp_gbs_naive << ","
               << r.serial_gflops << ","
               << r.omp_gflops << "\n";
        }
    };

    std::cout << "CSV (for plotting):\n";
    write_csv( std::cout );

    if ( !csv_path.empty() ) {
        std::ofstream csv{ csv_path };
        if ( !csv ) {
            std::cerr << "Error: could not open " << csv_path << " for writing.\n";
            return 1;
        }
        write_csv( csv );
        std::cout << "\nCSV also written to: " << csv_path << "\n";
    }

    return 0;
}
