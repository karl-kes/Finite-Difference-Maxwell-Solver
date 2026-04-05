#include "simulation.hpp"

#include "../Source/source.hpp"

#include <iostream>
#include <chrono>
#include <memory>
#include <cstddef>

Simulation::Simulation( Simulation_Config const &new_config )
: config_{ new_config }
, grid_{ config_ }
, output_{ "output" }
{ }

void Simulation::print_progress( std::size_t const current, std::size_t const total ) const {
    double const percent{ 100.0 * static_cast<double>( current ) / static_cast<double>( total ) };
    std::cout << "\rProgress: " << percent << "%" << std::flush;
}

void Simulation::initialize_sources() {
    // Constants:
    double const initial_time{ 0 };
    double const c{ 1.0 / std::sqrt( config().eps * config().mu ) };
    
    // Calculate usable grid space excluding the absorbing boundary:
    std::size_t const pml_t{ config().use_pml ? config().pml_thickness : 0 };
    std::size_t const usable_cells{ grid().Nx() > 2 * pml_t ? grid().Nx() - 2 * pml_t : 1 };

    // Enforce a strict minimum of 40 cells for high resolution:
    double const cells_per_wavelength{ std::max( 40.0, static_cast<double>( usable_cells ) * 0.2 ) };
    double const freq{ c / ( cells_per_wavelength * grid().dx() ) };
    double const amp{ 1.0 };

    // Grid points:
    std::size_t const half_z{ grid().Nz() / 2 };
    std::size_t const half_y{ grid().Ny() / 2 };
    std::size_t const half_x{ grid().Nx() / 2 };

    std::size_t const start_x{ 0 };
    std::size_t const end_x{ grid().Nx() };

    // AC current wire loop:
    grid().add_source( std::make_unique<AC_Current_Loop>(
        amp,
        freq,
        half_z
    ) );

    // // AC current concentric rings:
    // grid().add_source( std::make_unique<AC_Concentric_Rings>(
    //     amp,
    //     freq,
    //     half_z
    // ) );

    // // AC current wire along x:
    // grid().add_source( std::make_unique<Straight_Wire_X>(
    //     amp,
    //     freq,
    //     half_y,
    //     half_z,
    //     start_x,
    //     end_x
    // ) );

    // // Point source:
    // grid().add_source( std::make_unique<Point_Source>(
    //     amp,
    //     half_x,
    //     half_y,
    //     half_z
    // ) );

    // // Gaussian pulse:
    // grid().add_source( std::make_unique<Gaussian_Pulse>(
    //     amp,
    //     initial_time,
    //     freq,
    //     half_x,
    //     half_y,
    //     half_z
    // ) );
}

void Simulation::initialize() {
    output().initialize( grid() );
    initialize_sources();
}

void Simulation::run() {
    // Initialize loop variables:
    std::size_t const output_interval{ config().output_interval() };
    std::size_t const total_steps{ config().total_steps };

    std::cout << "<---- Maxwell Simulation ---->" << std::endl;
    auto const start_time{ std::chrono::high_resolution_clock::now() };

    // Simulation Loop:
    for ( std::size_t current_step{}; current_step < total_steps; ++current_step ) {
        grid().apply_sources( current_step );
        grid().step();

        if ( ( current_step % output_interval ) == 0 ) {
            output().write_field( grid() );
            print_progress( current_step, total_steps );
        }
    }

    // Ensure progress bar ends at 100%:
    std::cout << "\rProgress: " << 100.0 << "%" << std::flush;

    // End Outputs:
    output().finalize();

    // End Timer:
    auto const end_time{ std::chrono::high_resolution_clock::now() };
    auto const duration{ std::chrono::duration_cast<std::chrono::milliseconds>( end_time - start_time ) };

    // Report Results:  
    std::cout << "\n\nResults:\n";
    std::cout << "--------\n";
    std::cout << "Duration: " << duration.count() << " ms\n";
    std::cout << "Physical time: " << static_cast<double>( total_steps ) * grid().dt() << " s\n";
}