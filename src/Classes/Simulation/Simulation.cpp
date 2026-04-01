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
    double const c{ 1.0 / std::sqrt( config_.eps * config_.mu ) };
    double const cells_per_wavelength{ static_cast<double>( grid_.Nx() ) / 2.0 };
    double const freq{ c / ( cells_per_wavelength * grid_.dx() ) };
    double const amp{ 1.0 };

    // Grid points:
    std::size_t const half_z{ grid_.Nz() / 2 };
    std::size_t const half_y{ grid_.Ny() / 2 };
    std::size_t const half_x{ grid_.Nx() / 2 };

    std::size_t const start_x{ 0 };
    std::size_t const end_x{ grid_.Nx() };

    // AC current wire loop:
    grid_.add_source( std::make_unique<AC_Current_Loop>(
        amp,
        freq,
        half_z
    ) );

    // // AC current concentric rings:
    // grid_.add_source( std::make_unique<AC_Concentric_Rings>(
    //     amp,
    //     freq,
    //     half_z
    // ) );

    // // AC current wire along x:
    // grid_.add_source( std::make_unique<Straight_Wire_X>(
    //     amp,
    //     freq,
    //     half_y,
    //     half_z,
    //     start_x,
    //     end_x
    // ) );

    // // Point source:
    // grid_.add_source( std::make_unique<Point_Source>(
    //     amp,
    //     half_x,
    //     half_y,
    //     half_z
    // ) );

    // // Gaussian pulse:
    // grid_.add_source( std::make_unique<Gaussian_Pulse>(
    //     amp,
    //     initial_time,
    //     freq,
    //     half_x,
    //     half_y,
    //     half_z
    // ) );
}

void Simulation::initialize() {
    output_.initialize( grid_ );
    initialize_sources();
}

void Simulation::run() {
    std::cout << "<---- Maxwell Simulation ---->" << std::endl;

    // Run simulation and start timer:
    std::size_t const output_interval{ config_.output_interval() };
    auto const start_time{ std::chrono::high_resolution_clock::now() };

    // Simulation Loop:
    for ( std::size_t curr_time{1}; curr_time <= config_.total_steps; ++curr_time ) {
        grid_.apply_sources( curr_time );
        grid_.step();

        if ( ( curr_time % output_interval ) == 0 ) {
            output_.write_field( grid_ );
            print_progress( curr_time, config_.total_steps );
        }
    }

    // End Outputs:
    output_.finalize();

    // End Timer:
    auto const end_time{ std::chrono::high_resolution_clock::now() };
    auto const duration{ std::chrono::duration_cast<std::chrono::milliseconds>( end_time - start_time ) };

    // Report Results:  
    std::cout << "\n\nResults:\n";
    std::cout << "--------\n";
    std::cout << "Duration: " << duration.count() << " ms\n";
    std::cout << "Physical time: " << static_cast<double>( config_.total_steps ) * grid_.dt() << " s\n";
}