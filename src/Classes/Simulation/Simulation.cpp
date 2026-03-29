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

void Simulation::initialize() {
    output_.initialize();

    double const freq{ grid_.c() / ( 20.0 * grid_.dx() ) };
    grid_.add_source( std::make_unique<Straight_Wire_X>(
        10.0,                    // amplitude
        freq,                    // frequency
        config_.Ny / 2,          // y position
        config_.Nz / 2,          // z position
        config_.Nx / 4,          // x_start
        3 * config_.Nx / 4       // x_end
    ) );
}

void Simulation::run() {
    std::cout << "<-----Maxwell Simulation----->" << std::endl;

    // Run @ t = 0:
    grid_.apply_sources();
    grid_.step();

    // Run simulation and start timer:
    std::size_t const output_interval{ config_.output_interval() };
    auto const start_time{ std::chrono::high_resolution_clock::now() };

    // Simulation Loop:
    for ( std::size_t curr_time{1}; curr_time <= config_.total_steps; ++curr_time ) {
        grid_.apply_sources( curr_time );
        grid_.step();

        if ( ( curr_time % output_interval ) == 0 ) {
            output_.write_field( grid_, curr_time );
            print_progress( curr_time, config_.total_steps );
        }
    }

    // End Timer:
    auto const end_time{ std::chrono::high_resolution_clock::now() };
    auto const duration{ std::chrono::duration_cast<std::chrono::milliseconds>( end_time - start_time ) };

    // Report Results:  
    std::cout << "\n\nResults:\n";
    std::cout << "--------\n";
    std::cout << "Duration: " << duration.count() << " ms\n";
    std::cout << "Physical time: " << static_cast<double>( config_.total_steps ) * grid_.dt() << " s\n";
}