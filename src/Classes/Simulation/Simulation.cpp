#include "Simulation.hpp"

Simulation::Simulation( Simulation_Config const &new_config )
: config_{ new_config }
, grid_{ config_ }
, output_{ "output" }
{ }

void Simulation::initialize() {
    output_.initialize();

    // Add sources:
    // grid_.add_source( std::make_unique<Straight_Wire_X>(
    //     10.0,                           // amplitude
    //     1.0,                            // frequency
    //     config_.Ny / 2,                  // y position
    //     config_.Nz / 2,                  // z position
    //     config_.Nx / 4,                  // x start
    //     3 * config_.Nx / 4               // x end
    // ) );

    grid_.add_source( std::make_unique<Point_Source>(
        10.0,
        config_.Nx / 2,
        config_.Ny / 2,
        config_.Nz / 2
    ) );

    // grid_.add_source( std::make_unique<Gaussian_Pulse>(
    //     10.0,                   // Amplitude
    //     30 * 3 * grid_.dt(),     // Center Time
    //     30 * grid_.dt(),         // Width
    //     config_.Nx / 2,          // x
    //     config_.Ny / 2,          // y
    //     config_.Nz / 2           // z
    // ) );
}

void Simulation::run() {
    std::cout << "<-----Maxwell Simulation----->" << std::endl;

    // Run @ t = 0:
    grid_.apply_sources();
    grid_.step();

    // Run simulation and start timer:
    std::size_t const output__interval{ config_.output_interval() };
    auto const start_time{ std::chrono::high_resolution_clock::now() };

    // Simulation Loop:
    for ( std::size_t curr_time{1}; curr_time <= config_.total_time; ++curr_time ) {
        grid_.apply_sources( curr_time );
        grid_.step();

        if ( ( curr_time % config_.output_interval() ) == 0 ) {
            output_.write_field( grid_, curr_time );
            grid_.print_progress( curr_time, config_.total_time );
        }
    }

    // End Timer:
    auto const end_time{ std::chrono::high_resolution_clock::now() };
    auto const duration{ std::chrono::duration_cast<std::chrono::milliseconds>( end_time - start_time ) };

    // Report Results:  
    std::cout << "\n\nResults:\n";
    std::cout << "--------\n";
    std::cout << "Duration: " << duration.count() << " ms\n";
    std::cout << "Physical time: " << config_.total_time * grid_.dt() << " s\n";
}