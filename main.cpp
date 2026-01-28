#include "Classes/Config/config.hpp"
#include "Classes/Grid/grid.hpp"
#include "Classes/Source/source.hpp"
#include "Classes/Write_Output/output.hpp"

/* 
    To compile and run.

    For No Parallel (NOT RECOMMENDED):
    g++ -std=c++17 main.cpp Classes/Grid/*.cpp Classes/Source/*.cpp Classes/Write_Output/*.cpp -o main.exe
    ./main.exe
    ./render.py

    For Parallel (RECOMMENDED):
    g++ -std=c++17 -O3 -march=native -fopenmp main.cpp Classes/Grid/*.cpp Classes/Source/*.cpp Classes/Write_Output/*.cpp -o main.exe
    ./main.exe
    ./render.py
*/

void print_progress( double current, double total ) {
    double percent{ 100.0 * current / total };
    std::cout << "\rProgress: " << percent << "%" << std::flush;
}

int main() {
    // Configure:
    Simulation_Config config{};

    // Initialize:
    Grid grid( config );
    Output output{ "output" };
    output.initialize();

    // Add sources:
    grid.add_source( std::make_unique<Straight_Wire_X>(
        100000.0,                          // amplitude
        10.0,                            // frequency
        config.Ny / 2,                  // y position
        config.Nz / 2,                  // z position
        config.Nx / 4,                  // x start
        3 * config.Nx / 4               // x end
    ) );

    grid.add_source( std::make_unique<Point_Source>(
        100.0,
        config.Nx / 2,
        config.Ny / 2,
        config.Nz / 2
    ) );
    
    grid.apply_sources();
    grid.step();

    // Track Energy:
    double initial_energy{ grid.total_energy() };
    double max_energy{ initial_energy };

    // Run simulation and start timer:
    std::size_t output_interval{ config.output_interval() };
    auto start_time{ std::chrono::high_resolution_clock::now() };

    for ( std::size_t curr_time{}; curr_time <= config.total_time; ++curr_time ) {
        // grid.apply_sources( curr_time );
        grid.step();

        max_energy = std::max( grid.total_energy(), max_energy );

        if ( ( curr_time % output_interval ) == 0 ) {
            output.write_field( grid, Field::ELECTRIC, curr_time );
            output.write_field( grid, Field::MAGNETIC, curr_time );
            print_progress( curr_time, config.total_time );
        }
    }

    // End Timer:
    auto end_time{ std::chrono::high_resolution_clock::now() };
    auto duration{ std::chrono::duration_cast<std::chrono::milliseconds>( end_time - start_time ) };

    // Report results
    double energy_drift{ ( initial_energy == 0.0 ) ? 0.0 : ( 100.0 * ( max_energy - initial_energy ) / initial_energy ) };

    std::cout << "\n\n";
    std::cout << "Simulation Complete\n";
    std::cout << "-------------------\n";
    std::cout << "Duration: " << duration.count() << " ms\n";
    std::cout << "Physical time: " << config.total_time * grid.dt() << " s\n";
    std::cout << "Max energy drift: " << energy_drift << "%\n" << std::endl;

    return 0;
}