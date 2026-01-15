#include "Classes/Grid/grid.hpp"

/* 
    To compile and run.

    For No Parallel (NOT RECOMMENDED):
    g++ -std=c++17 main.cpp Classes/Grid/grid_constructor.cpp Classes/Grid/grid_getters.cpp Classes/Grid/grid_simulation.cpp Classes/Grid/grid_helpers.cpp -o main.exe
    ./main.exe
    ./render.py

    For Parallel (RECOMMENDED):
    g++ -std=c++17 main.cpp Classes/Grid/grid_constructor.cpp Classes/Grid/grid_getters.cpp Classes/Grid/grid_simulation.cpp Classes/Grid/grid_helpers.cpp -o main.exe -fopenmp
    ./main.exe
    ./render.py
*/

int main() {
    Grid grid{ config::Nx+1, config::Ny+1, config::Nz+1,
               config::dx, config::dy, config::dz,
               config::eps, config::mu,  };
    grid.create_directories();

    // grid.hard_source_inject( config::inject,
    //                          config::Nx/2, config::Ny/2, config::Nz/2 );

    // grid.dipole_antenna_inject( config::amp_one, config::amp_two,
    //                             config::freq_one, config::freq_two,
    //                             config::inject,
    //                             config::Nx/4, config::Ny/4, config::Nz/4 );

    double initial_energy{ grid.total_energy() };
    double max_energy_drift{ initial_energy };

    auto start{ std::chrono::high_resolution_clock::now() };
    for ( std::size_t curr_time{}; curr_time <= config::total_time; ++curr_time ) {
        grid.straight_wire_x( config::amp_one, config::freq_one, curr_time, config::Ny/2, config::Nz/2 );
        grid.step();

        max_energy_drift = std::max( grid.total_energy(), max_energy_drift );

        if ( curr_time % config::print_rate == 0 ) {
            grid.print_progress( curr_time, config::total_time );
            std::string t_sec{ std::to_string( curr_time ) };
            grid.vector_volume( "output/E/E" + t_sec + ".bin", config::E_field );
            grid.vector_volume( "output/B/B" + t_sec + ".bin", config::B_field );
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    max_energy_drift = ( 100.0 * ( max_energy_drift - initial_energy ) / initial_energy );

    std::cout << std::endl;
    std::cout << "Duration of Simulation: " << duration.count() << " ms" << std::endl;
    std::cout << "Physical Time Simulated: " << config::total_time * grid.dt() << " s" << std::endl;
    std::cout << "Max Energy Drift: " << max_energy_drift << "%\n" << std::endl;

    return 0;
}