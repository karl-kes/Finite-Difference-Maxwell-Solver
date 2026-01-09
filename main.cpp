#include "Classes/Grid/grid.hpp"

/* 
    To compile and run.

    For No Parallel (NOT RECOMMENDED):
    g++ -std=c++17 main.cpp Classes/Grid/grid.cpp -o main.exe
    ./main.exe
    ./render.py

    For Parallel (RECOMMENDED):
    g++ -std=c++17 main.cpp Classes/Grid/grid.cpp -o main.exe -fopenmp
    ./main.exe
    ./render.py
*/

int main() {
    Grid grid{ constant::Nx+1, constant::Ny+1, constant::Nz+1 };
    grid.create_directories();

    grid.hard_source_inject( constant::Nx/2, constant::Ny/2, constant::Nz/2, constant::inject );

    double initial_energy{ grid.total_energy() };
    double max_energy_drift{ initial_energy };

    auto start{ std::chrono::high_resolution_clock::now() };
    for ( int curr_time{}; curr_time <= constant::elapsed_time; ++curr_time ) {
        grid.step();

        max_energy_drift = std::max( grid.total_energy(), max_energy_drift );

        if ( curr_time % 10 == 0 ) {
            grid.print_progress( curr_time, constant::elapsed_time );
            std::string t_sec{ std::to_string( curr_time ) };
            grid.vector_volume( "output/E/E" + t_sec + ".bin", constant::E_field );
            grid.vector_volume( "output/B/B" + t_sec + ".bin", constant::B_field );
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    max_energy_drift = ( 100.0 * ( max_energy_drift - initial_energy ) / initial_energy );

    std::cout << std::endl;
    std::cout << "Duration of Simulation: " << duration.count() << " ms\n" << std::endl;
    std::cout << "Physical Time Simulated: " << constant::elapsed_time * grid.dt() << " s" << std::endl;
    std::cout << "Max Energy Drift: " << max_energy_drift << "%" << std::endl;

    return 0;
}