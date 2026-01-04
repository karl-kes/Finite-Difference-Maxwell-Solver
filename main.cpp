#include "Classes/Grid/grid.hpp"
#include "Classes/Timer/timer.hpp"

int main() {
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

    static constexpr int elapsed_time{ 1000 };

    static constexpr std::size_t Nx{ 10 }, Ny{ 10 }, Nz{ 10 };
    static constexpr double inject{ 1000.0 };

    static constexpr char B_field{ 'B' };
    static constexpr char E_field{ 'E' };

    Grid grid{ Nx+1, Ny+1, Nz+1 };

    grid.create_directories();
    grid.inject_source( Nx/2, Ny/2, Nz/2, inject );

    double initial_energy{ grid.total_energy() };
    double initial_div_B{ grid.div_B() };

    double max_energy_drift{ initial_energy };
    double max_div_B{ initial_div_B };

    auto start{ std::chrono::high_resolution_clock::now() };
    for ( int t{}; t <= elapsed_time; ++t ) {
        grid.step();

        // Max between current and previous max.
        max_energy_drift = std::max( grid.total_energy(), max_energy_drift );
        max_div_B = std::max( grid.div_B(), max_div_B );

        if ( t % 10 == 0 ) {
            std::string t_sec{ std::to_string( t ) };
            grid.vector_volume( "output/E/E" + t_sec + ".bin", E_field );
            grid.vector_volume( "output/B/B" + t_sec + ".bin", B_field );
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    double drift{ 100.0 * ( max_energy_drift - initial_energy ) / initial_energy };

    grid.output_final_metrics( elapsed_time, duration, drift, max_div_B );

    return 0;
}