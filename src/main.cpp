#include "Classes/Config/config.hpp"
#include "Classes/Grid/grid.hpp"
#include "Classes/Source/source.hpp"
#include "Classes/Write_Output/output.hpp"
#include "Classes/PML/PML.hpp"
#include "Classes/Simulation/Simulation.hpp"
#include "Classes/Validation/Validation.hpp"

/* 
    To compile and run:

    g++ -std=c++17 -O3 -march=native -fopenmp src/main.cpp src/Classes/Grid/*.cpp src/Classes/Source/*.cpp src/Classes/Write_Output/*.cpp src/Classes/PML/*.cpp src/Classes/Simulation/*.cpp src/Classes/Validation/*.cpp -o main.exe
    ./main.exe
    python src/render.py

    rm -r build
    mkdir build
    cd build
    cmake .. -G "MinGW Makefiles"
    mingw32-make -j
    .\main.exe
    python ../src/render.py
*/

int main() {
    Simulation_Config config{};
    Simulation sim{ config };

    // Validation:
    if ( config.run_validation ) {
        Plane_Wave_Test test{ config };
        Validation_Result result{ test.run() };
        test.print_report( result );
    }

    // Simulation:
    sim.initialize();
    sim.run();

    return 0;
}