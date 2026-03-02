#include "Classes/Config/config.hpp"
#include "Classes/Grid/grid.hpp"
#include "Classes/Source/source.hpp"
#include "Classes/Write_Output/output.hpp"
#include "Classes/PML/PML.hpp"
#include "Classes/Simulation/Simulation.hpp"
#include "Classes/Validation/Validation.hpp"

/* 
    To compile and run:

    # DEBUG MODE:
    rm -r build
    cmake -B build -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Debug
    cmake --build build
    ./build/main.exe

    # RELEASE MODE:
    rm -r build
    cmake -B build -G "MinGW Makefiles"
    cmake --build build
    ./build/main.exe

    # RENDER:
    python ./src/render.py
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