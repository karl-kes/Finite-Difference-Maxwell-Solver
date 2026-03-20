#include "Classes/Config/config.hpp"
#include "Classes/Grid/grid.hpp"
#include "Classes/Source/source.hpp"
#include "Classes/Write_Output/output.hpp"
#include "Classes/PML/PML.hpp"
#include "Classes/Simulation/Simulation.hpp"
#include "Classes/Validation/Validation.hpp"

int main() {
    // Config setup:
    std::string const config_path{ "config.cfg" };
    Simulation_Config config{ Simulation_Config::from_file( config_path ) };

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