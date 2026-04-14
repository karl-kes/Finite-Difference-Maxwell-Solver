#include "Classes/Config/config.hpp"
#include "Classes/Grid/grid.hpp"
#include "Classes/Source/source.hpp"
#include "Classes/Write_Output/output.hpp"
#include "Classes/PML/pml.hpp"
#include "Classes/Simulation/simulation.hpp"
#include "Classes/Validation/validation.hpp"

#include <iostream>
#include <string>
#include <stdexcept>

int main( int const argc, char const* const argv[] ) {
    // --- Parse --config and --help first; everything else is forwarded to
    //     Simulation_Config::apply_cli_overrides. ---
    std::string config_path{ "config.cfg" };

    for ( int i{ 1 }; i < argc; ++i ) {
        std::string const arg{ argv[i] };
        if ( arg == "-h" || arg == "--help" ) {
            Simulation_Config::print_cli_help( std::cout );
            return 0;
        }
        if ( arg == "--config" ) {
            if ( i + 1 >= argc ) {
                std::cerr << "Error: --config requires a path.\n";
                return 1;
            }
            config_path = argv[++i];
        }
    }

    // --- Load defaults + file, then apply CLI overrides, then recompute. ---
    Simulation_Config config{ Simulation_Config::from_file( config_path ) };

    std::vector<std::string> overridden{};
    try {
        overridden = config.apply_cli_overrides( argc, argv );
    } catch ( std::exception const &e ) {
        std::cerr << "Error parsing CLI: " << e.what() << "\n\n";
        Simulation_Config::print_cli_help( std::cerr );
        return 1;
    }

    if ( !overridden.empty() ) {
        config.compute_derived();
        std::cout << "[Config] CLI overrides applied:";
        for ( auto const &flag : overridden ) { std::cout << " " << flag; }
        std::cout << "\n";
    }

    Simulation sim{ config };

    if ( config.run_validation ) {
        Plane_Wave_Test test{ config };
        Validation_Result result{ test.run() };
        test.print_report( result );
    }

    sim.initialize();
    sim.run();

    return 0;
}