#pragma once

#include "../Config/config.hpp"
#include "../Grid/grid.hpp"
#include "../Source/source.hpp"
#include "../Write_Output/output.hpp"
#include "../PML/PML.hpp"

class Output;
class Simulation_Config;

class Simulation {
public:
    Simulation_Config config_;
    Grid grid_;
    Output output_;

    explicit Simulation( Simulation_Config const &new_config );

    void initialize();
    void run();
};