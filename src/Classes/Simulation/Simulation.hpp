#pragma once

#include "../Config/config.hpp"
#include "../Grid/grid.hpp"
#include "../Write_Output/output.hpp"

class Simulation {
private:
    void print_progress( double const current, double const total ) const;
    
public:
    Simulation_Config config_;
    Grid grid_;
    Output output_;

    explicit Simulation( Simulation_Config const &new_config );

    void initialize();
    void run();
};