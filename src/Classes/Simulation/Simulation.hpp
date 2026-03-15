#pragma once

#include "../Config/config.hpp"
#include "../Grid/grid.hpp"
#include "../Write_Output/output.hpp"

class Simulation {
private:
    void print_progress( std::size_t const current, std::size_t const total ) const;
    
public:
    Simulation_Config config_;
    Grid grid_;
    Output output_;

    explicit Simulation( Simulation_Config const &new_config );

    void initialize();
    void run();
};