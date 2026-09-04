#pragma once

#include "../config/config.hpp"
#include "../grid/grid.hpp"
#include "../write-output/output.hpp"

class Simulation {
private:
    void print_progress( std::size_t const current, std::size_t const total ) const;

    [[nodiscard]] Simulation_Config &config() { return config_; }
    [[nodiscard]] Simulation_Config const &config() const { return config_; }

    [[nodiscard]] Grid &grid() { return grid_; }
    [[nodiscard]] Grid const &grid() const { return grid_; }

    [[nodiscard]] Output &output() { return output_; }
    [[nodiscard]] Output const &output() const { return output_; }
    
public:
    Simulation_Config config_;
    Grid grid_;
    Output output_;

    explicit Simulation( Simulation_Config const &new_config );

    void initialize_sources();
    void initialize();
    void run();
};