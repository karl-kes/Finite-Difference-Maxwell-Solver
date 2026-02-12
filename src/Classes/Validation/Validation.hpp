#pragma once

#include "../Config/config.hpp"
#include "../Grid/grid.hpp"
#include "../Write_Output/output.hpp"

#include <cmath>
#include <iostream>
#include <iomanip>

struct Validation_Result {
    bool passed;
    double energy_drift_percent;
    double dispersion_percent;
    double phase_correlation;
};

class Plane_Wave_Test {
private:
    Simulation_Config config_;
    Grid grid_;
    Output output_;

    double wavelength_;
    double wavenumber_;

    std::size_t probe_x_, probe_y_, probe_z_;
    double initial_phase_;
    double phase_shift_per_step_;

    static double compute_wavelength( Grid const &grid, Simulation_Config const &cfg ) {
        std::size_t const usable_cells{ grid.Nx() - 2 * cfg.pml_thickness };
        double const wave_cells{ std::min( 20.0, usable_cells * 0.3 ) };
        return wave_cells * grid.dx();
    }

public:
    explicit Plane_Wave_Test( Simulation_Config const &cfg );

    void initialize();
    Validation_Result run( std::size_t const num_steps = 100 );
    void print_report( Validation_Result const &result ) const;
};