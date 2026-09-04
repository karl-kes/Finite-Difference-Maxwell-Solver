#pragma once

#include "../config/config.hpp"
#include "../grid/grid.hpp"
#include "../write-output/output.hpp"

#include <cstddef>

struct Validation_Result {
    bool passed;
    double rms_error_percent;
    double impedance_error_percent;
    double max_abs_error;
};

class Plane_Wave_Test {
private:
    Simulation_Config config_;
    Grid grid_;
    Output output_;

    double wavelength_;
    double wavenumber_;
    double c_;
    double eta_;

    std::size_t probe_x_, probe_y_, probe_z_;
    double initial_phase_;
    double phase_shift_per_step_;

    static double compute_wavelength( Grid const &grid, Simulation_Config const &cfg ) {
        std::size_t const usable_cells{ grid.Nx() - 2 * cfg.pml_thickness };
        double const wave_cells{ std::max( 40.0, std::floor( static_cast<double>( usable_cells ) * 0.2 ) ) };
        return wave_cells * grid.dx();
    }

public:
    explicit Plane_Wave_Test( Simulation_Config const &cfg );

    void initialize();
    Validation_Result run( std::size_t num_steps = 0 );
    void print_report( Validation_Result const &result ) const;
};