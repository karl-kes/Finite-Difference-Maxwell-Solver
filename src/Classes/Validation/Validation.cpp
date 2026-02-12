#include "Validation.hpp"

Plane_Wave_Test::Plane_Wave_Test( Simulation_Config const &cfg )
: config_{ cfg }
, grid_{ config_ }
, output_{ "validation_output" }
, wavelength_{ compute_wavelength( grid_, config_ ) }
, wavenumber_{ 2.0 * config::PI / wavelength_ }
, probe_x_{ grid_.Nx() / 2 }
, probe_y_{ grid_.Ny() / 2 }
, probe_z_{ grid_.Nz() / 2 }
, initial_phase_{ wavenumber_ * probe_x_ * grid_.dx() }
, phase_shift_per_step_{ wavenumber_ * grid_.c() * grid_.dt() }
{ }

void Plane_Wave_Test::initialize() {
    std::size_t const margin{ grid_.Nx() / 10 };

    #pragma omp parallel for collapse( 3 )
    for ( std::size_t z = margin; z < grid_.Nz() - margin; ++z ) {
        for ( std::size_t y = margin; y < grid_.Ny() - margin; ++y ) {
            for ( std::size_t x = margin; x < grid_.Nx() - margin; ++x ) {
                double const phase{ wavenumber_ * x * grid_.dx() };
                grid_.Ey( x, y, z ) = std::sin( phase );
                grid_.Bz( x, y, z ) = std::sin( phase ) / grid_.c();
            }
        }
    }
}

Validation_Result Plane_Wave_Test::run( std::size_t const num_steps ) {
    initialize();

    double const initial_energy{ grid_.total_energy() };

    double sum_expected{};
    double sum_actual{};
    double sum_expected_sq{};
    double sum_actual_sq{};
    double sum_product{};
    double total_phase_error{};

    for ( std::size_t t = 0; t < num_steps; ++t ) {
        double const expected_Ey{ std::sin( initial_phase_ - phase_shift_per_step_ * t ) };
        double const actual_Ey{ grid_.Ey( probe_x_, probe_y_, probe_z_ ) };

        sum_expected += expected_Ey;
        sum_actual += actual_Ey;
        sum_expected_sq += expected_Ey * expected_Ey;
        sum_actual_sq += actual_Ey * actual_Ey;
        sum_product += expected_Ey * actual_Ey;
        total_phase_error += std::abs( actual_Ey - expected_Ey );

        grid_.step();
    }

    double const final_energy{ grid_.total_energy() };

    // Calculate Metrics:
    double const energy_drift{ 100.0 * std::abs( final_energy - initial_energy ) / initial_energy };

    double const n{ static_cast<double>( num_steps ) };
    double const correlation_num{ n * sum_product - sum_expected * sum_actual };
    double const correlation_den{ 
        std::sqrt( ( n * sum_expected_sq - sum_expected * sum_expected ) *
                   ( n * sum_actual_sq - sum_actual * sum_actual ) ) 
    };
    double const correlation{ correlation_den > 1e-10 ? ( correlation_num / correlation_den ) : 0.0 };

    // Dispersion: average phase error relative to expected phase shift
    double const avg_phase_error{ total_phase_error / n };
    double const total_expected_shift{ phase_shift_per_step_ * n };
    double const dispersion{ 100.0 * avg_phase_error / std::max( total_expected_shift, 1e-10 ) };

    // Pass criteria
    bool const passed{ 
        energy_drift < 5.0 &&      // <5% energy drift
        correlation > 0.99 &&      // >99% correlation
        dispersion < 10.0          // <10% dispersion
    };

    return Validation_Result{ passed, energy_drift, dispersion, correlation };
}

void Plane_Wave_Test::print_report( Validation_Result const &result ) const {
    std::cout << "\n<-----Plane Wave Validation Test----->\n";
    std::cout << "Grid:       " << config_.Nx << " x " << config_.Ny << " x " << config_.Nz << "\n";
    std::cout << "Wavelength: " << wavelength_ << " (" << wavelength_ / grid_.dx() << " cells)\n";
    std::cout << "CFL Factor: " << config_.cfl_factor << "\n\n";
    std::cout << "Results:\n";
    std::cout << "--------\n";
    std::cout << std::fixed << std::setprecision( 5 );
    std::cout << "Energy Drift:      " << result.energy_drift_percent << " %\n";
    std::cout << "Phase Correlation: " << result.phase_correlation << "\n";
    std::cout << "Dispersion Error:  " << result.dispersion_percent << " %\n";
    std::cout << "\n";
    std::cout << "Status: " << ( result.passed ? "PASSED" : "FAILED" ) << "\n";
    std::cout << std::endl;
}