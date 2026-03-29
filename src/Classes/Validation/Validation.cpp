#include "Validation.hpp"

#include <numbers>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>

Plane_Wave_Test::Plane_Wave_Test( Simulation_Config const &cfg )
: config_{ cfg }
, grid_{ config_ }
, output_{ "validation_output" }
, wavelength_{ compute_wavelength( grid_, config_ ) }
, wavenumber_{ 2.0 * std::numbers::pi / wavelength_ }
, probe_x_{ grid_.Nx() / 2 }
, probe_y_{ grid_.Ny() / 2 }
, probe_z_{ grid_.Nz() / 2 }
, initial_phase_{ wavenumber_ * static_cast<double>( probe_x_ ) * grid_.dx() }
, phase_shift_per_step_{ wavenumber_ * grid_.c() * grid_.dt() }
{ }

void Plane_Wave_Test::initialize() {
    std::size_t const margin{ grid_.Nx() / 8 };

    double* RESTRICT Ey{ grid_.Ey_ptr() };
    double* RESTRICT Bz{ grid_.Bz_ptr() };

    std::size_t Nx_end{grid_.Nx() - margin};
    std::size_t Ny_end{grid_.Ny() - margin};
    std::size_t Nz_end{grid_.Nz() - margin};

    double const inv_c{ 1.0 / grid_.c() };
    double const dx_local{ grid_.dx() };

    #pragma omp parallel for collapse( 2 )
    for ( std::size_t z = margin; z < Nz_end; ++z ) {
        for ( std::size_t y = margin; y < Ny_end; ++y ) {

            #pragma omp simd
            for ( std::size_t x = margin; x < Nx_end; ++x ) {
                double const phase{ wavenumber_ * static_cast<double>( x ) * dx_local };
                std::size_t const i{ grid_.idx(x,y,z) };

                Ey[i] = std::sin( phase );
                Bz[i] = std::sin( phase ) * inv_c;
            }
        }
    }
}

Validation_Result Plane_Wave_Test::run( std::size_t num_steps ) {
    initialize();

    if ( num_steps == 0 ) {
        double const margin_dist{ static_cast<double>( grid_.Nx() / 10 ) * grid_.dx() };
        double const travel_time{ 0.5 * margin_dist / grid_.c() };
        num_steps = std::max( std::size_t{10},
                              static_cast<std::size_t>( travel_time / grid_.dt() ) );
    }

    double const initial_energy{ grid_.total_energy() };

    double sum_expected{};
    double sum_actual{};
    double sum_expected_sq{};
    double sum_actual_sq{};
    double sum_product{};
    double total_phase_error{};

    std::size_t const probe_i{ grid_.idx(probe_x_, probe_y_, probe_z_) };
    double const* RESTRICT Ey{ grid_.Ey_ptr() };

    for ( std::size_t t = 0; t < num_steps; ++t ) {
        double const expected_Ey{ std::sin( initial_phase_ - phase_shift_per_step_ * static_cast<double>( t ) ) };
        double const actual_Ey{ Ey[probe_i] };

        sum_expected += expected_Ey;
        sum_actual += actual_Ey;
        sum_expected_sq += expected_Ey * expected_Ey;
        sum_actual_sq += actual_Ey * actual_Ey;
        sum_product += expected_Ey * actual_Ey;
        total_phase_error += std::abs( actual_Ey - expected_Ey );

        grid_.step();
    }

    double const final_energy{ grid_.total_energy() };

    double const energy_drift{ 100.0 * std::abs( final_energy - initial_energy ) / initial_energy };

    double const n{ static_cast<double>( num_steps ) };
    double const correlation_num{ n * sum_product - sum_expected * sum_actual };
    double const correlation_den{ 
        std::sqrt( ( n * sum_expected_sq - sum_expected * sum_expected ) *
                   ( n * sum_actual_sq - sum_actual * sum_actual ) ) 
    };
    double const correlation{ correlation_den > 1e-10 ? ( correlation_num / correlation_den ) : 0.0 };

    double const avg_phase_error{ total_phase_error / n };
    double const total_expected_shift{ phase_shift_per_step_ * n };
    double const dispersion{ 100.0 * avg_phase_error / std::max( total_expected_shift, 1e-10 ) };

    bool const passed{ 
        energy_drift < 5.0 &&
        correlation > 0.99 &&
        dispersion < 10.0
    };

    return Validation_Result{ passed, energy_drift, dispersion, correlation };
}

void Plane_Wave_Test::print_report( Validation_Result const &result ) const {
    std::cout << "\n<-----Plane Wave Validation Test----->\n";
    std::cout << "Grid:       " << config_.Nx << " x " << config_.Ny << " x " << config_.Nz << "\n";
    std::cout << "Wavelength: " << wavelength_ << " (" << wavelength_ / grid_.dx() << " cells)\n";
    std::cout << "CFL Factor: " << config_.cfl_factor << "\n";
    std::cout << "dt:         " << grid_.dt() << "\n\n";
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