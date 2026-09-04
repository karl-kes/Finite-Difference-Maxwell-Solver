#include "validation.hpp"

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
, c_{ 1.0 / std::sqrt( config_.mu * config_.eps ) }
, eta_{ std::sqrt( config_.mu / config_.eps ) }
, probe_x_{ grid_.Nx() / 2 }
, probe_y_{ grid_.Ny() / 2 }
, probe_z_{ grid_.Nz() / 2 }
, initial_phase_{ wavenumber_ * static_cast<double>( probe_x_ ) * grid_.dx() }
, phase_shift_per_step_{ wavenumber_ * c_ * grid_.dt() }
{ }

void Plane_Wave_Test::initialize() {
    double* RESTRICT Ey{ grid_.Ey_ptr() };
    double* RESTRICT Hz{ grid_.Hz_ptr() };

    std::size_t const Nx_local{ grid_.Nx() };
    std::size_t const Ny_local{ grid_.Ny() };
    std::size_t const Nz_local{ grid_.Nz() };

    double const inv_eta{ 1.0 / eta_ };
    double const dx_local{ grid_.dx() };
    
    // Use the 3D Yee discrete dispersion relation to calculate the exact numerical frequency
    double const cour{ c_ * grid_.dt() / dx_local };
    double const arg{ cour * std::sin( wavenumber_ * dx_local * 0.5 ) };
    double const num_omega_dt{ 2.0 * std::asin( arg ) };

    // The spatial shift places Hz a half-cell ahead. 
    // The temporal shift places Hz a half-step backward in time, using the exact numerical frequency.
    double const spatial_phase_shift{ wavenumber_ * dx_local * 0.5 };
    double const temporal_phase_shift{ num_omega_dt * 0.5 };
    double const total_hz_shift{ spatial_phase_shift + temporal_phase_shift };

    // Initialize the entire grid. No margins.
    #pragma omp parallel for collapse( 2 )
    for ( std::size_t z = 0; z < Nz_local; ++z ) {
        for ( std::size_t y = 0; y < Ny_local; ++y ) {

            #pragma omp simd
            for ( std::size_t x = 0; x < Nx_local; ++x ) {
                double const phase_x{ wavenumber_ * static_cast<double>( x ) * dx_local };
                std::size_t const i{ grid_.idx(x,y,z) };

                Ey[i] = std::sin( phase_x );
                Hz[i] = std::sin( phase_x + total_hz_shift ) * inv_eta;
            }
        }
    }
}

Validation_Result Plane_Wave_Test::run( std::size_t num_steps ) {
    initialize();

    if ( num_steps == 0 ) {
        // Run until right before the fastest boundary artifact hits the probe
        double const safe_dist{ static_cast<double>( grid_.Nx() / 2 ) * grid_.dx() * 0.8 };
        num_steps = static_cast<std::size_t>( safe_dist / c_ / grid_.dt() );
        num_steps = std::max( std::size_t{5}, num_steps );
    }

    std::size_t const probe_i{ grid_.idx(probe_x_, probe_y_, probe_z_) };
    double const* RESTRICT Ey{ grid_.Ey_ptr() };
    double const* RESTRICT Hz{ grid_.Hz_ptr() };

    double sum_sq_error{};
    double sum_sq_expected{};
    
    double max_Ey_amplitude{};
    double max_Hz_amplitude{};
    double max_abs_error{};

    // Calculate the expected numerical phase step from the discrete dispersion relation
    double const cour{ c_ * grid_.dt() / grid_.dx() };
    double const num_omega_dt{ 2.0 * std::asin( cour * std::sin( wavenumber_ * grid_.dx() * 0.5 ) ) };

    for ( std::size_t t = 1; t <= num_steps; ++t ) {
        grid_.step();
        
        // Compare the simulated wave against the exact mathematical expectation for a discrete Yee grid
        double const expected_Ey{ std::sin( initial_phase_ - num_omega_dt * static_cast<double>( t ) ) };
        double const actual_Ey{ Ey[probe_i] };
        double const actual_Hz{ Hz[probe_i] };

        double const diff{ std::abs( actual_Ey - expected_Ey ) };
        sum_sq_error += diff * diff;
        sum_sq_expected += expected_Ey * expected_Ey;
        
        max_Ey_amplitude = std::max( max_Ey_amplitude, std::abs( actual_Ey ) );
        max_Hz_amplitude = std::max( max_Hz_amplitude, std::abs( actual_Hz ) );
        max_abs_error = std::max( max_abs_error, diff );
    }

    double const rms_error{ std::sqrt( sum_sq_error / static_cast<double>( num_steps ) ) };
    double const rms_expected{ std::sqrt( sum_sq_expected / static_cast<double>( num_steps ) ) };
    double const rms_percent{ ( rms_expected > 1e-10 ) ? ( 100.0 * rms_error / rms_expected ) : 100.0 };

    double impedance_percent{ 100.0 };
    if ( max_Hz_amplitude > 1e-10 ) {
        double const simulated_eta{ max_Ey_amplitude / max_Hz_amplitude };
        impedance_percent = 100.0 * std::abs( simulated_eta - eta_ ) / eta_;
    }

    bool const passed{ 
        rms_percent < 5.0 && 
        impedance_percent < 5.0 
    };

    return Validation_Result{ passed, rms_percent, impedance_percent, max_abs_error };
}

void Plane_Wave_Test::print_report( Validation_Result const &result ) const {
    std::cout << "\n<--- Plane Wave Validation Test --->\n";
    std::cout << "Grid:       " << config_.Nx << " x " << config_.Ny << " x " << config_.Nz << "\n";
    std::cout << "Wavelength: " << wavelength_ << " (" << wavelength_ / grid_.dx() << " cells)\n";
    std::cout << "CFL Factor: " << config_.cfl_factor << "\n";
    std::cout << "dt:         " << grid_.dt() << "\n\n";
    std::cout << "Results:\n";
    std::cout << std::fixed << std::setprecision( 5 );
    std::cout << "RMS Wave Error:    " << result.rms_error_percent << " %\n";
    std::cout << "Impedance Error:   " << result.impedance_error_percent << " %\n";
    std::cout << "Max Abs Error:     " << result.max_abs_error << "\n";
    std::cout << "\n";
    std::cout << "Status: " << ( result.passed ? "PASSED" : "FAILED" ) << "\n";
    std::cout << std::endl;
}