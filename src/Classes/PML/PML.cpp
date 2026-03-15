#include "PML.hpp"

#include <omp.h>

PML::PML( Simulation_Config const &config )
: thickness_{ config.use_pml ? config.pml_thickness : 0 }
, Nx_{ config.Nx + 1 }
, Ny_{ config.Ny + 1 }
, Nz_{ config.Nz + 1 }
, Nx_padded_{ AlignedSoA<double>::round_up( Nx_ ) }
, Ny_padded_{ AlignedSoA<double>::round_up( Ny_ ) }
, Nz_padded_{ AlignedSoA<double>::round_up( Nz_ ) }
, order_{ config.pml_order }
, sigma_max_{ config.pml_sigma_max }
, kappa_max_{ config.pml_kappa_max }
, alpha_max_{ config.pml_alpha_max }
, coeffs_{}
, psi_{}
, psi_face_x_{ thickness_ * Ny_ * Nz_ }
, psi_face_y_{ Nx_ * thickness_ * Nz_ }
, psi_face_z_{ Nx_ * Ny_ * thickness_ }
{
    if ( !config.use_pml || thickness_ == 0 ) return;

    // Coefficient block: 18 arrays, each of length thickness (padded by AlignedSoA):
    coeffs_ = AlignedSoA<double>{ thickness_, NUM_COEFF_ARRAYS_ };

    // Psi block: 12 arrays, each needs 2 * max_face_size elements (lo + hi face):
    std::size_t const max_face{ std::max({ psi_face_x_, psi_face_y_, psi_face_z_ }) };
    psi_ = AlignedSoA<double>{ 2 * max_face, NUM_PSI_ARRAYS_ };

    // Compute grading coefficients for each PML layer:
    double const d{ static_cast<double>( thickness_ ) };

    double* RESTRICT b_Ex{ b_Ex_ptr() }; double* RESTRICT c_Ex{ c_Ex_ptr() }; double* RESTRICT kappa_Ex{ kappa_Ex_ptr() };
    double* RESTRICT b_Bx{ b_Bx_ptr() }; double* RESTRICT c_Bx{ c_Bx_ptr() }; double* RESTRICT kappa_Bx{ kappa_Bx_ptr() };
    
    double* RESTRICT b_Ey{ b_Ey_ptr() }; double* RESTRICT c_Ey{ c_Ey_ptr() }; double* RESTRICT kappa_Ey{ kappa_Ey_ptr() };
    double* RESTRICT b_By{ b_By_ptr() }; double* RESTRICT c_By{ c_By_ptr() }; double* RESTRICT kappa_By{ kappa_By_ptr() };
    
    double* RESTRICT b_Ez{ b_Ez_ptr() }; double* RESTRICT c_Ez{ c_Ez_ptr() }; double* RESTRICT kappa_Ez{ kappa_Ez_ptr() };
    double* RESTRICT b_Bz{ b_Bz_ptr() }; double* RESTRICT c_Bz{ c_Bz_ptr() }; double* RESTRICT kappa_Bz{ kappa_Bz_ptr() };

    for ( std::size_t i{ 0 }; i < thickness_; ++i ) {
        double const depth_E{ ( d - static_cast<double>( i ) ) / d };
        double const sig_E{ sigma( depth_E ) };
        double const kap_E{ kappa( depth_E ) };
        double const alp_E{ alpha( depth_E ) };

        double const depth_B{ ( d - ( static_cast<double>( i ) + 0.5 ) ) / d };
        double const sig_B{ sigma( depth_B ) };
        double const kap_B{ kappa( depth_B ) };
        double const alp_B{ alpha( depth_B ) };

        compute_coefficients( sig_E, kap_E, alp_E, config.dt, config.eps, b_Ex[i], c_Ex[i] );
        compute_coefficients( sig_B, kap_B, alp_B, config.dt, config.eps, b_Bx[i], c_Bx[i] );
        kappa_Ex[i] = kap_E;
        kappa_Bx[i] = kap_B;

        compute_coefficients( sig_E, kap_E, alp_E, config.dt, config.eps, b_Ey[i], c_Ey[i] );
        compute_coefficients( sig_B, kap_B, alp_B, config.dt, config.eps, b_By[i], c_By[i] );
        kappa_Ey[i] = kap_E;
        kappa_By[i] = kap_B;

        compute_coefficients( sig_E, kap_E, alp_E, config.dt, config.eps, b_Ez[i], c_Ez[i] );
        compute_coefficients( sig_B, kap_B, alp_B, config.dt, config.eps, b_Bz[i], c_Bz[i] );
        kappa_Ez[i] = kap_E;
        kappa_Bz[i] = kap_B;
    }
}

void PML::update_B_psi(
    double* RESTRICT Ex, double* RESTRICT Ey, double* RESTRICT Ez,
    double* RESTRICT Bx, double* RESTRICT By, double* RESTRICT Bz,
    double const dt, double const dx, double const dy, double const dz ) {
    if ( !is_active() ) return;

    std::size_t const t{ thickness_ };
    std::size_t const face_x{ psi_face_x_ };
    std::size_t const face_y{ psi_face_y_ };
    std::size_t const face_z{ psi_face_z_ };

    double const* RESTRICT b_Bx{ b_Bx_ptr() };
    double const* RESTRICT c_Bx{ c_Bx_ptr() };
    double const* RESTRICT b_By{ b_By_ptr() };
    double const* RESTRICT c_By{ c_By_ptr() };
    double const* RESTRICT b_Bz{ b_Bz_ptr() };
    double const* RESTRICT c_Bz{ c_Bz_ptr() };

    double* RESTRICT psi_Byx{ psi_Byx_ptr() };
    double* RESTRICT psi_Bzx{ psi_Bzx_ptr() };
    double* RESTRICT psi_Bxy{ psi_Bxy_ptr() };
    double* RESTRICT psi_Bzy{ psi_Bzy_ptr() };
    double* RESTRICT psi_Bxz{ psi_Bxz_ptr() };
    double* RESTRICT psi_Byz{ psi_Byz_ptr() };

    // x-faces (low and high):
    #pragma omp parallel for collapse( 2 )
    for ( std::size_t z = 0; z < Nz_ - 1; ++z ) {
        for ( std::size_t y = 0; y < Ny_ - 1; ++y ) {
            for ( std::size_t d = 0; d < t; ++d ) {
                // Low x-face:
                std::size_t const x_lo{ d };
                std::size_t const gi_lo{ idx( x_lo, y, z ) };
                std::size_t const pi_lo{ psi_idx_x( d, y, z ) };

                double const dEy_dx_lo{ ( Ey[idx( x_lo + 1, y, z )] - Ey[gi_lo] ) / dx };
                double const dEz_dx_lo{ ( Ez[idx( x_lo + 1, y, z )] - Ez[gi_lo] ) / dx };

                psi_Byx[pi_lo] = b_Bx[d] * psi_Byx[pi_lo] + c_Bx[d] * dEy_dx_lo;
                psi_Bzx[pi_lo] = b_Bx[d] * psi_Bzx[pi_lo] + c_Bx[d] * dEz_dx_lo;

                Bz[gi_lo] -= dt * psi_Byx[pi_lo];
                By[gi_lo] += dt * psi_Bzx[pi_lo];

                // High x-face:
                std::size_t const x_hi{ Nx_ - 1 - t + d };
                std::size_t const gi_hi{ idx( x_hi, y, z ) };
                std::size_t const pi_hi{ face_x + psi_idx_x( d, y, z ) };

                double const dEy_dx_hi{ ( Ey[idx( x_hi + 1, y, z )] - Ey[gi_hi] ) / dx };
                double const dEz_dx_hi{ ( Ez[idx( x_hi + 1, y, z )] - Ez[gi_hi] ) / dx };

                psi_Byx[pi_hi] = b_Bx[t - 1 - d] * psi_Byx[pi_hi] + c_Bx[t - 1 - d] * dEy_dx_hi;
                psi_Bzx[pi_hi] = b_Bx[t - 1 - d] * psi_Bzx[pi_hi] + c_Bx[t - 1 - d] * dEz_dx_hi;

                Bz[gi_hi] -= dt * psi_Byx[pi_hi];
                By[gi_hi] += dt * psi_Bzx[pi_hi];
            }
        }
    }

    // y-faces (low and high):
    #pragma omp parallel for collapse( 2 )
    for ( std::size_t z = 0; z < Nz_ - 1; ++z ) {
        for ( std::size_t x = 0; x < Nx_ - 1; ++x ) {
            for ( std::size_t d = 0; d < t; ++d ) {
                // Low y-face:
                std::size_t const y_lo{ d };
                std::size_t const gi_lo{ idx( x, y_lo, z ) };
                std::size_t const pi_lo{ psi_idx_y( x, d, z ) };

                double const dEx_dy_lo{ ( Ex[idx( x, y_lo + 1, z )] - Ex[gi_lo] ) / dy };
                double const dEz_dy_lo{ ( Ez[idx( x, y_lo + 1, z )] - Ez[gi_lo] ) / dy };

                psi_Bxy[pi_lo] = b_By[d] * psi_Bxy[pi_lo] + c_By[d] * dEx_dy_lo;
                psi_Bzy[pi_lo] = b_By[d] * psi_Bzy[pi_lo] + c_By[d] * dEz_dy_lo;

                Bz[gi_lo] += dt * psi_Bxy[pi_lo];
                Bx[gi_lo] -= dt * psi_Bzy[pi_lo];

                // High y-face:
                std::size_t const y_hi{ Ny_ - 1 - t + d };
                std::size_t const gi_hi{ idx( x, y_hi, z ) };
                std::size_t const pi_hi{ face_y + psi_idx_y( x, d, z ) };

                double const dEx_dy_hi{ ( Ex[idx( x, y_hi + 1, z )] - Ex[gi_hi] ) / dy };
                double const dEz_dy_hi{ ( Ez[idx( x, y_hi + 1, z )] - Ez[gi_hi] ) / dy };

                psi_Bxy[pi_hi] = b_By[t - 1 - d] * psi_Bxy[pi_hi] + c_By[t - 1 - d] * dEx_dy_hi;
                psi_Bzy[pi_hi] = b_By[t - 1 - d] * psi_Bzy[pi_hi] + c_By[t - 1 - d] * dEz_dy_hi;

                Bz[gi_hi] += dt * psi_Bxy[pi_hi];
                Bx[gi_hi] -= dt * psi_Bzy[pi_hi];
            }
        }
    }

    // z-faces (low and high):
    #pragma omp parallel for collapse( 2 )
    for ( std::size_t y = 0; y < Ny_ - 1; ++y ) {
        for ( std::size_t x = 0; x < Nx_ - 1; ++x ) {
            for ( std::size_t d = 0; d < t; ++d ) {
                // Low z-face:
                std::size_t const z_lo{ d };
                std::size_t const gi_lo{ idx( x, y, z_lo ) };
                std::size_t const pi_lo{ psi_idx_z( x, y, d ) };

                double const dEx_dz_lo{ ( Ex[idx( x, y, z_lo + 1 )] - Ex[gi_lo] ) / dz };
                double const dEy_dz_lo{ ( Ey[idx( x, y, z_lo + 1 )] - Ey[gi_lo] ) / dz };

                psi_Bxz[pi_lo] = b_Bz[d] * psi_Bxz[pi_lo] + c_Bz[d] * dEx_dz_lo;
                psi_Byz[pi_lo] = b_Bz[d] * psi_Byz[pi_lo] + c_Bz[d] * dEy_dz_lo;

                Bx[gi_lo] += dt * psi_Byz[pi_lo];
                By[gi_lo] -= dt * psi_Bxz[pi_lo];

                // High z-face:
                std::size_t const z_hi{ Nz_ - 1 - t + d };
                std::size_t const gi_hi{ idx( x, y, z_hi ) };
                std::size_t const pi_hi{ face_z + psi_idx_z( x, y, d ) };

                double const dEx_dz_hi{ ( Ex[idx( x, y, z_hi + 1 )] - Ex[gi_hi] ) / dz };
                double const dEy_dz_hi{ ( Ey[idx( x, y, z_hi + 1 )] - Ey[gi_hi] ) / dz };

                psi_Bxz[pi_hi] = b_Bz[t - 1 - d] * psi_Bxz[pi_hi] + c_Bz[t - 1 - d] * dEx_dz_hi;
                psi_Byz[pi_hi] = b_Bz[t - 1 - d] * psi_Byz[pi_hi] + c_Bz[t - 1 - d] * dEy_dz_hi;

                Bx[gi_hi] += dt * psi_Byz[pi_hi];
                By[gi_hi] -= dt * psi_Bxz[pi_hi];
            }
        }
    }
}

void PML::update_E_psi(
    double* RESTRICT Ex, double* RESTRICT Ey, double* RESTRICT Ez,
    double* RESTRICT Bx, double* RESTRICT By, double* RESTRICT Bz,
    double const dt, double const dx, double const dy, double const dz,
    double const c_sq ) {
    if ( !is_active() ) return;

    std::size_t const t{ thickness_ };
    std::size_t const face_x{ psi_face_x_ };
    std::size_t const face_y{ psi_face_y_ };
    std::size_t const face_z{ psi_face_z_ };

    double const* RESTRICT b_Ex{ b_Ex_ptr() };
    double const* RESTRICT c_Ex{ c_Ex_ptr() };
    double const* RESTRICT b_Ey{ b_Ey_ptr() };
    double const* RESTRICT c_Ey{ c_Ey_ptr() };
    double const* RESTRICT b_Ez{ b_Ez_ptr() };
    double const* RESTRICT c_Ez{ c_Ez_ptr() };

    double* RESTRICT psi_Eyx{ psi_Eyx_ptr() };
    double* RESTRICT psi_Ezx{ psi_Ezx_ptr() };
    double* RESTRICT psi_Exy{ psi_Exy_ptr() };
    double* RESTRICT psi_Ezy{ psi_Ezy_ptr() };
    double* RESTRICT psi_Exz{ psi_Exz_ptr() };
    double* RESTRICT psi_Eyz{ psi_Eyz_ptr() };

    // x-faces (low and high):
    #pragma omp parallel for collapse( 2 )
    for ( std::size_t z = 1; z < Nz_; ++z ) {
        for ( std::size_t y = 1; y < Ny_; ++y ) {
            for ( std::size_t d = 0; d < t; ++d ) {
                // Low x-face (starts at x=1 for E update):
                std::size_t const x_lo{ d + 1 };
                if ( x_lo >= Nx_ ) continue;
                std::size_t const gi_lo{ idx( x_lo, y, z ) };
                std::size_t const pi_lo{ psi_idx_x( d, y, z ) };

                double const dBy_dx_lo{ ( By[gi_lo] - By[idx( x_lo - 1, y, z )] ) / dx };
                double const dBz_dx_lo{ ( Bz[gi_lo] - Bz[idx( x_lo - 1, y, z )] ) / dx };

                psi_Eyx[pi_lo] = b_Ex[d] * psi_Eyx[pi_lo] + c_Ex[d] * dBy_dx_lo;
                psi_Ezx[pi_lo] = b_Ex[d] * psi_Ezx[pi_lo] + c_Ex[d] * dBz_dx_lo;

                Ez[gi_lo] += dt * c_sq * psi_Eyx[pi_lo];
                Ey[gi_lo] -= dt * c_sq * psi_Ezx[pi_lo];

                // High x-face:
                std::size_t const x_hi{ Nx_ - t + d };
                if ( x_hi >= Nx_ ) continue;
                std::size_t const gi_hi{ idx( x_hi, y, z ) };
                std::size_t const pi_hi{ face_x + psi_idx_x( d, y, z ) };

                double const dBy_dx_hi{ ( By[gi_hi] - By[idx( x_hi - 1, y, z )] ) / dx };
                double const dBz_dx_hi{ ( Bz[gi_hi] - Bz[idx( x_hi - 1, y, z )] ) / dx };

                psi_Eyx[pi_hi] = b_Ex[t - 1 - d] * psi_Eyx[pi_hi] + c_Ex[t - 1 - d] * dBy_dx_hi;
                psi_Ezx[pi_hi] = b_Ex[t - 1 - d] * psi_Ezx[pi_hi] + c_Ex[t - 1 - d] * dBz_dx_hi;

                Ez[gi_hi] += dt * c_sq * psi_Eyx[pi_hi];
                Ey[gi_hi] -= dt * c_sq * psi_Ezx[pi_hi];
            }
        }
    }

    // y-faces (low and high):
    #pragma omp parallel for collapse( 2 )
    for ( std::size_t z = 1; z < Nz_; ++z ) {
        for ( std::size_t x = 1; x < Nx_; ++x ) {
            for ( std::size_t d = 0; d < t; ++d ) {
                // Low y-face:
                std::size_t const y_lo{ d + 1 };
                if ( y_lo >= Ny_ ) continue;
                std::size_t const gi_lo{ idx( x, y_lo, z ) };
                std::size_t const pi_lo{ psi_idx_y( x, d, z ) };

                double const dBx_dy_lo{ ( Bx[gi_lo] - Bx[idx( x, y_lo - 1, z )] ) / dy };
                double const dBz_dy_lo{ ( Bz[gi_lo] - Bz[idx( x, y_lo - 1, z )] ) / dy };

                psi_Exy[pi_lo] = b_Ey[d] * psi_Exy[pi_lo] + c_Ey[d] * dBx_dy_lo;
                psi_Ezy[pi_lo] = b_Ey[d] * psi_Ezy[pi_lo] + c_Ey[d] * dBz_dy_lo;

                Ez[gi_lo] -= dt * c_sq * psi_Exy[pi_lo];
                Ex[gi_lo] += dt * c_sq * psi_Ezy[pi_lo];

                // High y-face:
                std::size_t const y_hi{ Ny_ - t + d };
                if ( y_hi >= Ny_ ) continue;
                std::size_t const gi_hi{ idx( x, y_hi, z ) };
                std::size_t const pi_hi{ face_y + psi_idx_y( x, d, z ) };

                double const dBx_dy_hi{ ( Bx[gi_hi] - Bx[idx( x, y_hi - 1, z )] ) / dy };
                double const dBz_dy_hi{ ( Bz[gi_hi] - Bz[idx( x, y_hi - 1, z )] ) / dy };

                psi_Exy[pi_hi] = b_Ey[t - 1 - d] * psi_Exy[pi_hi] + c_Ey[t - 1 - d] * dBx_dy_hi;
                psi_Ezy[pi_hi] = b_Ey[t - 1 - d] * psi_Ezy[pi_hi] + c_Ey[t - 1 - d] * dBz_dy_hi;

                Ez[gi_hi] -= dt * c_sq * psi_Exy[pi_hi];
                Ex[gi_hi] += dt * c_sq * psi_Ezy[pi_hi];
            }
        }
    }

    // z-faces (low and high):
    #pragma omp parallel for collapse( 2 )
    for ( std::size_t y = 1; y < Ny_; ++y ) {
        for ( std::size_t x = 1; x < Nx_; ++x ) {
            for ( std::size_t d = 0; d < t; ++d ) {
                // Low z-face:
                std::size_t const z_lo{ d + 1 };
                if ( z_lo >= Nz_ ) continue;
                std::size_t const gi_lo{ idx( x, y, z_lo ) };
                std::size_t const pi_lo{ psi_idx_z( x, y, d ) };

                double const dBx_dz_lo{ ( Bx[gi_lo] - Bx[idx( x, y, z_lo - 1 )] ) / dz };
                double const dBy_dz_lo{ ( By[gi_lo] - By[idx( x, y, z_lo - 1 )] ) / dz };

                psi_Exz[pi_lo] = b_Ez[d] * psi_Exz[pi_lo] + c_Ez[d] * dBx_dz_lo;
                psi_Eyz[pi_lo] = b_Ez[d] * psi_Eyz[pi_lo] + c_Ez[d] * dBy_dz_lo;

                Ey[gi_lo] += dt * c_sq * psi_Exz[pi_lo];
                Ex[gi_lo] -= dt * c_sq * psi_Eyz[pi_lo];

                // High z-face:
                std::size_t const z_hi{ Nz_ - t + d };
                if ( z_hi >= Nz_ ) continue;
                std::size_t const gi_hi{ idx( x, y, z_hi ) };
                std::size_t const pi_hi{ face_z + psi_idx_z( x, y, d ) };

                double const dBx_dz_hi{ ( Bx[gi_hi] - Bx[idx( x, y, z_hi - 1 )] ) / dz };
                double const dBy_dz_hi{ ( By[gi_hi] - By[idx( x, y, z_hi - 1 )] ) / dz };

                psi_Exz[pi_hi] = b_Ez[t - 1 - d] * psi_Exz[pi_hi] + c_Ez[t - 1 - d] * dBx_dz_hi;
                psi_Eyz[pi_hi] = b_Ez[t - 1 - d] * psi_Eyz[pi_hi] + c_Ez[t - 1 - d] * dBy_dz_hi;

                Ey[gi_hi] += dt * c_sq * psi_Exz[pi_hi];
                Ex[gi_hi] -= dt * c_sq * psi_Eyz[pi_hi];
            }
        }
    }
}