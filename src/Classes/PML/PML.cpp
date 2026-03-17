#include "PML.hpp"

#include <omp.h>
#include <algorithm>

PML::PML( Simulation_Config const &config )
: thickness_{ config.use_pml ? config.pml_thickness : 0 }
, Nx_{ config.Nx + 1 }
, Ny_{ config.Ny + 1 }
, Nz_{ config.Nz + 1 }
, Nx_padded_{ AlignedSoA<double>::round_up( Nx_ ) }
, Ny_padded_{ AlignedSoA<double>::round_up( Ny_ ) }
, order_{ config.pml_order }
, sigma_max_{ config.pml_sigma_max }
, kappa_max_{ config.pml_kappa_max }
, alpha_max_{ config.pml_alpha_max }
, coeffs_{}
, psi_{}
, psi_face_x_{ thickness_ * Ny_ * Nz_ }
, psi_face_y_{ Nx_ * thickness_ * Nz_ }
, psi_face_z_{ Nx_ * Ny_ * thickness_ } {
    if ( !config.use_pml || thickness_ == 0 ) return;

    coeffs_ = AlignedSoA<double>{ thickness_, NUM_COEFF_ARRAYS_ };

    std::size_t const max_face{ std::max({ psi_face_x_, psi_face_y_, psi_face_z_ }) };
    psi_ = AlignedSoA<double>{ 2 * max_face, NUM_PSI_ARRAYS_ };

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
        double const depth_B_clamped{ std::max( depth_B, 0.0 ) };
        double const sig_B{ sigma( depth_B_clamped ) };
        double const kap_B{ kappa( depth_B_clamped ) };
        double const alp_B{ alpha( depth_B_clamped ) };

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

// B-field PML

void PML::update_B_psi(
    double* RESTRICT Ex, double* RESTRICT Ey, double* RESTRICT Ez,
    double* RESTRICT Bx, double* RESTRICT By, double* RESTRICT Bz,
    double const dt, double const dx, double const dy, double const dz ) {
    if ( !is_active() ) return;

    std::size_t const t{ thickness_ };
    std::size_t const face_x{ psi_face_x_ };
    std::size_t const face_y{ psi_face_y_ };
    std::size_t const face_z{ psi_face_z_ };

    std::size_t const Sx{ 1 };
    std::size_t const Sy{ Nx_padded_ };
    std::size_t const Sz{ Nx_padded_ * Ny_padded_ };

    double const inv_dx{ 1.0 / dx };
    double const inv_dy{ 1.0 / dy };
    double const inv_dz{ 1.0 / dz };

    double const* RESTRICT bBx{ b_Bx_ptr() };
    double const* RESTRICT cBx{ c_Bx_ptr() };
    double const* RESTRICT bBy{ b_By_ptr() };
    double const* RESTRICT cBy{ c_By_ptr() };
    double const* RESTRICT bBz{ b_Bz_ptr() };
    double const* RESTRICT cBz{ c_Bz_ptr() };

    double* RESTRICT pByx{ psi_Byx_ptr() };
    double* RESTRICT pBzx{ psi_Bzx_ptr() };
    double* RESTRICT pBxy{ psi_Bxy_ptr() };
    double* RESTRICT pBzy{ psi_Bzy_ptr() };
    double* RESTRICT pBxz{ psi_Bxz_ptr() };
    double* RESTRICT pByz{ psi_Byz_ptr() };

    std::size_t const nx{ Nx_ };
    std::size_t const ny{ Ny_ };
    std::size_t const nz{ Nz_ };

    std::size_t const x_hi_base{ nx - 1 - t };

    // x-faces
    #pragma omp parallel for collapse( 2 )
    for ( std::size_t z = 0; z < nz - 1; ++z ) {
        for ( std::size_t y = 0; y < ny - 1; ++y ) {
            std::size_t const yz_grid_base{ y * Sy + z * Sz };
            std::size_t const yz_psi_base{ t * ( y + ny * z ) };

            for ( std::size_t d = 0; d < t; ++d ) {
                // Low face: x = d
                std::size_t const gi_lo{ yz_grid_base + d };
                std::size_t const pi_lo{ yz_psi_base + d };

                double const dEy_lo{ ( Ey[gi_lo + Sx] - Ey[gi_lo] ) * inv_dx };
                double const dEz_lo{ ( Ez[gi_lo + Sx] - Ez[gi_lo] ) * inv_dx };

                double const b_lo{ bBx[d] };
                double const c_lo{ cBx[d] };

                pByx[pi_lo] = b_lo * pByx[pi_lo] + c_lo * dEy_lo;
                pBzx[pi_lo] = b_lo * pBzx[pi_lo] + c_lo * dEz_lo;

                Bz[gi_lo] -= dt * pByx[pi_lo];
                By[gi_lo] += dt * pBzx[pi_lo];

                // High face: x = x_hi_base + d
                std::size_t const gi_hi{ yz_grid_base + x_hi_base + d };
                std::size_t const pi_hi{ face_x + pi_lo };

                double const dEy_hi{ ( Ey[gi_hi + Sx] - Ey[gi_hi] ) * inv_dx };
                double const dEz_hi{ ( Ez[gi_hi + Sx] - Ez[gi_hi] ) * inv_dx };

                double const b_hi{ bBx[t - 1 - d] };
                double const c_hi{ cBx[t - 1 - d] };

                pByx[pi_hi] = b_hi * pByx[pi_hi] + c_hi * dEy_hi;
                pBzx[pi_hi] = b_hi * pBzx[pi_hi] + c_hi * dEz_hi;

                Bz[gi_hi] -= dt * pByx[pi_hi];
                By[gi_hi] += dt * pBzx[pi_hi];
            }
        }
    }

    // y-faces
    std::size_t const y_hi_base{ ny - 1 - t };

    #pragma omp parallel for collapse( 2 )
    for ( std::size_t z = 0; z < nz - 1; ++z ) {
        for ( std::size_t x = 0; x < nx - 1; ++x ) {
            std::size_t const xz_grid_base{ x + z * Sz };
            std::size_t const xz_psi_base{ x + nx * t * z };

            for ( std::size_t d = 0; d < t; ++d ) {
                // Low face: y = d
                std::size_t const gi_lo{ xz_grid_base + d * Sy };
                std::size_t const pi_lo{ xz_psi_base + nx * d };

                double const dEx_lo{ ( Ex[gi_lo + Sy] - Ex[gi_lo] ) * inv_dy };
                double const dEz_lo{ ( Ez[gi_lo + Sy] - Ez[gi_lo] ) * inv_dy };

                double const b_lo{ bBy[d] };
                double const c_lo{ cBy[d] };

                pBxy[pi_lo] = b_lo * pBxy[pi_lo] + c_lo * dEx_lo;
                pBzy[pi_lo] = b_lo * pBzy[pi_lo] + c_lo * dEz_lo;

                Bz[gi_lo] += dt * pBxy[pi_lo];
                Bx[gi_lo] -= dt * pBzy[pi_lo];

                // High face: y = y_hi_base + d
                std::size_t const gi_hi{ xz_grid_base + ( y_hi_base + d ) * Sy };
                std::size_t const pi_hi{ face_y + pi_lo };

                double const dEx_hi{ ( Ex[gi_hi + Sy] - Ex[gi_hi] ) * inv_dy };
                double const dEz_hi{ ( Ez[gi_hi + Sy] - Ez[gi_hi] ) * inv_dy };

                double const b_hi{ bBy[t - 1 - d] };
                double const c_hi{ cBy[t - 1 - d] };

                pBxy[pi_hi] = b_hi * pBxy[pi_hi] + c_hi * dEx_hi;
                pBzy[pi_hi] = b_hi * pBzy[pi_hi] + c_hi * dEz_hi;

                Bz[gi_hi] += dt * pBxy[pi_hi];
                Bx[gi_hi] -= dt * pBzy[pi_hi];
            }
        }
    }

    // z-faces
    std::size_t const z_hi_base{ nz - 1 - t };

    #pragma omp parallel for collapse( 2 )
    for ( std::size_t y = 0; y < ny - 1; ++y ) {
        for ( std::size_t x = 0; x < nx - 1; ++x ) {
            std::size_t const xy_grid_base{ x + y * Sy };
            std::size_t const xy_psi_base{ x + nx * y };

            for ( std::size_t d = 0; d < t; ++d ) {
                // Low face: z = d
                std::size_t const gi_lo{ xy_grid_base + d * Sz };
                std::size_t const pi_lo{ xy_psi_base + nx * ny * d };

                double const dEx_lo{ ( Ex[gi_lo + Sz] - Ex[gi_lo] ) * inv_dz };
                double const dEy_lo{ ( Ey[gi_lo + Sz] - Ey[gi_lo] ) * inv_dz };

                double const b_lo{ bBz[d] };
                double const c_lo{ cBz[d] };

                pBxz[pi_lo] = b_lo * pBxz[pi_lo] + c_lo * dEx_lo;
                pByz[pi_lo] = b_lo * pByz[pi_lo] + c_lo * dEy_lo;

                Bx[gi_lo] += dt * pByz[pi_lo];
                By[gi_lo] -= dt * pBxz[pi_lo];

                // High face: z = z_hi_base + d
                std::size_t const gi_hi{ xy_grid_base + ( z_hi_base + d ) * Sz };
                std::size_t const pi_hi{ face_z + pi_lo };

                double const dEx_hi{ ( Ex[gi_hi + Sz] - Ex[gi_hi] ) * inv_dz };
                double const dEy_hi{ ( Ey[gi_hi + Sz] - Ey[gi_hi] ) * inv_dz };

                double const b_hi{ bBz[t - 1 - d] };
                double const c_hi{ cBz[t - 1 - d] };

                pBxz[pi_hi] = b_hi * pBxz[pi_hi] + c_hi * dEx_hi;
                pByz[pi_hi] = b_hi * pByz[pi_hi] + c_hi * dEy_hi;

                Bx[gi_hi] += dt * pByz[pi_hi];
                By[gi_hi] -= dt * pBxz[pi_hi];
            }
        }
    }
}

// E-field PML

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

    std::size_t const Sx{ 1 };
    std::size_t const Sy{ Nx_padded_ };
    std::size_t const Sz{ Nx_padded_ * Ny_padded_ };

    double const inv_dx{ 1.0 / dx };
    double const inv_dy{ 1.0 / dy };
    double const inv_dz{ 1.0 / dz };
    double const dt_csq{ dt * c_sq };

    double const* RESTRICT bEx{ b_Ex_ptr() };
    double const* RESTRICT cEx{ c_Ex_ptr() };
    double const* RESTRICT bEy{ b_Ey_ptr() };
    double const* RESTRICT cEy{ c_Ey_ptr() };
    double const* RESTRICT bEz{ b_Ez_ptr() };
    double const* RESTRICT cEz{ c_Ez_ptr() };

    double* RESTRICT pEyx{ psi_Eyx_ptr() };
    double* RESTRICT pEzx{ psi_Ezx_ptr() };
    double* RESTRICT pExy{ psi_Exy_ptr() };
    double* RESTRICT pEzy{ psi_Ezy_ptr() };
    double* RESTRICT pExz{ psi_Exz_ptr() };
    double* RESTRICT pEyz{ psi_Eyz_ptr() };

    std::size_t const nx{ Nx_ };
    std::size_t const ny{ Ny_ };
    std::size_t const nz{ Nz_ };

    std::size_t const d_max_lo_x{ std::min( t, nx - 2 ) };
    std::size_t const d_max_hi_x{ ( t >= 2 ) ? t - 1 : std::size_t{0} };
    std::size_t const x_hi_base{ nx - t };

    std::size_t const d_max_lo_y{ std::min( t, ny - 2 ) };
    std::size_t const d_max_hi_y{ ( t >= 2 ) ? t - 1 : std::size_t{0} };
    std::size_t const y_hi_base{ ny - t };

    std::size_t const d_max_lo_z{ std::min( t, nz - 2 ) };
    std::size_t const d_max_hi_z{ ( t >= 2 ) ? t - 1 : std::size_t{0} };
    std::size_t const z_hi_base{ nz - t };

    // x-faces
    #pragma omp parallel for collapse( 2 )
    for ( std::size_t z = 1; z < nz - 1; ++z ) {
        for ( std::size_t y = 1; y < ny - 1; ++y ) {
            std::size_t const yz_grid_base{ y * Sy + z * Sz };
            std::size_t const yz_psi_base{ t * ( y + ny * z ) };

            for ( std::size_t d = 0; d < d_max_lo_x; ++d ) {
                std::size_t const gi{ yz_grid_base + d + 1 };
                std::size_t const pi{ yz_psi_base + d };

                double const dBy{ ( By[gi] - By[gi - Sx] ) * inv_dx };
                double const dBz{ ( Bz[gi] - Bz[gi - Sx] ) * inv_dx };

                pEyx[pi] = bEx[d] * pEyx[pi] + cEx[d] * dBy;
                pEzx[pi] = bEx[d] * pEzx[pi] + cEx[d] * dBz;

                Ez[gi] += dt_csq * pEyx[pi];
                Ey[gi] -= dt_csq * pEzx[pi];
            }

            for ( std::size_t d = 0; d < d_max_hi_x; ++d ) {
                std::size_t const gi{ yz_grid_base + x_hi_base + d };
                std::size_t const pi{ face_x + yz_psi_base + d };

                double const dBy{ ( By[gi] - By[gi - Sx] ) * inv_dx };
                double const dBz{ ( Bz[gi] - Bz[gi - Sx] ) * inv_dx };

                std::size_t const rd{ t - 1 - d };

                pEyx[pi] = bEx[rd] * pEyx[pi] + cEx[rd] * dBy;
                pEzx[pi] = bEx[rd] * pEzx[pi] + cEx[rd] * dBz;

                Ez[gi] += dt_csq * pEyx[pi];
                Ey[gi] -= dt_csq * pEzx[pi];
            }
        }
    }

    // y-faces
    #pragma omp parallel for collapse( 2 )
    for ( std::size_t z = 1; z < nz - 1; ++z ) {
        for ( std::size_t x = 1; x < nx - 1; ++x ) {
            std::size_t const xz_grid_base{ x + z * Sz };
            std::size_t const xz_psi_base{ x + nx * t * z };

            for ( std::size_t d = 0; d < d_max_lo_y; ++d ) {
                std::size_t const gi{ xz_grid_base + ( d + 1 ) * Sy };
                std::size_t const pi{ xz_psi_base + nx * d };

                double const dBx{ ( Bx[gi] - Bx[gi - Sy] ) * inv_dy };
                double const dBz{ ( Bz[gi] - Bz[gi - Sy] ) * inv_dy };

                pExy[pi] = bEy[d] * pExy[pi] + cEy[d] * dBx;
                pEzy[pi] = bEy[d] * pEzy[pi] + cEy[d] * dBz;

                Ez[gi] -= dt_csq * pExy[pi];
                Ex[gi] += dt_csq * pEzy[pi];
            }

            for ( std::size_t d = 0; d < d_max_hi_y; ++d ) {
                std::size_t const gi{ xz_grid_base + ( y_hi_base + d ) * Sy };
                std::size_t const pi{ face_y + xz_psi_base + nx * d };

                double const dBx{ ( Bx[gi] - Bx[gi - Sy] ) * inv_dy };
                double const dBz{ ( Bz[gi] - Bz[gi - Sy] ) * inv_dy };

                std::size_t const rd{ t - 1 - d };

                pExy[pi] = bEy[rd] * pExy[pi] + cEy[rd] * dBx;
                pEzy[pi] = bEy[rd] * pEzy[pi] + cEy[rd] * dBz;

                Ez[gi] -= dt_csq * pExy[pi];
                Ex[gi] += dt_csq * pEzy[pi];
            }
        }
    }

    // z-faces
    #pragma omp parallel for collapse( 2 )
    for ( std::size_t y = 1; y < ny - 1; ++y ) {
        for ( std::size_t x = 1; x < nx - 1; ++x ) {
            std::size_t const xy_grid_base{ x + y * Sy };
            std::size_t const xy_psi_base{ x + nx * y };

            for ( std::size_t d = 0; d < d_max_lo_z; ++d ) {
                std::size_t const gi{ xy_grid_base + ( d + 1 ) * Sz };
                std::size_t const pi{ xy_psi_base + nx * ny * d };

                double const dBx{ ( Bx[gi] - Bx[gi - Sz] ) * inv_dz };
                double const dBy{ ( By[gi] - By[gi - Sz] ) * inv_dz };

                pExz[pi] = bEz[d] * pExz[pi] + cEz[d] * dBx;
                pEyz[pi] = bEz[d] * pEyz[pi] + cEz[d] * dBy;

                Ey[gi] += dt_csq * pExz[pi];
                Ex[gi] -= dt_csq * pEyz[pi];
            }

            for ( std::size_t d = 0; d < d_max_hi_z; ++d ) {
                std::size_t const gi{ xy_grid_base + ( z_hi_base + d ) * Sz };
                std::size_t const pi{ face_z + xy_psi_base + nx * ny * d };

                double const dBx{ ( Bx[gi] - Bx[gi - Sz] ) * inv_dz };
                double const dBy{ ( By[gi] - By[gi - Sz] ) * inv_dz };

                std::size_t const rd{ t - 1 - d };

                pExz[pi] = bEz[rd] * pExz[pi] + cEz[rd] * dBx;
                pEyz[pi] = bEz[rd] * pEyz[pi] + cEz[rd] * dBy;

                Ey[gi] += dt_csq * pExz[pi];
                Ex[gi] -= dt_csq * pEyz[pi];
            }
        }
    }
}