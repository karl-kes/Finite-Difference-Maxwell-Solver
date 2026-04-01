#include "pml.hpp"

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
, sigma_max_x_{ config.pml_sigma_max_x }
, sigma_max_y_{ config.pml_sigma_max_y }
, sigma_max_z_{ config.pml_sigma_max_z }
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
    double* RESTRICT b_Hx{ b_Hx_ptr() }; double* RESTRICT c_Hx{ c_Hx_ptr() }; double* RESTRICT kappa_Hx{ kappa_Hx_ptr() };
    double* RESTRICT b_Ey{ b_Ey_ptr() }; double* RESTRICT c_Ey{ c_Ey_ptr() }; double* RESTRICT kappa_Ey{ kappa_Ey_ptr() };
    double* RESTRICT b_Hy{ b_Hy_ptr() }; double* RESTRICT c_Hy{ c_Hy_ptr() }; double* RESTRICT kappa_Hy{ kappa_Hy_ptr() };
    double* RESTRICT b_Ez{ b_Ez_ptr() }; double* RESTRICT c_Ez{ c_Ez_ptr() }; double* RESTRICT kappa_Ez{ kappa_Ez_ptr() };
    double* RESTRICT b_Hz{ b_Hz_ptr() }; double* RESTRICT c_Hz{ c_Hz_ptr() }; double* RESTRICT kappa_Hz{ kappa_Hz_ptr() };

    for ( std::size_t i{ 0 }; i < thickness_; ++i ) {
        double const depth_E{ ( d - static_cast<double>( i ) ) / d };
        double const kap_E{ kappa( depth_E ) };
        double const alp_E{ alpha( depth_E ) };

        double const depth_H{ ( d - ( static_cast<double>( i ) + 0.5 ) ) / d };
        double const depth_H_clamped{ std::max( depth_H, 0.0 ) };
        double const kap_H{ kappa( depth_H_clamped ) };
        double const alp_H{ alpha( depth_H_clamped ) };

        // x-axis coefficients:
        compute_coefficients( sigma_x(depth_E), kap_E, alp_E, config.dt, config.eps, b_Ex[i], c_Ex[i] );
        compute_coefficients( sigma_x(depth_H_clamped), kap_H, alp_H, config.dt, config.eps, b_Hx[i], c_Hx[i] );
        kappa_Ex[i] = kap_E;
        kappa_Hx[i] = kap_H;

        // y-axis coefficients:
        compute_coefficients( sigma_y(depth_E), kap_E, alp_E, config.dt, config.eps, b_Ey[i], c_Ey[i] );
        compute_coefficients( sigma_y(depth_H_clamped), kap_H, alp_H, config.dt, config.eps, b_Hy[i], c_Hy[i] );
        kappa_Ey[i] = kap_E;
        kappa_Hy[i] = kap_H;

        // z-axis coefficients:
        compute_coefficients( sigma_z(depth_E), kap_E, alp_E, config.dt, config.eps, b_Ez[i], c_Ez[i] );
        compute_coefficients( sigma_z(depth_H_clamped), kap_H, alp_H, config.dt, config.eps, b_Hz[i], c_Hz[i] );
        kappa_Ez[i] = kap_E;
        kappa_Hz[i] = kap_H;
    }
}

void PML::update_H_psi(
    double* RESTRICT Ex, double* RESTRICT Ey, double* RESTRICT Ez,
    double* RESTRICT Hx, double* RESTRICT Hy, double* RESTRICT Hz,
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

    double const* RESTRICT bHx{ b_Hx_ptr() };
    double const* RESTRICT cHx{ c_Hx_ptr() };
    double const* RESTRICT bHy{ b_Hy_ptr() };
    double const* RESTRICT cHy{ c_Hy_ptr() };
    double const* RESTRICT bHz{ b_Hz_ptr() };
    double const* RESTRICT cHz{ c_Hz_ptr() };

    double* RESTRICT pHyx{ psi_Hyx_ptr() };
    double* RESTRICT pHzx{ psi_Hzx_ptr() };
    double* RESTRICT pHxy{ psi_Hxy_ptr() };
    double* RESTRICT pHzy{ psi_Hzy_ptr() };
    double* RESTRICT pHxz{ psi_Hxz_ptr() };
    double* RESTRICT pHyz{ psi_Hyz_ptr() };

    ASSUME_ALIGNED(bHx, SIMD_BYTES);
    ASSUME_ALIGNED(cHx, SIMD_BYTES);
    ASSUME_ALIGNED(bHy, SIMD_BYTES);
    ASSUME_ALIGNED(cHy, SIMD_BYTES);
    ASSUME_ALIGNED(bHz, SIMD_BYTES);
    ASSUME_ALIGNED(cHz, SIMD_BYTES);

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
                std::size_t const gi_lo{ yz_grid_base + d };
                std::size_t const pi_lo{ yz_psi_base + d };

                double const dEy_lo{ ( Ey[gi_lo + Sx] - Ey[gi_lo] ) * inv_dx };
                double const dEz_lo{ ( Ez[gi_lo + Sx] - Ez[gi_lo] ) * inv_dx };

                double const b_lo{ bHx[d] };
                double const c_lo{ cHx[d] };

                pHyx[pi_lo] = b_lo * pHyx[pi_lo] + c_lo * dEy_lo;
                pHzx[pi_lo] = b_lo * pHzx[pi_lo] + c_lo * dEz_lo;

                Hz[gi_lo] -= dt * pHyx[pi_lo];
                Hy[gi_lo] += dt * pHzx[pi_lo];

                std::size_t const gi_hi{ yz_grid_base + x_hi_base + d };
                std::size_t const pi_hi{ face_x + pi_lo };

                double const dEy_hi{ ( Ey[gi_hi + Sx] - Ey[gi_hi] ) * inv_dx };
                double const dEz_hi{ ( Ez[gi_hi + Sx] - Ez[gi_hi] ) * inv_dx };

                double const b_hi{ bHx[t - 1 - d] };
                double const c_hi{ cHx[t - 1 - d] };

                pHyx[pi_hi] = b_hi * pHyx[pi_hi] + c_hi * dEy_hi;
                pHzx[pi_hi] = b_hi * pHzx[pi_hi] + c_hi * dEz_hi;

                Hz[gi_hi] -= dt * pHyx[pi_hi];
                Hy[gi_hi] += dt * pHzx[pi_hi];
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
                std::size_t const gi_lo{ xz_grid_base + d * Sy };
                std::size_t const pi_lo{ xz_psi_base + nx * d };

                double const dEx_lo{ ( Ex[gi_lo + Sy] - Ex[gi_lo] ) * inv_dy };
                double const dEz_lo{ ( Ez[gi_lo + Sy] - Ez[gi_lo] ) * inv_dy };

                double const b_lo{ bHy[d] };
                double const c_lo{ cHy[d] };

                pHxy[pi_lo] = b_lo * pHxy[pi_lo] + c_lo * dEx_lo;
                pHzy[pi_lo] = b_lo * pHzy[pi_lo] + c_lo * dEz_lo;

                Hz[gi_lo] += dt * pHxy[pi_lo];
                Hx[gi_lo] -= dt * pHzy[pi_lo];

                std::size_t const gi_hi{ xz_grid_base + ( y_hi_base + d ) * Sy };
                std::size_t const pi_hi{ face_y + pi_lo };

                double const dEx_hi{ ( Ex[gi_hi + Sy] - Ex[gi_hi] ) * inv_dy };
                double const dEz_hi{ ( Ez[gi_hi + Sy] - Ez[gi_hi] ) * inv_dy };

                double const b_hi{ bHy[t - 1 - d] };
                double const c_hi{ cHy[t - 1 - d] };

                pHxy[pi_hi] = b_hi * pHxy[pi_hi] + c_hi * dEx_hi;
                pHzy[pi_hi] = b_hi * pHzy[pi_hi] + c_hi * dEz_hi;

                Hz[gi_hi] += dt * pHxy[pi_hi];
                Hx[gi_hi] -= dt * pHzy[pi_hi];
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
                std::size_t const gi_lo{ xy_grid_base + d * Sz };
                std::size_t const pi_lo{ xy_psi_base + nx * ny * d };

                double const dEx_lo{ ( Ex[gi_lo + Sz] - Ex[gi_lo] ) * inv_dz };
                double const dEy_lo{ ( Ey[gi_lo + Sz] - Ey[gi_lo] ) * inv_dz };

                double const b_lo{ bHz[d] };
                double const c_lo{ cHz[d] };

                pHxz[pi_lo] = b_lo * pHxz[pi_lo] + c_lo * dEx_lo;
                pHyz[pi_lo] = b_lo * pHyz[pi_lo] + c_lo * dEy_lo;

                Hx[gi_lo] += dt * pHyz[pi_lo];
                Hy[gi_lo] -= dt * pHxz[pi_lo];

                std::size_t const gi_hi{ xy_grid_base + ( z_hi_base + d ) * Sz };
                std::size_t const pi_hi{ face_z + pi_lo };

                double const dEx_hi{ ( Ex[gi_hi + Sz] - Ex[gi_hi] ) * inv_dz };
                double const dEy_hi{ ( Ey[gi_hi + Sz] - Ey[gi_hi] ) * inv_dz };

                double const b_hi{ bHz[t - 1 - d] };
                double const c_hi{ cHz[t - 1 - d] };

                pHxz[pi_hi] = b_hi * pHxz[pi_hi] + c_hi * dEx_hi;
                pHyz[pi_hi] = b_hi * pHyz[pi_hi] + c_hi * dEy_hi;

                Hx[gi_hi] += dt * pHyz[pi_hi];
                Hy[gi_hi] -= dt * pHxz[pi_hi];
            }
        }
    }
}

void PML::update_E_psi(
    double* RESTRICT Ex, double* RESTRICT Ey, double* RESTRICT Ez,
    double* RESTRICT Hx, double* RESTRICT Hy, double* RESTRICT Hz,
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

    ASSUME_ALIGNED(bEx, SIMD_BYTES);
    ASSUME_ALIGNED(cEx, SIMD_BYTES);
    ASSUME_ALIGNED(bEy, SIMD_BYTES);
    ASSUME_ALIGNED(cEy, SIMD_BYTES);
    ASSUME_ALIGNED(bEz, SIMD_BYTES);
    ASSUME_ALIGNED(cEz, SIMD_BYTES);

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

                double const dHy{ ( Hy[gi] - Hy[gi - Sx] ) * inv_dx };
                double const dHz{ ( Hz[gi] - Hz[gi - Sx] ) * inv_dx };

                pEyx[pi] = bEx[d] * pEyx[pi] + cEx[d] * dHy;
                pEzx[pi] = bEx[d] * pEzx[pi] + cEx[d] * dHz;

                Ez[gi] += dt * pEyx[pi];
                Ey[gi] -= dt * pEzx[pi];
            }

            for ( std::size_t d = 0; d < d_max_hi_x; ++d ) {
                std::size_t const gi{ yz_grid_base + x_hi_base + d };
                std::size_t const pi{ face_x + yz_psi_base + d };

                double const dHy{ ( Hy[gi] - Hy[gi - Sx] ) * inv_dx };
                double const dHz{ ( Hz[gi] - Hz[gi - Sx] ) * inv_dx };

                std::size_t const rd{ t - 1 - d };

                pEyx[pi] = bEx[rd] * pEyx[pi] + cEx[rd] * dHy;
                pEzx[pi] = bEx[rd] * pEzx[pi] + cEx[rd] * dHz;

                Ez[gi] += dt * pEyx[pi];
                Ey[gi] -= dt * pEzx[pi];
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

                double const dHx{ ( Hx[gi] - Hx[gi - Sy] ) * inv_dy };
                double const dHz{ ( Hz[gi] - Hz[gi - Sy] ) * inv_dy };

                pExy[pi] = bEy[d] * pExy[pi] + cEy[d] * dHx;
                pEzy[pi] = bEy[d] * pEzy[pi] + cEy[d] * dHz;

                Ez[gi] -= dt * pExy[pi];
                Ex[gi] += dt * pEzy[pi];
            }

            for ( std::size_t d = 0; d < d_max_hi_y; ++d ) {
                std::size_t const gi{ xz_grid_base + ( y_hi_base + d ) * Sy };
                std::size_t const pi{ face_y + xz_psi_base + nx * d };

                double const dHx{ ( Hx[gi] - Hx[gi - Sy] ) * inv_dy };
                double const dHz{ ( Hz[gi] - Hz[gi - Sy] ) * inv_dy };

                std::size_t const rd{ t - 1 - d };

                pExy[pi] = bEy[rd] * pExy[pi] + cEy[rd] * dHx;
                pEzy[pi] = bEy[rd] * pEzy[pi] + cEy[rd] * dHz;

                Ez[gi] -= dt * pExy[pi];
                Ex[gi] += dt * pEzy[pi];
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

                double const dHx{ ( Hx[gi] - Hx[gi - Sz] ) * inv_dz };
                double const dHy{ ( Hy[gi] - Hy[gi - Sz] ) * inv_dz };

                pExz[pi] = bEz[d] * pExz[pi] + cEz[d] * dHx;
                pEyz[pi] = bEz[d] * pEyz[pi] + cEz[d] * dHy;

                Ey[gi] += dt * pExz[pi];
                Ex[gi] -= dt * pEyz[pi];
            }

            for ( std::size_t d = 0; d < d_max_hi_z; ++d ) {
                std::size_t const gi{ xy_grid_base + ( z_hi_base + d ) * Sz };
                std::size_t const pi{ face_z + xy_psi_base + nx * ny * d };

                double const dHx{ ( Hx[gi] - Hx[gi - Sz] ) * inv_dz };
                double const dHy{ ( Hy[gi] - Hy[gi - Sz] ) * inv_dz };

                std::size_t const rd{ t - 1 - d };

                pExz[pi] = bEz[rd] * pExz[pi] + cEz[rd] * dHx;
                pEyz[pi] = bEz[rd] * pEyz[pi] + cEz[rd] * dHy;

                Ey[gi] += dt * pExz[pi];
                Ex[gi] -= dt * pEyz[pi];
            }
        }
    }
}