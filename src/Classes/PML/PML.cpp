#include "PML.hpp"

PML::PML( Simulation_Config const &config )
: thickness_{ config.use_pml ? config.pml_thickness : 0 }
, Nx_{ config.Nx + 1 }
, Ny_{ config.Ny + 1 }
, Nz_{ config.Nz + 1 }
, order_{ config.pml_order }
, sigma_max_{ config.pml_sigma_max }
, kappa_max_{ config.pml_kappa_max }
, alpha_max_{ config.pml_alpha_max } {
    if ( !config.use_pml ) return;

    auto alloc_coeff = [&]() { return std::make_unique<double[]>( thickness_ ); };

    b_Ex_ = alloc_coeff();  c_Ex_ = alloc_coeff();  kappa_Ex_ = alloc_coeff();
    b_Bx_ = alloc_coeff();  c_Bx_ = alloc_coeff();  kappa_Bx_ = alloc_coeff();
    b_Ey_ = alloc_coeff();  c_Ey_ = alloc_coeff();  kappa_Ey_ = alloc_coeff();
    b_By_ = alloc_coeff();  c_By_ = alloc_coeff();  kappa_By_ = alloc_coeff();
    b_Ez_ = alloc_coeff();  c_Ez_ = alloc_coeff();  kappa_Ez_ = alloc_coeff();
    b_Bz_ = alloc_coeff();  c_Bz_ = alloc_coeff();  kappa_Bz_ = alloc_coeff();

    // Compute grading coefficients for each PML layer:
    // E-fields are at integer positions, B-fields at half-integer.
    double const d{ static_cast<double>( thickness_ ) };

    for ( std::size_t i = 0; i < thickness_; ++i ) {
        // E-field: depth at integer position (cell edge):
        double const depth_E{ ( d - static_cast<double>( i ) ) / d };
        double const sig_E{ sigma( depth_E ) };
        double const kap_E{ kappa( depth_E ) };
        double const alp_E{ alpha( depth_E ) };

        // B-field: depth at half-integer position (cell center):
        double const depth_B{ ( d - ( static_cast<double>( i ) + 0.5 ) ) / d };
        double const sig_B{ sigma( depth_B ) };
        double const kap_B{ kappa( depth_B ) };
        double const alp_B{ alpha( depth_B ) };

        // All directions use the same grading (isotropic grid assumed):
        compute_coefficients( sig_E, kap_E, alp_E, config.dt, config.eps,
                              b_Ex_[i], c_Ex_[i] );
        compute_coefficients( sig_B, kap_B, alp_B, config.dt, config.eps,
                              b_Bx_[i], c_Bx_[i] );

        kappa_Ex_[i] = kap_E;
        kappa_Bx_[i] = kap_B;

        compute_coefficients( sig_E, kap_E, alp_E, config.dt, config.eps,
                              b_Ey_[i], c_Ey_[i] );
        compute_coefficients( sig_B, kap_B, alp_B, config.dt, config.eps,
                              b_By_[i], c_By_[i] );

        kappa_Ey_[i] = kap_E;
        kappa_By_[i] = kap_B;

        compute_coefficients( sig_E, kap_E, alp_E, config.dt, config.eps,
                              b_Ez_[i], c_Ez_[i] );
        compute_coefficients( sig_B, kap_B, alp_B, config.dt, config.eps,
                              b_Bz_[i], c_Bz_[i] );

        kappa_Ez_[i] = kap_E;
        kappa_Bz_[i] = kap_B;
    }

    // Allocate psi arrays (two faces per direction, stored contiguously):
    auto alloc_psi_x = [&]() { return std::make_unique<double[]>( 2 * psi_size_x() ); };
    auto alloc_psi_y = [&]() { return std::make_unique<double[]>( 2 * psi_size_y() ); };
    auto alloc_psi_z = [&]() { return std::make_unique<double[]>( 2 * psi_size_z() ); };

    psi_Eyx_ = alloc_psi_x();  psi_Ezx_ = alloc_psi_x();
    psi_Exy_ = alloc_psi_y();  psi_Ezy_ = alloc_psi_y();
    psi_Exz_ = alloc_psi_z();  psi_Eyz_ = alloc_psi_z();

    psi_Byx_ = alloc_psi_x();  psi_Bzx_ = alloc_psi_x();
    psi_Bxy_ = alloc_psi_y();  psi_Bzy_ = alloc_psi_y();
    psi_Bxz_ = alloc_psi_z();  psi_Byz_ = alloc_psi_z();
}

void PML::update_B_psi( double *Ex, double *Ey, double *Ez,
                         double *Bx, double *By, double *Bz,
                         double const dt, double const dx, double const dy, double const dz ) {
    if ( !is_active() ) return;

    std::size_t const t{ thickness_ };
    std::size_t const face_x{ psi_size_x() };
    std::size_t const face_y{ psi_size_y() };
    std::size_t const face_z{ psi_size_z() };

    // x-faces (low and high):
    // B update uses E derivatives; x-derivative terms need PML correction.
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

                psi_Byx_[pi_lo] = b_Bx_[d] * psi_Byx_[pi_lo] + c_Bx_[d] * dEy_dx_lo;
                psi_Bzx_[pi_lo] = b_Bx_[d] * psi_Bzx_[pi_lo] + c_Bx_[d] * dEz_dx_lo;

                Bz[gi_lo] -= dt * psi_Byx_[pi_lo];
                By[gi_lo] += dt * psi_Bzx_[pi_lo];

                // High x-face:
                std::size_t const x_hi{ Nx_ - 1 - t + d };
                std::size_t const gi_hi{ idx( x_hi, y, z ) };
                std::size_t const pi_hi{ face_x + psi_idx_x( d, y, z ) };

                double const dEy_dx_hi{ ( Ey[idx( x_hi + 1, y, z )] - Ey[gi_hi] ) / dx };
                double const dEz_dx_hi{ ( Ez[idx( x_hi + 1, y, z )] - Ez[gi_hi] ) / dx };

                psi_Byx_[pi_hi] = b_Bx_[t - 1 - d] * psi_Byx_[pi_hi] + c_Bx_[t - 1 - d] * dEy_dx_hi;
                psi_Bzx_[pi_hi] = b_Bx_[t - 1 - d] * psi_Bzx_[pi_hi] + c_Bx_[t - 1 - d] * dEz_dx_hi;

                Bz[gi_hi] -= dt * psi_Byx_[pi_hi];
                By[gi_hi] += dt * psi_Bzx_[pi_hi];
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

                psi_Bxy_[pi_lo] = b_By_[d] * psi_Bxy_[pi_lo] + c_By_[d] * dEx_dy_lo;
                psi_Bzy_[pi_lo] = b_By_[d] * psi_Bzy_[pi_lo] + c_By_[d] * dEz_dy_lo;

                Bz[gi_lo] += dt * psi_Bxy_[pi_lo];
                Bx[gi_lo] -= dt * psi_Bzy_[pi_lo];

                // High y-face:
                std::size_t const y_hi{ Ny_ - 1 - t + d };
                std::size_t const gi_hi{ idx( x, y_hi, z ) };
                std::size_t const pi_hi{ face_y + psi_idx_y( x, d, z ) };

                double const dEx_dy_hi{ ( Ex[idx( x, y_hi + 1, z )] - Ex[gi_hi] ) / dy };
                double const dEz_dy_hi{ ( Ez[idx( x, y_hi + 1, z )] - Ez[gi_hi] ) / dy };

                psi_Bxy_[pi_hi] = b_By_[t - 1 - d] * psi_Bxy_[pi_hi] + c_By_[t - 1 - d] * dEx_dy_hi;
                psi_Bzy_[pi_hi] = b_By_[t - 1 - d] * psi_Bzy_[pi_hi] + c_By_[t - 1 - d] * dEz_dy_hi;

                Bz[gi_hi] += dt * psi_Bxy_[pi_hi];
                Bx[gi_hi] -= dt * psi_Bzy_[pi_hi];
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

                psi_Bxz_[pi_lo] = b_Bz_[d] * psi_Bxz_[pi_lo] + c_Bz_[d] * dEx_dz_lo;
                psi_Byz_[pi_lo] = b_Bz_[d] * psi_Byz_[pi_lo] + c_Bz_[d] * dEy_dz_lo;

                Bx[gi_lo] += dt * psi_Byz_[pi_lo];
                By[gi_lo] -= dt * psi_Bxz_[pi_lo];

                // High z-face:
                std::size_t const z_hi{ Nz_ - 1 - t + d };
                std::size_t const gi_hi{ idx( x, y, z_hi ) };
                std::size_t const pi_hi{ face_z + psi_idx_z( x, y, d ) };

                double const dEx_dz_hi{ ( Ex[idx( x, y, z_hi + 1 )] - Ex[gi_hi] ) / dz };
                double const dEy_dz_hi{ ( Ey[idx( x, y, z_hi + 1 )] - Ey[gi_hi] ) / dz };

                psi_Bxz_[pi_hi] = b_Bz_[t - 1 - d] * psi_Bxz_[pi_hi] + c_Bz_[t - 1 - d] * dEx_dz_hi;
                psi_Byz_[pi_hi] = b_Bz_[t - 1 - d] * psi_Byz_[pi_hi] + c_Bz_[t - 1 - d] * dEy_dz_hi;

                Bx[gi_hi] += dt * psi_Byz_[pi_hi];
                By[gi_hi] -= dt * psi_Bxz_[pi_hi];
            }
        }
    }
}

void PML::update_E_psi( double *Ex, double *Ey, double *Ez,
                        double *Bx, double *By, double *Bz,
                        double const dt, double const dx, double const dy, double const dz,
                        double const c_sq ) {
    if ( !is_active() ) return;

    std::size_t const t{ thickness_ };
    std::size_t const face_x{ psi_size_x() };
    std::size_t const face_y{ psi_size_y() };
    std::size_t const face_z{ psi_size_z() };

    // x-faces (low and high):
    // E update uses B derivatives; x-derivative terms need PML correction.
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

                psi_Eyx_[pi_lo] = b_Ex_[d] * psi_Eyx_[pi_lo] + c_Ex_[d] * dBy_dx_lo;
                psi_Ezx_[pi_lo] = b_Ex_[d] * psi_Ezx_[pi_lo] + c_Ex_[d] * dBz_dx_lo;

                Ez[gi_lo] += dt * c_sq * psi_Eyx_[pi_lo];
                Ey[gi_lo] -= dt * c_sq * psi_Ezx_[pi_lo];

                // High x-face:
                std::size_t const x_hi{ Nx_ - t + d };
                if ( x_hi >= Nx_ ) continue;
                std::size_t const gi_hi{ idx( x_hi, y, z ) };
                std::size_t const pi_hi{ face_x + psi_idx_x( d, y, z ) };

                double const dBy_dx_hi{ ( By[gi_hi] - By[idx( x_hi - 1, y, z )] ) / dx };
                double const dBz_dx_hi{ ( Bz[gi_hi] - Bz[idx( x_hi - 1, y, z )] ) / dx };

                psi_Eyx_[pi_hi] = b_Ex_[t - 1 - d] * psi_Eyx_[pi_hi] + c_Ex_[t - 1 - d] * dBy_dx_hi;
                psi_Ezx_[pi_hi] = b_Ex_[t - 1 - d] * psi_Ezx_[pi_hi] + c_Ex_[t - 1 - d] * dBz_dx_hi;

                Ez[gi_hi] += dt * c_sq * psi_Eyx_[pi_hi];
                Ey[gi_hi] -= dt * c_sq * psi_Ezx_[pi_hi];
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

                psi_Exy_[pi_lo] = b_Ey_[d] * psi_Exy_[pi_lo] + c_Ey_[d] * dBx_dy_lo;
                psi_Ezy_[pi_lo] = b_Ey_[d] * psi_Ezy_[pi_lo] + c_Ey_[d] * dBz_dy_lo;

                Ez[gi_lo] -= dt * c_sq * psi_Exy_[pi_lo];
                Ex[gi_lo] += dt * c_sq * psi_Ezy_[pi_lo];

                // High y-face:
                std::size_t const y_hi{ Ny_ - t + d };
                if ( y_hi >= Ny_ ) continue;
                std::size_t const gi_hi{ idx( x, y_hi, z ) };
                std::size_t const pi_hi{ face_y + psi_idx_y( x, d, z ) };

                double const dBx_dy_hi{ ( Bx[gi_hi] - Bx[idx( x, y_hi - 1, z )] ) / dy };
                double const dBz_dy_hi{ ( Bz[gi_hi] - Bz[idx( x, y_hi - 1, z )] ) / dy };

                psi_Exy_[pi_hi] = b_Ey_[t - 1 - d] * psi_Exy_[pi_hi] + c_Ey_[t - 1 - d] * dBx_dy_hi;
                psi_Ezy_[pi_hi] = b_Ey_[t - 1 - d] * psi_Ezy_[pi_hi] + c_Ey_[t - 1 - d] * dBz_dy_hi;

                Ez[gi_hi] -= dt * c_sq * psi_Exy_[pi_hi];
                Ex[gi_hi] += dt * c_sq * psi_Ezy_[pi_hi];
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

                psi_Exz_[pi_lo] = b_Ez_[d] * psi_Exz_[pi_lo] + c_Ez_[d] * dBx_dz_lo;
                psi_Eyz_[pi_lo] = b_Ez_[d] * psi_Eyz_[pi_lo] + c_Ez_[d] * dBy_dz_lo;

                Ey[gi_lo] += dt * c_sq * psi_Exz_[pi_lo];
                Ex[gi_lo] -= dt * c_sq * psi_Eyz_[pi_lo];

                // High z-face:
                std::size_t const z_hi{ Nz_ - t + d };
                if ( z_hi >= Nz_ ) continue;
                std::size_t const gi_hi{ idx( x, y, z_hi ) };
                std::size_t const pi_hi{ face_z + psi_idx_z( x, y, d ) };

                double const dBx_dz_hi{ ( Bx[gi_hi] - Bx[idx( x, y, z_hi - 1 )] ) / dz };
                double const dBy_dz_hi{ ( By[gi_hi] - By[idx( x, y, z_hi - 1 )] ) / dz };

                psi_Exz_[pi_hi] = b_Ez_[t - 1 - d] * psi_Exz_[pi_hi] + c_Ez_[t - 1 - d] * dBx_dz_hi;
                psi_Eyz_[pi_hi] = b_Ez_[t - 1 - d] * psi_Eyz_[pi_hi] + c_Ez_[t - 1 - d] * dBy_dz_hi;

                Ey[gi_hi] += dt * c_sq * psi_Exz_[pi_hi];
                Ex[gi_hi] -= dt * c_sq * psi_Eyz_[pi_hi];
            }
        }
    }
}