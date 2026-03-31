#include "grid.hpp"

#include "../Source/source.hpp"

#include <cmath>
#include <stdexcept>
#include <omp.h>
#include <utility>

Grid::Grid( Simulation_Config const &config )
: Nx_{ config.Nx + 1 }
, Ny_{ config.Ny + 1 }
, Nz_{ config.Nz + 1 }
, Nx_padded_{ AlignedSoA<double>::round_up( Nx_ ) }
, Ny_padded_{ AlignedSoA<double>::round_up( Ny_ ) }
, Nz_padded_{ AlignedSoA<double>::round_up( Nz_ ) }
, dx_{ config.dx }
, dy_{ config.dy }
, dz_{ config.dz }
, dt_{ config.dt }
, data_{ Nx_padded_ * Ny_padded_ * Nz_padded_, NUM_ARRAYS_ }
, pml_{ config } {
    std::size_t const total{ Nx_padded_ * Ny_padded_ * Nz_padded_ };

    // Default initialize eps:
    std::fill_n( eps_x_ptr(), total, config.eps );
    std::fill_n( eps_y_ptr(), total, config.eps );
    std::fill_n( eps_z_ptr(), total, config.eps );

    // Default initialize mu:
    std::fill_n( mu_x_ptr(), total, config.mu );
    std::fill_n( mu_y_ptr(), total, config.mu );
    std::fill_n( mu_z_ptr(), total, config.mu );

    // sig arrays are already zero from AlignedSoA zero-init.

    // Bake hot-loop coeffs:
    #pragma omp simd
    for ( std::size_t i = 0; i < total; ++i ) {
        // Losses in all components:
        double const loss_x{ sig_x_ptr()[i] * dt_ / ( 2.0 * eps_x_ptr()[i] ) };
        double const loss_y{ sig_y_ptr()[i] * dt_ / ( 2.0 * eps_y_ptr()[i] ) };
        double const loss_z{ sig_z_ptr()[i] * dt_ / ( 2.0 * eps_z_ptr()[i] ) };

        // C_a constant:
        Ca_x_ptr()[i] = ( 1.0 - loss_x ) / ( 1.0 + loss_x );
        Ca_y_ptr()[i] = ( 1.0 - loss_y ) / ( 1.0 + loss_y );
        Ca_z_ptr()[i] = ( 1.0 - loss_z ) / ( 1.0 + loss_z );

        double const denom_x{ eps_x_ptr()[i] / dt_ + sig_x_ptr()[i] / 2.0 };
        double const denom_y{ eps_y_ptr()[i] / dt_ + sig_y_ptr()[i] / 2.0 };
        double const denom_z{ eps_z_ptr()[i] / dt_ + sig_z_ptr()[i] / 2.0 };
        
        // C_b constant:
        Cb_x_ptr()[i] = 1.0 / denom_x;
        Cb_y_ptr()[i] = 1.0 / denom_y;
        Cb_z_ptr()[i] = 1.0 / denom_z;

        // D_b constant:
        Db_x_ptr()[i] = dt_ / mu_x_ptr()[i];
        Db_y_ptr()[i] = dt_ / mu_y_ptr()[i];
        Db_z_ptr()[i] = dt_ / mu_z_ptr()[i];
    }
}

Grid::~Grid() = default;

void Grid::add_source( std::unique_ptr<Source> source ) {
    sources_.push_back( std::move( source ) );
}

void Grid::apply_sources( std::size_t const time_step ) {
    for ( auto const &source : sources_ ) {
        source->apply( *this, time_step );
    }
}

void Grid::update_H() {
    double* RESTRICT Hx{ Hx_ptr() };
    double* RESTRICT Hy{ Hy_ptr() };
    double* RESTRICT Hz{ Hz_ptr() };

    double* RESTRICT Ex{ Ex_ptr() };
    double* RESTRICT Ey{ Ey_ptr() };
    double* RESTRICT Ez{ Ez_ptr() };

    double* RESTRICT Db_x{ Db_x_ptr() };
    double* RESTRICT Db_y{ Db_y_ptr() };
    double* RESTRICT Db_z{ Db_z_ptr() };

    ASSUME_ALIGNED(Ex, SIMD_BYTES);
    ASSUME_ALIGNED(Ey, SIMD_BYTES);
    ASSUME_ALIGNED(Ez, SIMD_BYTES);

    ASSUME_ALIGNED(Hx, SIMD_BYTES);
    ASSUME_ALIGNED(Hy, SIMD_BYTES);
    ASSUME_ALIGNED(Hz, SIMD_BYTES);

    ASSUME_ALIGNED(Db_x, SIMD_BYTES);
    ASSUME_ALIGNED(Db_y, SIMD_BYTES);
    ASSUME_ALIGNED(Db_z, SIMD_BYTES);

    double const inv_dx{ 1.0 / dx_ };
    double const inv_dy{ 1.0 / dy_ };
    double const inv_dz{ 1.0 / dz_ };

    std::size_t const nx{ Nx_ };
    std::size_t const ny{ Ny_ };
    std::size_t const nz{ Nz_ };

    std::size_t const Sy{ Nx_padded_ };
    std::size_t const Sz{ Nx_padded_ * Ny_padded_ };

    #pragma omp parallel for collapse( 2 ) schedule( static )
    for ( std::size_t z = 0; z < nz - 1; ++z ) {
        for ( std::size_t y = 0; y < ny - 1; ++y ) {
            std::size_t const base{ y * Sy + z * Sz };

            #pragma omp simd
            for ( std::size_t x = 0; x < nx - 1; ++x ) {
                std::size_t const i{ base + x };

                double const curl_x{
                    ( Ez[i + Sy] - Ez[i] ) * inv_dy
                  - ( Ey[i + Sz] - Ey[i] ) * inv_dz
                };
                double const curl_y{
                    ( Ex[i + Sz] - Ex[i] ) * inv_dz
                  - ( Ez[i + 1]  - Ez[i] ) * inv_dx
                };
                double const curl_z{
                    ( Ey[i + 1]  - Ey[i] ) * inv_dx
                  - ( Ex[i + Sy] - Ex[i] ) * inv_dy
                };

                Hx[i] -= Db_x[i] * curl_x;
                Hy[i] -= Db_y[i] * curl_y;
                Hz[i] -= Db_z[i] * curl_z;
            }
        }
    }

    pml_.update_H_psi(
        Ex_ptr(), Ey_ptr(), Ez_ptr(),
        Hx_ptr(), Hy_ptr(), Hz_ptr(),
        dt_, dx_, dy_, dz_ 
    );    
}

void Grid::update_E() {
    double* RESTRICT Hx{ Hx_ptr() };
    double* RESTRICT Hy{ Hy_ptr() };
    double* RESTRICT Hz{ Hz_ptr() };

    double* RESTRICT Ex{ Ex_ptr() };
    double* RESTRICT Ey{ Ey_ptr() };
    double* RESTRICT Ez{ Ez_ptr() };

    double* RESTRICT Jx{ Jx_ptr() };
    double* RESTRICT Jy{ Jy_ptr() };
    double* RESTRICT Jz{ Jz_ptr() };

    double* RESTRICT Ca_x{ Ca_x_ptr() };
    double* RESTRICT Ca_y{ Ca_y_ptr() };
    double* RESTRICT Ca_z{ Ca_z_ptr() };

    double* RESTRICT Cb_x{ Cb_x_ptr() };
    double* RESTRICT Cb_y{ Cb_y_ptr() };
    double* RESTRICT Cb_z{ Cb_z_ptr() };

    ASSUME_ALIGNED(Ex, SIMD_BYTES);
    ASSUME_ALIGNED(Ey, SIMD_BYTES);
    ASSUME_ALIGNED(Ez, SIMD_BYTES);

    ASSUME_ALIGNED(Hx, SIMD_BYTES);
    ASSUME_ALIGNED(Hy, SIMD_BYTES);
    ASSUME_ALIGNED(Hz, SIMD_BYTES);

    ASSUME_ALIGNED(Jx, SIMD_BYTES);
    ASSUME_ALIGNED(Jy, SIMD_BYTES);
    ASSUME_ALIGNED(Jz, SIMD_BYTES);

    ASSUME_ALIGNED(Ca_x, SIMD_BYTES);
    ASSUME_ALIGNED(Ca_y, SIMD_BYTES);
    ASSUME_ALIGNED(Ca_z, SIMD_BYTES);

    ASSUME_ALIGNED(Cb_x, SIMD_BYTES);
    ASSUME_ALIGNED(Cb_y, SIMD_BYTES);
    ASSUME_ALIGNED(Cb_z, SIMD_BYTES);

    double const inv_dx{ 1.0 / dx_ };
    double const inv_dy{ 1.0 / dy_ };
    double const inv_dz{ 1.0 / dz_ };

    std::size_t const nx{ Nx_ };
    std::size_t const ny{ Ny_ };
    std::size_t const nz{ Nz_ };

    std::size_t const Sy{ Nx_padded_ };
    std::size_t const Sz{ Nx_padded_ * Ny_padded_ };

    #pragma omp parallel for collapse( 2 ) schedule( static )
    for ( std::size_t z = 1; z < nz - 1; ++z ) {
        for ( std::size_t y = 1; y < ny - 1; ++y ) {
            std::size_t const base{ y * Sy + z * Sz };

            #pragma omp simd
            for ( std::size_t x = 1; x < nx - 1; ++x ) {
                std::size_t const i{ base + x };

                double const curl_x{
                    ( Hz[i] - Hz[i - Sy] ) * inv_dy
                  - ( Hy[i] - Hy[i - Sz] ) * inv_dz
                };
                double const curl_y{
                    ( Hx[i] - Hx[i - Sz] ) * inv_dz
                  - ( Hz[i] - Hz[i - 1]  ) * inv_dx
                };
                double const curl_z{
                    ( Hy[i] - Hy[i - 1]  ) * inv_dx
                  - ( Hx[i] - Hx[i - Sy] ) * inv_dy
                };

                Ex[i] = Ca_x[i] * Ex[i] + Cb_x[i] * ( curl_x - Jx[i] );
                Ey[i] = Ca_y[i] * Ey[i] + Cb_y[i] * ( curl_y - Jy[i] );
                Ez[i] = Ca_z[i] * Ez[i] + Cb_z[i] * ( curl_z - Jz[i] );
            }
        }
    }

    pml_.update_E_psi(
        Ex_ptr(), Ey_ptr(), Ez_ptr(),
        Hx_ptr(), Hy_ptr(), Hz_ptr(),
        dt_, dx_, dy_, dz_
    );
}

void Grid::step() {
    update_H();
    update_E();
}

double Grid::field(
    Field const field, Component const component,
    std::size_t const x, std::size_t const y, std::size_t const z ) const {
    std::size_t const i{ idx(x,y,z) };
    if ( field == Field::ELECTRIC ) {
        switch ( component ) {
            case Component::X: return Ex_ptr()[i];
            case Component::Y: return Ey_ptr()[i];
            case Component::Z: return Ez_ptr()[i];
        }
    } else if ( field == Field::MAGNETIC ) {
        switch ( component ) {
            case Component::X: return Hx_ptr()[i];
            case Component::Y: return Hy_ptr()[i];
            case Component::Z: return Hz_ptr()[i];
        }
    }
    throw std::invalid_argument{ "Invalid field or component specifier" };
}

double &Grid::field(
    Field const field, Component const component,
    std::size_t const x, std::size_t const y, std::size_t const z ) {
    std::size_t const i{ idx(x,y,z) };
    if ( field == Field::ELECTRIC ) {
        switch ( component ) {
            case Component::X: return Ex_ptr()[i];
            case Component::Y: return Ey_ptr()[i];
            case Component::Z: return Ez_ptr()[i];
        }
    } else if ( field == Field::MAGNETIC ) {
        switch ( component ) {
            case Component::X: return Hx_ptr()[i];
            case Component::Y: return Hy_ptr()[i];
            case Component::Z: return Hz_ptr()[i];
        }
    }
    throw std::invalid_argument{ "Invalid field or component specifier" };
}

double Grid::field_magnitude(
    Field const field,
    std::size_t const x, std::size_t const y, std::size_t const z ) const {
    double Fx{ this->field( field, Component::X, x, y, z ) };
    double Fy{ this->field( field, Component::Y, x, y, z ) };
    double Fz{ this->field( field, Component::Z, x, y, z ) };

    return std::sqrt( Fx*Fx + Fy*Fy + Fz*Fz );
}

double Grid::e_energy() const {
    double energy{};
    double const dV{ dx_ * dy_ * dz_ };

    double const* RESTRICT Ex{ Ex_ptr() };
    double const* RESTRICT Ey{ Ey_ptr() };
    double const* RESTRICT Ez{ Ez_ptr() };

    double const* RESTRICT eps_x{ eps_x_ptr() };
    double const* RESTRICT eps_y{ eps_y_ptr() };
    double const* RESTRICT eps_z{ eps_z_ptr() };

    std::size_t const Sy{ Nx_padded_ };
    std::size_t const Sz{ Nx_padded_ * Ny_padded_ };

    #pragma omp parallel for collapse( 2 ) reduction( +:energy )
    for ( std::size_t z = 1; z < Nz_ - 1; ++z ) {
        for ( std::size_t y = 1; y < Ny_ - 1; ++y ) {
            std::size_t const base{ y * Sy + z * Sz };

            for ( std::size_t x = 1; x < Nx_ - 1; ++x ) {
                std::size_t const i{ base + x };

                double const Ex_avg{ ( Ex[i] + Ex[i+1] ) * 0.5 };
                double const Ey_avg{ ( Ey[i] + Ey[i+Sy] ) * 0.5 };
                double const Ez_avg{ ( Ez[i] + Ez[i+Sz] ) * 0.5 };

                energy += eps_x[i] * Ex_avg * Ex_avg
                        + eps_y[i] * Ey_avg * Ey_avg
                        + eps_z[i] * Ez_avg * Ez_avg;
            }
        }
    }
    return 0.5 * energy * dV;
}

double Grid::h_energy() const {
    double energy{};
    double const dV{ dx_ * dy_ * dz_ };

    double const* RESTRICT Hx{ Hx_ptr() };
    double const* RESTRICT Hy{ Hy_ptr() };
    double const* RESTRICT Hz{ Hz_ptr() };

    double const* RESTRICT mu_x{ mu_x_ptr() };
    double const* RESTRICT mu_y{ mu_y_ptr() };
    double const* RESTRICT mu_z{ mu_z_ptr() };

    std::size_t const Sy{ Nx_padded_ };
    std::size_t const Sz{ Nx_padded_ * Ny_padded_ };

    #pragma omp parallel for collapse( 2 ) reduction( +:energy )
    for ( std::size_t z = 1; z < Nz_ - 1; ++z ) {
        for ( std::size_t y = 1; y < Ny_ - 1; ++y ) {
            std::size_t const base{ y * Sy + z * Sz };

            for ( std::size_t x = 1; x < Nx_ - 1; ++x ) {
                std::size_t const i{ base + x };

                double const Hx_avg{ ( Hx[i] + Hx[i+Sy] + Hx[i+Sz] + Hx[i+Sy+Sz] ) * 0.25 };
                double const Hy_avg{ ( Hy[i] + Hy[i+1]  + Hy[i+Sz] + Hy[i+1+Sz]  ) * 0.25 };
                double const Hz_avg{ ( Hz[i] + Hz[i+1]  + Hz[i+Sy] + Hz[i+1+Sy]  ) * 0.25 };

                energy += mu_x[i] * Hx_avg * Hx_avg
                        + mu_y[i] * Hy_avg * Hy_avg
                        + mu_z[i] * Hz_avg * Hz_avg;
            }
        }
    }
    return 0.5 * energy * dV;
}

double Grid::total_energy() const {
    return e_energy() + h_energy();
}