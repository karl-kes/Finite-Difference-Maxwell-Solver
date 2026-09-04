#include "grid.hpp"

#include "../source/source.hpp"

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
    double const cfg_eps{ config.eps };
    double const cfg_mu{ config.mu };

    double* RESTRICT eps_x{ eps_x_ptr() }; ASSUME_ALIGNED(eps_x, SIMD_BYTES);
    double* RESTRICT eps_y{ eps_y_ptr() }; ASSUME_ALIGNED(eps_y, SIMD_BYTES);
    double* RESTRICT eps_z{ eps_z_ptr() }; ASSUME_ALIGNED(eps_z, SIMD_BYTES);

    double* RESTRICT mu_x{ mu_x_ptr() }; ASSUME_ALIGNED(mu_x, SIMD_BYTES);
    double* RESTRICT mu_y{ mu_y_ptr() }; ASSUME_ALIGNED(mu_y, SIMD_BYTES);
    double* RESTRICT mu_z{ mu_z_ptr() }; ASSUME_ALIGNED(mu_z, SIMD_BYTES);

    // Default initialize mu and eps:
    for ( std::size_t i = 0; i < total; ++i ) {
        eps_x[i] = cfg_eps;
        eps_y[i] = cfg_eps;
        eps_z[i] = cfg_eps;

        mu_x[i] = cfg_mu;
        mu_y[i] = cfg_mu;
        mu_z[i] = cfg_mu;
    }

    // sig arrays are already zero from AlignedSoA zero-init.

    // Adjust coefficients:
    bake_coefficients();
}

Grid::~Grid() = default;

void Grid::bake_coefficients() {
    std::size_t const total{ Nx_padded() * Ny_padded() * Nz_padded() };
    double const dt_local{ dt() };

    double const* RESTRICT eps_x{ eps_x_ptr() }; ASSUME_ALIGNED(eps_x, SIMD_BYTES);
    double const* RESTRICT eps_y{ eps_y_ptr() }; ASSUME_ALIGNED(eps_y, SIMD_BYTES);
    double const* RESTRICT eps_z{ eps_z_ptr() }; ASSUME_ALIGNED(eps_z, SIMD_BYTES);

    double const* RESTRICT sig_x{ sig_x_ptr() }; ASSUME_ALIGNED(sig_x, SIMD_BYTES);
    double const* RESTRICT sig_y{ sig_y_ptr() }; ASSUME_ALIGNED(sig_y, SIMD_BYTES);
    double const* RESTRICT sig_z{ sig_z_ptr() }; ASSUME_ALIGNED(sig_z, SIMD_BYTES);

    double const* RESTRICT mu_x{ mu_x_ptr() }; ASSUME_ALIGNED(mu_x, SIMD_BYTES);
    double const* RESTRICT mu_y{ mu_y_ptr() }; ASSUME_ALIGNED(mu_y, SIMD_BYTES);
    double const* RESTRICT mu_z{ mu_z_ptr() }; ASSUME_ALIGNED(mu_z, SIMD_BYTES);

    double* RESTRICT Ca_x{ Ca_x_ptr() }; ASSUME_ALIGNED(Ca_x, SIMD_BYTES);
    double* RESTRICT Ca_y{ Ca_y_ptr() }; ASSUME_ALIGNED(Ca_y, SIMD_BYTES);
    double* RESTRICT Ca_z{ Ca_z_ptr() }; ASSUME_ALIGNED(Ca_z, SIMD_BYTES);

    double* RESTRICT Cb_x{ Cb_x_ptr() }; ASSUME_ALIGNED(Cb_x, SIMD_BYTES);
    double* RESTRICT Cb_y{ Cb_y_ptr() }; ASSUME_ALIGNED(Cb_y, SIMD_BYTES);
    double* RESTRICT Cb_z{ Cb_z_ptr() }; ASSUME_ALIGNED(Cb_z, SIMD_BYTES);

    double* RESTRICT Db_x{ Db_x_ptr() }; ASSUME_ALIGNED(Db_x, SIMD_BYTES);
    double* RESTRICT Db_y{ Db_y_ptr() }; ASSUME_ALIGNED(Db_y, SIMD_BYTES);
    double* RESTRICT Db_z{ Db_z_ptr() }; ASSUME_ALIGNED(Db_z, SIMD_BYTES);

    #pragma omp simd
    for ( std::size_t i = 0; i < total; ++i ) {
        // Losses in all components:
        double const loss_x{ sig_x[i] * dt_local / ( 2.0 * eps_x[i] ) };
        double const loss_y{ sig_y[i] * dt_local / ( 2.0 * eps_y[i] ) };
        double const loss_z{ sig_z[i] * dt_local / ( 2.0 * eps_z[i] ) };

        // C_a constant:
        Ca_x[i] = ( 1.0 - loss_x ) / ( 1.0 + loss_x );
        Ca_y[i] = ( 1.0 - loss_y ) / ( 1.0 + loss_y );
        Ca_z[i] = ( 1.0 - loss_z ) / ( 1.0 + loss_z );

        double const denom_x{ eps_x[i] / dt_local + sig_x[i] / 2.0 };
        double const denom_y{ eps_y[i] / dt_local + sig_y[i] / 2.0 };
        double const denom_z{ eps_z[i] / dt_local + sig_z[i] / 2.0 };
        
        // C_b constant:
        Cb_x[i] = 1.0 / denom_x;
        Cb_y[i] = 1.0 / denom_y;
        Cb_z[i] = 1.0 / denom_z;

        // D_b constant:
        Db_x[i] = dt_local / mu_x[i];
        Db_y[i] = dt_local / mu_y[i];
        Db_z[i] = dt_local / mu_z[i];
    }
}

void Grid::add_source( std::unique_ptr<Source> source ) {
    sources().emplace_back( std::move( source ) );
}

void Grid::apply_sources( std::size_t const time_step ) {
    for ( auto const &source : sources() ) {
        source->apply( *this, time_step );
    }
}

void Grid::update_H() {
    double* RESTRICT Hx{ Hx_ptr() }; ASSUME_ALIGNED(Hx, SIMD_BYTES);
    double* RESTRICT Hy{ Hy_ptr() }; ASSUME_ALIGNED(Hy, SIMD_BYTES);
    double* RESTRICT Hz{ Hz_ptr() }; ASSUME_ALIGNED(Hz, SIMD_BYTES);

    double* RESTRICT Ex{ Ex_ptr() }; ASSUME_ALIGNED(Ex, SIMD_BYTES);
    double* RESTRICT Ey{ Ey_ptr() }; ASSUME_ALIGNED(Ey, SIMD_BYTES);
    double* RESTRICT Ez{ Ez_ptr() }; ASSUME_ALIGNED(Ez, SIMD_BYTES);

    double* RESTRICT Db_x{ Db_x_ptr() }; ASSUME_ALIGNED(Db_x, SIMD_BYTES);
    double* RESTRICT Db_y{ Db_y_ptr() }; ASSUME_ALIGNED(Db_y, SIMD_BYTES);
    double* RESTRICT Db_z{ Db_z_ptr() }; ASSUME_ALIGNED(Db_z, SIMD_BYTES);

    double const inv_dx{ 1.0 / dx() };
    double const inv_dy{ 1.0 / dy() };
    double const inv_dz{ 1.0 / dz() };

    std::size_t const Nx_local{ Nx() };
    std::size_t const Ny_local{ Ny() };
    std::size_t const Nz_local{ Nz() };

    std::size_t const Sy{ Nx_padded() };
    std::size_t const Sz{ Nx_padded() * Ny_padded() };

    #pragma omp parallel for collapse( 2 ) schedule( static )
    for ( std::size_t z = 0; z < Nz_local - 1; ++z ) {
        for ( std::size_t y = 0; y < Ny_local - 1; ++y ) {
            std::size_t const base{ y * Sy + z * Sz };

            #pragma omp simd
            for ( std::size_t x = 0; x < Nx_local - 1; ++x ) {
                std::size_t const i{ base + x };
                std::size_t const i_x{ i + 1 };
                std::size_t const i_y{ i + Sy };
                std::size_t const i_z{ i + Sz };

                double const curl_x{
                    ( Ez[i_y] - Ez[i] ) * inv_dy
                  - ( Ey[i_z] - Ey[i] ) * inv_dz
                };
                double const curl_y{
                    ( Ex[i_z] - Ex[i] ) * inv_dz
                  - ( Ez[i_x] - Ez[i] ) * inv_dx
                };
                double const curl_z{
                    ( Ey[i_x] - Ey[i] ) * inv_dx
                  - ( Ex[i_y] - Ex[i] ) * inv_dy
                };

                Hx[i] -= Db_x[i] * curl_x;
                Hy[i] -= Db_y[i] * curl_y;
                Hz[i] -= Db_z[i] * curl_z;
            }
        }
    }

    pml().update_H_psi(
        Ex_ptr(), Ey_ptr(), Ez_ptr(),
        Hx_ptr(), Hy_ptr(), Hz_ptr(),
        Db_x_ptr(), Db_y_ptr(), Db_z_ptr(),
        dx(), dy(), dz() 
    );    
}

void Grid::update_E() {
    double* RESTRICT Ex{ Ex_ptr() }; ASSUME_ALIGNED(Ex, SIMD_BYTES);
    double* RESTRICT Ey{ Ey_ptr() }; ASSUME_ALIGNED(Ey, SIMD_BYTES);
    double* RESTRICT Ez{ Ez_ptr() }; ASSUME_ALIGNED(Ez, SIMD_BYTES);

    double* RESTRICT Hx{ Hx_ptr() }; ASSUME_ALIGNED(Hx, SIMD_BYTES);
    double* RESTRICT Hy{ Hy_ptr() }; ASSUME_ALIGNED(Hy, SIMD_BYTES);
    double* RESTRICT Hz{ Hz_ptr() }; ASSUME_ALIGNED(Hz, SIMD_BYTES);

    double* RESTRICT Jx{ Jx_ptr() }; ASSUME_ALIGNED(Jx, SIMD_BYTES);
    double* RESTRICT Jy{ Jy_ptr() }; ASSUME_ALIGNED(Jy, SIMD_BYTES);
    double* RESTRICT Jz{ Jz_ptr() }; ASSUME_ALIGNED(Jz, SIMD_BYTES);

    double* RESTRICT Ca_x{ Ca_x_ptr() }; ASSUME_ALIGNED(Ca_x, SIMD_BYTES);
    double* RESTRICT Ca_y{ Ca_y_ptr() }; ASSUME_ALIGNED(Ca_y, SIMD_BYTES);
    double* RESTRICT Ca_z{ Ca_z_ptr() }; ASSUME_ALIGNED(Ca_z, SIMD_BYTES);

    double* RESTRICT Cb_x{ Cb_x_ptr() }; ASSUME_ALIGNED(Cb_x, SIMD_BYTES);
    double* RESTRICT Cb_y{ Cb_y_ptr() }; ASSUME_ALIGNED(Cb_y, SIMD_BYTES);
    double* RESTRICT Cb_z{ Cb_z_ptr() }; ASSUME_ALIGNED(Cb_z, SIMD_BYTES);

    double const inv_dx{ 1.0 / dx() };
    double const inv_dy{ 1.0 / dy() };
    double const inv_dz{ 1.0 / dz() };

    std::size_t const Nx_local{ Nx() };
    std::size_t const Ny_local{ Ny() };
    std::size_t const Nz_local{ Nz() };

    std::size_t const Sy{ Nx_padded() };
    std::size_t const Sz{ Nx_padded() * Ny_padded() };

    #pragma omp parallel for collapse( 2 ) schedule( static )
    for ( std::size_t z = 1; z < Nz_local - 1; ++z ) {
        for ( std::size_t y = 1; y < Ny_local - 1; ++y ) {
            std::size_t const base{ y * Sy + z * Sz };

            #pragma omp simd
            for ( std::size_t x = 1; x < Nx_local - 1; ++x ) {
                std::size_t const i{ base + x };
                std::size_t const i_x{ i - 1 };
                std::size_t const i_y{ i - Sy };
                std::size_t const i_z{ i - Sz };

                double const curl_x{
                    ( Hz[i] - Hz[i_y] ) * inv_dy
                  - ( Hy[i] - Hy[i_z] ) * inv_dz
                };
                double const curl_y{
                    ( Hx[i] - Hx[i_z] ) * inv_dz
                  - ( Hz[i] - Hz[i_x] ) * inv_dx
                };
                double const curl_z{
                    ( Hy[i] - Hy[i_x] ) * inv_dx
                  - ( Hx[i] - Hx[i_y] ) * inv_dy
                };

                Ex[i] = Ca_x[i] * Ex[i] + Cb_x[i] * ( curl_x - Jx[i] );
                Ey[i] = Ca_y[i] * Ey[i] + Cb_y[i] * ( curl_y - Jy[i] );
                Ez[i] = Ca_z[i] * Ez[i] + Cb_z[i] * ( curl_z - Jz[i] );
            }
        }
    }

    pml().update_E_psi(
        Ex_ptr(), Ey_ptr(), Ez_ptr(),
        Hx_ptr(), Hy_ptr(), Hz_ptr(),
        Cb_x_ptr(), Cb_y_ptr(), Cb_z_ptr(),
        dx(), dy(), dz()
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
    double const dV{ dx() * dy() * dz() };

    double const* RESTRICT Ex{ Ex_ptr() }; ASSUME_ALIGNED(Ex, SIMD_BYTES);
    double const* RESTRICT Ey{ Ey_ptr() }; ASSUME_ALIGNED(Ey, SIMD_BYTES);
    double const* RESTRICT Ez{ Ez_ptr() }; ASSUME_ALIGNED(Ez, SIMD_BYTES);

    double const* RESTRICT eps_x{ eps_x_ptr() }; ASSUME_ALIGNED(eps_x, SIMD_BYTES);
    double const* RESTRICT eps_y{ eps_y_ptr() }; ASSUME_ALIGNED(eps_y, SIMD_BYTES);
    double const* RESTRICT eps_z{ eps_z_ptr() }; ASSUME_ALIGNED(eps_z, SIMD_BYTES);

    std::size_t const Nx_local{ Nx() };
    std::size_t const Ny_local{ Ny() };
    std::size_t const Nz_local{ Nz() };

    std::size_t const Sy{ Nx_padded() };
    std::size_t const Sz{ Nx_padded() * Ny_padded() };

    #pragma omp parallel for collapse( 2 ) reduction( +:energy )
    for ( std::size_t z = 1; z < Nz_local - 1; ++z ) {
        for ( std::size_t y = 1; y < Ny_local - 1; ++y ) {
            std::size_t const base{ y * Sy + z * Sz };

            for ( std::size_t x = 1; x < Nx_local - 1; ++x ) {
                std::size_t const i{ base + x };
                std::size_t const i_x{ i + 1 };
                std::size_t const i_y{ i + Sy };
                std::size_t const i_z{ i + Sz };

                double const Ex_avg{ ( Ex[i] + Ex[i_x] ) * 0.5 };
                double const Ey_avg{ ( Ey[i] + Ey[i_y] ) * 0.5 };
                double const Ez_avg{ ( Ez[i] + Ez[i_z] ) * 0.5 };

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
    double const dV{ dx() * dy() * dz() };

    double const* RESTRICT Hx{ Hx_ptr() }; ASSUME_ALIGNED(Hx, SIMD_BYTES);
    double const* RESTRICT Hy{ Hy_ptr() }; ASSUME_ALIGNED(Hy, SIMD_BYTES);
    double const* RESTRICT Hz{ Hz_ptr() }; ASSUME_ALIGNED(Hz, SIMD_BYTES);

    double const* RESTRICT mu_x{ mu_x_ptr() }; ASSUME_ALIGNED(mu_x, SIMD_BYTES);
    double const* RESTRICT mu_y{ mu_y_ptr() }; ASSUME_ALIGNED(mu_y, SIMD_BYTES);
    double const* RESTRICT mu_z{ mu_z_ptr() }; ASSUME_ALIGNED(mu_z, SIMD_BYTES);

    std::size_t const Nx_local{ Nx() };
    std::size_t const Ny_local{ Ny() };
    std::size_t const Nz_local{ Nz() };

    std::size_t const Sy{ Nx_padded() };
    std::size_t const Sz{ Nx_padded() * Ny_padded() };

    #pragma omp parallel for collapse( 2 ) reduction( +:energy )
    for ( std::size_t z = 1; z < Nz_local - 1; ++z ) {
        for ( std::size_t y = 1; y < Ny_local - 1; ++y ) {
            std::size_t const base{ y * Sy + z * Sz };

            for ( std::size_t x = 1; x < Nx_local - 1; ++x ) {
                std::size_t const i{ base + x };
                std::size_t const i_x{ i + 1 };
                std::size_t const i_y{ i + Sy };
                std::size_t const i_z{ i + Sz };

                std::size_t const i_x_y{ i_x + Sy };
                std::size_t const i_x_z{ i_x + Sz };
                std::size_t const i_y_z{ i_y + Sz };

                double const Hx_avg{ ( Hx[i] + Hx[i_y] + Hx[i_z] + Hx[i_y_z] ) * 0.25 };
                double const Hy_avg{ ( Hy[i] + Hy[i_x] + Hy[i_z] + Hy[i_x_z] ) * 0.25 };
                double const Hz_avg{ ( Hz[i] + Hz[i_x] + Hz[i_y] + Hz[i_x_y] ) * 0.25 };

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