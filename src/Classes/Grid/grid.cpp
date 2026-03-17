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
, N_{ Nx_*Ny_*Nz_ }
, dx_{ config.dx }
, dy_{ config.dy }
, dz_{ config.dz }
, eps_{ config.eps }
, mu_{ config.mu }
, c_{ config.c }
, c_sq_{ config.c * config.c }
, dt_{ config.dt }
, fields_{ Nx_padded_ * Ny_padded_ * Nz_padded_, NUM_GRID_ARRAYS_ }
, pml_{ config }
{ }

Grid::~Grid() = default;

void Grid::add_source( std::unique_ptr<Source> source ) {
    sources_.push_back( std::move( source ) );
}

void Grid::apply_sources( std::size_t const time_step ) {
    for ( auto const &source : sources_ ) {
        source->apply( *this, time_step );
    }
}

void Grid::update_B() {
    double* RESTRICT Bx{ Bx_ptr() };
    double* RESTRICT By{ By_ptr() };
    double* RESTRICT Bz{ Bz_ptr() };
    double* RESTRICT Ex{ Ex_ptr() };
    double* RESTRICT Ey{ Ey_ptr() };
    double* RESTRICT Ez{ Ez_ptr() };

    double const dt_local{ dt_ };
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

                Bx[i] -= dt_local * curl_x;
                By[i] -= dt_local * curl_y;
                Bz[i] -= dt_local * curl_z;
            }
        }
    }
    pml_.update_B_psi( Ex_ptr(), Ey_ptr(), Ez_ptr(),
                       Bx_ptr(), By_ptr(), Bz_ptr(),
                       dt_, dx_, dy_, dz_ );
}

void Grid::update_E() {
    double* RESTRICT Bx{ Bx_ptr() };
    double* RESTRICT By{ By_ptr() };
    double* RESTRICT Bz{ Bz_ptr() };
    double* RESTRICT Ex{ Ex_ptr() };
    double* RESTRICT Ey{ Ey_ptr() };
    double* RESTRICT Ez{ Ez_ptr() };
    double* RESTRICT Jx{ Jx_ptr() };
    double* RESTRICT Jy{ Jy_ptr() };
    double* RESTRICT Jz{ Jz_ptr() };

    double const dt_local{ dt_ };
    double const csq{ c_sq_ };
    double const inv_eps{ 1.0 / eps_ };
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
                    ( Bz[i] - Bz[i - Sy] ) * inv_dy
                  - ( By[i] - By[i - Sz] ) * inv_dz
                };
                double const curl_y{
                    ( Bx[i] - Bx[i - Sz] ) * inv_dz
                  - ( Bz[i] - Bz[i - 1]  ) * inv_dx
                };
                double const curl_z{
                    ( By[i] - By[i - 1]  ) * inv_dx
                  - ( Bx[i] - Bx[i - Sy] ) * inv_dy
                };

                Ex[i] += dt_local * ( csq * curl_x - Jx[i] * inv_eps );
                Ey[i] += dt_local * ( csq * curl_y - Jy[i] * inv_eps );
                Ez[i] += dt_local * ( csq * curl_z - Jz[i] * inv_eps );
            }
        }
    }
    pml_.update_E_psi( Ex_ptr(), Ey_ptr(), Ez_ptr(),
                       Bx_ptr(), By_ptr(), Bz_ptr(),
                       dt_, dx_, dy_, dz_, c_sq_ );
}

void Grid::step() {
    update_B();
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
            case Component::X: return Bx_ptr()[i];
            case Component::Y: return By_ptr()[i];
            case Component::Z: return Bz_ptr()[i];
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
            case Component::X: return Bx_ptr()[i];
            case Component::Y: return By_ptr()[i];
            case Component::Z: return Bz_ptr()[i];
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

double Grid::total_energy() const {
    double energy{};
    double const dV{ dx_ * dy_ * dz_ };

    double const* RESTRICT Ex{ Ex_ptr() };
    double const* RESTRICT Ey{ Ey_ptr() };
    double const* RESTRICT Ez{ Ez_ptr() };
    double const* RESTRICT Bx{ Bx_ptr() };
    double const* RESTRICT By{ By_ptr() };
    double const* RESTRICT Bz{ Bz_ptr() };

    double const inv_mu{ 1.0 / mu_ };
    double const eps_local{ eps_ };

    std::size_t const Sy{ Nx_padded_ };
    std::size_t const Sz{ Nx_padded_ * Ny_padded_ };

    #pragma omp parallel for collapse( 2 ) reduction( +:energy )
    for ( std::size_t z = 0; z < Nz_; ++z ) {
        for ( std::size_t y = 0; y < Ny_; ++y ) {
            std::size_t const base{ y * Sy + z * Sz };

            for ( std::size_t x = 0; x < Nx_; ++x ) {
                std::size_t const i{ base + x };

                double const E_sq{ Ex[i]*Ex[i] + Ey[i]*Ey[i] + Ez[i]*Ez[i] };
                double const B_sq{ Bx[i]*Bx[i] + By[i]*By[i] + Bz[i]*Bz[i] };

                energy += 0.5 * ( eps_local * E_sq + B_sq * inv_mu );
            }
        }
    }
    return energy * dV;
}

double Grid::source_power() const {
    double power{};
    double const dV{ dx_ * dy_ * dz_ };

    double const* RESTRICT Ex{ Ex_ptr() };
    double const* RESTRICT Ey{ Ey_ptr() };
    double const* RESTRICT Ez{ Ez_ptr() };
    double const* RESTRICT Jx{ Jx_ptr() };
    double const* RESTRICT Jy{ Jy_ptr() };
    double const* RESTRICT Jz{ Jz_ptr() };

    std::size_t const Sy{ Nx_padded_ };
    std::size_t const Sz{ Nx_padded_ * Ny_padded_ };

    #pragma omp parallel for collapse( 2 ) reduction( +:power )
    for ( std::size_t z = 0; z < Nz_; ++z ) {
        for ( std::size_t y = 0; y < Ny_; ++y ) {
            std::size_t const base{ y * Sy + z * Sz };

            for ( std::size_t x = 0; x < Nx_; ++x ) {
                std::size_t const i{ base + x };

                power -= Jx[i] * Ex[i] + Jy[i] * Ey[i] + Jz[i] * Ez[i];
            }
        }
    }
    return power * dV;
}