#include "source.hpp"

#include <numbers>
#include <cmath>

void Straight_Wire_X::apply( Grid &grid, std::size_t const time_step ) {
    double const omega{ 2.0 * std::numbers::pi * frequency_ };
    double const current{ amplitude_ * std::sin( omega * static_cast<double>( time_step ) * grid.dt() ) };

    double* RESTRICT Jx{ grid.Jx_ptr() };

    for ( std::size_t x{ x_start_ }; x <= x_end_; ++x ) {
        std::size_t const i{ grid.idx(x,y_,z_) };
        Jx[i] = current;
    }
}

void AC_Current_Loop::apply( Grid &grid, std::size_t const time_step ) {
    double const omega{ 2.0 * std::numbers::pi * frequency_ };
    double const current{ amplitude_ * std::sin( omega * static_cast<double>( time_step ) * grid.dt() ) };

    double* RESTRICT Jx{ grid.Jx_ptr() };
    double* RESTRICT Jy{ grid.Jy_ptr() };

    std::size_t const x_f{ 3 * grid.Nx() / 4 };
    std::size_t const x_i{ grid.Nx() / 4 };

    std::size_t const y_f{ 3 * grid.Ny() / 4 };
    std::size_t const y_i{ grid.Ny() / 4 };

    // Bottom edge: +x at y = y_i
    for ( std::size_t x{ x_i }; x < x_f; ++x ) {
        Jx[grid.idx( x, y_i, z_ )] = current;
    }

    // Top edge: -x at y = y_f
    for ( std::size_t x{ x_i }; x < x_f; ++x ) {
        Jx[grid.idx( x, y_f, z_ )] = -current;
    }

    // Right edge: +y at x = x_f
    for ( std::size_t y{ y_i }; y < y_f; ++y ) {
        Jy[grid.idx( x_f, y, z_ )] = current;
    }

    // Left edge: -y at x = x_i
    for ( std::size_t y{ y_i }; y < y_f; ++y ) {
        Jy[grid.idx( x_i, y, z_ )] = -current;
    }
}

void Point_Source::apply( Grid &grid, std::size_t const time_step ) {
    ( void )time_step;
    std::size_t const i{ grid.idx(x_,y_,z_) };

    double* RESTRICT Jx{ grid.Jx_ptr() };
    double* RESTRICT Jy{ grid.Jy_ptr() };
    double* RESTRICT Jz{ grid.Jz_ptr() };

    Jx[i] = value_;
    Jy[i] = value_;
    Jz[i] = value_;
}

void Gaussian_Pulse::apply( Grid& grid, std::size_t const time_step ) {
    double const t{ static_cast<double>( time_step ) * grid.dt() };
    double const exponent{
        -0.5 * ( ( t - t_0_ ) / width_ ) * ( ( t - t_0_ ) / width_ )
    };
    std::size_t const i{ grid.idx(x_,y_,z_) };

    double* RESTRICT Jz{ grid.Jz_ptr() };
    
    Jz[i] = amplitude_ * std::exp( exponent );
}