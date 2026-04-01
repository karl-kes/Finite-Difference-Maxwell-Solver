#include "source.hpp"

#include <numbers>
#include <cmath>

void Straight_Wire_X::apply( Grid &grid, std::size_t const time_step ) {
    double const omega{ 2.0 * std::numbers::pi * frequency_ };
    double const current{ amplitude_ * std::sin( omega * static_cast<double>( time_step ) * grid.dt() ) };

    double* RESTRICT Jx{ grid.Jx_ptr() };

    for ( std::size_t x{ x_start_ }; x < x_end_; ++x ) {
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

    // Bottom edge and top edge:
    for ( std::size_t x{ x_i }; x < x_f; ++x ) {
        Jx[grid.idx( x, y_i, z_ )] = current;
        Jx[grid.idx( x, y_f, z_ )] = -current;
    }

    // Right edge and left edge:
    for ( std::size_t y{ y_i }; y < y_f; ++y ) {
        Jy[grid.idx( x_f, y, z_ )] = current;
        Jy[grid.idx( x_i, y, z_ )] = -current;
    }
}

void AC_Concentric_Rings::apply( Grid &grid, std::size_t const time_step ) {
    double const omega{ 2.0 * std::numbers::pi * frequency_ };
    double const current{ amplitude_ * std::sin( omega * static_cast<double>( time_step ) * grid.dt() ) };

    double* RESTRICT Jx{ grid.Jx_ptr() };
    double* RESTRICT Jy{ grid.Jy_ptr() };

    // Define the outermost boundary of the largest rings
    std::size_t const x_outer_min{ grid.Nx() / 4 };
    std::size_t const x_outer_max{ 3 * grid.Nx() / 4 };

    std::size_t const y_outer_min{ grid.Ny() / 4 };
    std::size_t const y_outer_max{ 3 * grid.Ny() / 4 };

    // Total number of rings
    std::size_t const num_rings{ 3 };

    // Calculate the distance to shrink inward for each successive ring.
    // Dividing by (2 * numrings) ensures the innermost ring doesn't collapse into a single point:
    std::size_t const dx{ (x_outer_max - x_outer_min) / (2 * num_rings) };
    std::size_t const dy{ (y_outer_max - y_outer_min) / (2 * num_rings) };

    // Draw each ring, starting from the outside ( i = 0 ) and moving inwards
    for ( std::size_t i{}; i < num_rings; ++i ) {
        
        // Calculate the specific bounds for the current ring
        std::size_t const cur_x_min{ x_outer_min + i * dx };
        std::size_t const cur_x_max{ x_outer_max - i * dx };

        std::size_t const cur_y_min{ y_outer_min + i * dy };
        std::size_t const cur_y_max{ y_outer_max - i * dy };

        // Bottom edge (+x) and Top edge (-x)
        for ( std::size_t x{ cur_x_min }; x < cur_x_max; ++x ) {
            Jx[grid.idx( x, cur_y_min, z_ )] = current;
            Jx[grid.idx( x, cur_y_max, z_ )] = -current;
        }

        // Left edge (-y) and Right edge (+y)
        for ( std::size_t y{ cur_y_min }; y < cur_y_max; ++y ) {
            Jy[grid.idx( cur_x_min, y, z_ )] = -current;
            Jy[grid.idx( cur_x_max, y, z_ )] = current;
        }
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
    std::size_t const i{ grid.idx(x_,y_,z_) };
    double* RESTRICT Jz{ grid.Jz_ptr() };

    double const t{ static_cast<double>( time_step ) * grid.dt() };
    double const factor{ ( t - t_0_ ) / width_ };
    double const exponent{ -0.5 * factor * factor };
    
    Jz[i] = amplitude_ * std::exp( exponent );
}