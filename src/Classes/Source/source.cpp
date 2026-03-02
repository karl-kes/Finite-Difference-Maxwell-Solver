#include "../Grid/grid.hpp"
#include "source.hpp"

void Straight_Wire_X::apply( Grid &grid, double const time_step ) {
    double const omega{ 2.0 * config::PI * frequency_ };
    double const current{ amplitude_ * std::sin( omega * time_step * grid.dt() ) };

    double* RESTRICT Jx{ grid.Jx_() };

    for ( std::size_t x{ x_start_ }; x <= x_end_; ++x ) {
        std::size_t const i{ grid.idx(x,y_,z_) };
        Jx[i] = current;
    }
}

void Point_Source::apply( Grid &grid, double const time_step ) {
    ( void )time_step;
    std::size_t const i{ grid.idx(x_,y_,z_) };

    double* RESTRICT Jx{ grid.Jx_() };
    double* RESTRICT Jy{ grid.Jy_() };
    double* RESTRICT Jz{ grid.Jz_() };

    Jx[i] = value_;
    Jy[i] = value_;
    Jz[i] = value_;
}

void Gaussian_Pulse::apply( Grid& grid, double const time_step ) {
    double const exponent{
        -0.5 * ( ( time_step * grid.dt() - t_0_ ) / width_ ) * ( ( time_step * grid.dt() - t_0_ ) / width_ )
    };
    std::size_t const i{ grid.idx(x_,y_,z_) };

    double* RESTRICT Jz{ grid.Jz_() };
    
    Jz[i] = amplitude_ * std::exp( exponent );
}