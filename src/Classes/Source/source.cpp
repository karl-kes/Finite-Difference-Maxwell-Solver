#include "../Grid/grid.hpp"
#include "source.hpp"

void Straight_Wire_X::apply( Grid &grid, double const time_step ) {
    double const omega{ 2.0 * config::PI * frequency_ };
    double const current{ amplitude_ * std::sin( omega * time_step * grid.dt() ) };

    for ( std::size_t x{ x_start_ }; x <= x_end_; ++x ) {
        grid.Jx(x,y_,z_) = current;
    }
}

void Point_Source::apply( Grid &grid, double const time_step ) {
    ( void )time_step;

    grid.Jx(x_,y_,z_) = value_;
    grid.Jy(x_,y_,z_) = value_;
    grid.Jz(x_,y_,z_) = value_;
}

void Gaussian_Pulse::apply( Grid& grid, double const time_step ) {
    double const exponent{ -0.5 * ( ( time_step * grid.dt() - t_0_ ) / width_ ) * ( ( time_step * grid.dt() - t_0_ ) / width_ ) };
    
    grid.Jz(x_,y_,z_) = amplitude_ * std::exp( exponent );
}