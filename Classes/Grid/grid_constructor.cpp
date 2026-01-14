#include "grid.hpp"

// Constructor:
Grid::Grid( std::size_t new_Nx, std::size_t new_Ny, std::size_t new_Nz,
            double new_dx, double new_dy, double new_dz,
            double new_eps, double new_mu ):
Nx_{ new_Nx }, Ny_{ new_Ny }, Nz_{ new_Nz },
dx_{ new_dx }, dy_{ new_dy }, dz_{ new_dz },
eps_{ new_eps }, mu_{ new_mu },
c_{ 1.0 / std::sqrt( mu() * eps() ) },
dt_{ config::cfl_factor / ( c() * std::sqrt( 1.0/(dx()*dx()) + 1.0/(dy()*dy()) + 1.0/(dz()*dz()) ) ) } {
    
    std::size_t const grid_size{ Nx() * Ny() * Nz() };
    Ex_ = std::make_unique<double[]>( grid_size );
    Ey_ = std::make_unique<double[]>( grid_size );
    Ez_ = std::make_unique<double[]>( grid_size );

    Bx_ = std::make_unique<double[]>( grid_size );
    By_ = std::make_unique<double[]>( grid_size );
    Bz_ = std::make_unique<double[]>( grid_size );

    Jx_ = std::make_unique<double[]>( grid_size );
    Jy_ = std::make_unique<double[]>( grid_size );
    Jz_ = std::make_unique<double[]>( grid_size );
}