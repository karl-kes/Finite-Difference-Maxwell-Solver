#include "grid.hpp"

// Helpers:
void Grid::print_progress( int curr_time, int total_time ) const {
    double percent{ 100.0 * curr_time / total_time };
    std::cout << "\rProgress: " << percent << "%" << std::flush;
    
    if ( curr_time == total_time ) {
        std::cout << std::endl;
    }
}

// Finds 3D index
std::size_t Grid::idx( std::size_t const x, std::size_t const y, std::size_t const z ) const {
    return x + Nx() * ( y + Ny() * z );
}

// Curls in X, Y, Z
double Grid::curl_x( double const Y_0, double const Y_1,
                     double const Z_0, double const Z_1 ) const {
    double dZdY{ ( Z_1 - Z_0 ) / dy() };
    double dYdZ{ ( Y_1 - Y_0 ) / dz() };

    return ( dZdY - dYdZ );
}

double Grid::curl_y( double const X_0, double const X_1,
                     double const Z_0, double const Z_1 ) const {
    double dXdZ{ ( X_1 - X_0 ) / dz() };
    double dZdX{ ( Z_1 - Z_0 ) / dx() };

    return ( dXdZ - dZdX );
}

double Grid::curl_z( double const Y_0, double const Y_1,
                     double const X_0, double const X_1 ) const {
    double dYdX{ ( Y_1 - Y_0 ) / dx() };
    double dXdY{ ( X_1 - X_0 ) / dy() };

    return ( dYdX - dXdY );
}

double Grid::total_energy() const {
    double energy{};
    double dV{ dx() * dy() * dz() };

    #pragma omp parallel for collapse( 2 ) reduction( +:energy )
    for ( std::size_t z = 0; z < Nz(); ++z ) {
        for ( std::size_t y = 0; y < Ny(); ++y ) {
            #pragma omp simd
            for ( std::size_t x = 0; x < Nx(); ++x ) {
                std::size_t i = idx(x,y,z);

                double E_sq{ Ex_[i]*Ex_[i] + Ey_[i]*Ey_[i] + Ez_[i]*Ez_[i] };
                double B_sq{ Bx_[i]*Bx_[i] + By_[i]*By_[i] + Bz_[i]*Bz_[i] };

                // Electromagnetic Energy Density:
                // 1/2 * e_0 * E^2 + 1/(2mu_0) * B^2
                energy += 0.5 * ( eps() * E_sq + B_sq / mu() );
            }
        }
    }
    return energy * dV;
}

void Grid::create_directories() const {
    // Clear previous and create new output folder.
    std::filesystem::remove_all("output");
    std::filesystem::create_directories("output/E");
    std::filesystem::create_directories("output/B");
}