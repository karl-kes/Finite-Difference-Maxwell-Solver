#include "grid.hpp"

// Constructor:
Grid::Grid( std::size_t new_Nx, std::size_t new_Ny, std::size_t new_Nz,
            double new_dx, double new_dy, double new_dz,
            double new_eps, double new_mu ):
Nx_{ new_Nx }, Ny_{ new_Ny }, Nz_{ new_Nz },
dx_{ new_dx }, dy_{ new_dy }, dz_{ new_dz },
eps_{ new_eps }, mu_{ new_mu } {

    c_ = 1.0 / std::sqrt( mu() * eps() );
    dt_ = 0.1 / ( c() * std::sqrt( 1.0/(dx()*dx()) + 1.0/(dy()*dy()) + 1.0/(dz()*dz()) ) );
    
    std::size_t const grid_size{ Nx_ * Ny_ * Nz_ };
    Ex_ = std::make_unique<double[]>( grid_size );
    Ey_ = std::make_unique<double[]>( grid_size );
    Ez_ = std::make_unique<double[]>( grid_size );
    Bx_ = std::make_unique<double[]>( grid_size );
    By_ = std::make_unique<double[]>( grid_size );
    Bz_ = std::make_unique<double[]>( grid_size );
}

// System Simulation:
void Grid::update_B() {
    // dB/dt = -curl( E )
    #pragma omp parallel for collapse( 2 )
    for ( std::size_t z = 0; z < Nz() - 1; ++z ) {
        for ( std::size_t y = 0; y < Ny() - 1; ++y ) {
            #pragma omp simd
            for ( std::size_t x = 0; x < Nx() - 1; ++x ) {
                // ∂Bx/∂t = -(∂Ez/∂y - ∂Ey/∂z)
                Bx_[idx(x,y,z)] -= dt() * curl_X( Ey_[idx(x,y,z)], Ey_[idx(x,y,z+1)],
                                                  Ez_[idx(x,y,z)], Ez_[idx(x,y+1,z)] );

                // ∂By/∂t = -(∂Ex/∂z - ∂Ez/∂x)
                By_[idx(x,y,z)] -= dt() * curl_Y( Ex_[idx(x,y,z)], Ex_[idx(x,y,z+1)],
                                                  Ez_[idx(x,y,z)], Ez_[idx(x+1,y,z)] );

                // ∂Bz/∂t = -(∂Ey/∂x - ∂Ex/∂y)
                Bz_[idx(x,y,z)] -= dt() * curl_Z( Ey_[idx(x,y,z)], Ey_[idx(x+1,y,z)],
                                                  Ex_[idx(x,y,z)], Ex_[idx(x,y+1,z)] );
            }
        }
    }
}
void Grid::update_E() {
    // dE/dt = c*c * curl(B)
    #pragma omp parallel for collapse( 2 )
    for ( std::size_t z = 1; z < Nz(); ++z ) {
        for ( std::size_t y = 1; y < Ny(); ++y ) {
            #pragma omp simd
            for ( std::size_t x = 1; x < Nx(); ++x ) {
                // ∂Ex/∂t = c*c * (∂Ez/∂y - ∂Ey/∂z)
                Ex_[idx(x,y,z)] += dt() * c_sq() * curl_X( By_[idx(x,y,z-1)], By_[idx(x,y,z)],
                                                           Bz_[idx(x,y-1,z)], Bz_[idx(x,y,z)] );

                // ∂Ey/∂t = c*c * (∂Ex/∂z - ∂Ez/∂x)
                Ey_[idx(x,y,z)] += dt() * c_sq() * curl_Y( Bx_[idx(x,y,z-1)], Bx_[idx(x,y,z)],
                                                           Bz_[idx(x-1,y,z)], Bz_[idx(x,y,z)] );

                // ∂Ez/∂t = c*c * (∂Ex/∂y - ∂Ey/∂x)
                Ez_[idx(x,y,z)] += dt() * c_sq() * curl_Z( By_[idx(x-1,y,z)], By_[idx(x,y,z)],
                                                           Bx_[idx(x,y-1,z)], Bx_[idx(x,y,z)] );
            }
        }
    }
}
void Grid::step() {
    update_B();
    update_E();
}
void Grid::inject_source( std::size_t const x,
                          std::size_t const y,
                          std::size_t const z, 
                          double const value ) {
    Ez_[idx(x,y,z)] += value;
}
void Grid::component_slice( std::size_t const z,
                            std::string const &file_name,
                            char const field,
                            char const component ) {
    std::ofstream file( file_name );
    for ( std::size_t y{}; y < Ny(); ++ y ) {
        for ( std::size_t x{}; x< Nx(); ++x ) {
            file << get_field( field, component, x, y, z );
            
            if ( x < Nx() - 1 ) {
                file << ",";
            }
        }
        file << "\n";
    }
    file.close();
}
void Grid::magnitude_slice( std::size_t const z,
                            std::string const &file_name,
                            char const field ) {
    std::ofstream file ( file_name );
    for ( std::size_t y{}; y < Ny(); ++y ) {
        for ( std::size_t x{}; x < Nx(); ++ x ) {
            double mag{ field_mag( field, x, y, z )};
            file << mag;

            if ( x < Nx() - 1 ) {
                file << ",";
            }
        }
        file << "\n";
    }
    file.close();
}
void Grid::magnitude_volume( std::string const &file_name, char const field ) {
    std::ofstream file( file_name );

    for ( std::size_t z = 0; z < Nz(); ++z ) {
        for ( std::size_t y = 0; y < Ny(); ++y ) {
            for ( std::size_t x = 0; x < Nx(); ++x ) {
                double mag{ field_mag( field, x, y, z ) };

                file << x << "," << y << "," << z << "," << mag << "\n";
            }
        }
    }
    file.close();
}
void Grid::vector_volume( std::string const &file_name, char const field ) {
    std::ofstream file( file_name, std::ios::binary | std::ios::out );

    uint64_t dimension[3]{ (uint64_t)(Nx()-1), (uint64_t)(Ny()-1), (uint64_t)(Nz()-1) };
    file.write( reinterpret_cast<const char*>( dimension ), sizeof( dimension ) );

    std::vector<double> buffer{};
    buffer.reserve( (Nx()-1) * (Ny()-1) * 4 );

    for ( std::size_t z = 0; z < Nz() - 1; ++z ) {
        buffer.clear();
        for ( std::size_t y = 0; y < Ny() - 1; ++y ) {
            for ( std::size_t x = 0; x < Nx() - 1; ++x ) {
                double Fx_avg{ 0.5 * ( get_field( field, 'x', x, y, z ) + get_field( field, 'x', x+1, y, z ) ) };
                double Fy_avg{ 0.5 * ( get_field( field, 'y', x, y, z ) + get_field( field, 'y', x, y+1, z ) ) };
                double Fz_avg{ 0.5 * ( get_field( field, 'z', x, y, z ) + get_field( field, 'z', x, y, z+1 ) ) };

                buffer.push_back( Fx_avg );
                buffer.push_back( Fy_avg );
                buffer.push_back( Fz_avg );
                buffer.push_back( std::sqrt( Fx_avg*Fx_avg + Fy_avg*Fy_avg + Fz_avg*Fz_avg ) );
            }
        }
        file.write( reinterpret_cast<const char*>( buffer.data() ), buffer.size() * sizeof( double ) );
    }
    file.close();
}

// Getters:
// Differentials
double Grid::dx() const {
    return dx_;
}
double Grid::dy() const {
    return dy_;
}
double Grid::dz() const {
    return dz_;
}
double Grid::dt() const {
    return dt_;
}
// Speed, Mu, Epsilon
double Grid::c() const {
    return c_;
}
double Grid::c_sq() const {
    return c_*c_;
}
double Grid::eps() const {
    return eps_;
}
double Grid::mu() const {
    return mu_;
}
// Dimensions
std::size_t Grid::Nx() const {
    return Nx_;
}
std::size_t Grid::Ny() const {
    return Ny_;
}
std::size_t Grid::Nz() const {
    return Nz_;
}
// Fields
double Grid::get_field( char const field,
                        char const component,
                        std::size_t const x,
                        std::size_t const y,
                        std::size_t const z ) const {
    std::size_t index{ idx(x,y,z) };

    if ( std::tolower( field ) == 'e' ) {
        switch ( std::tolower( component ) ) {
            case 'x': return Ex_[index];
            case 'y': return Ey_[index];
            case 'z': return Ez_[index];
        }
    } else if ( std::tolower( field ) == 'b' ) {
        switch ( std::tolower( component ) ) {
            case 'x': return Bx_[index];
            case 'y': return By_[index];
            case 'z': return Bz_[index];
        }
    }
    std::cout << "ERROR! CHECK PARAMETERS!" << std::endl;
    return 0.0;
}
double Grid::field_mag( char const field,
                        std::size_t const x,
                        std::size_t const y,
                        std::size_t const z ) const {
    double Fx{ get_field( field, 'x', x, y, z ) };
    double Fy{ get_field( field, 'y', x, y, z ) };
    double Fz{ get_field( field, 'z', x, y, z ) };

    return std::sqrt( Fx*Fx + Fy*Fy + Fz*Fz );
}

// Helpers:
// Finds 3D index
std::size_t Grid::idx( std::size_t const x,
                       std::size_t const y,
                       std::size_t const z ) const {
    return x + Nx() * ( y + Ny() * z );
}
// Curls in X, Y, Z
double Grid::curl_X( double const Y_0, double const Y_1,
                     double const Z_0, double const Z_1 ) const {
    double dZdY{ ( Z_1 - Z_0 ) / dy() };
    double dYdZ{ ( Y_1 - Y_0 ) / dz() };

    return ( dZdY - dYdZ );
}
double Grid::curl_Y( double const X_0, double const X_1,
                     double const Z_0, double const Z_1 ) const {
    double dXdZ{ ( X_1 - X_0 ) / dz() };
    double dZdX{ ( Z_1 - Z_0 ) / dx() };

    return ( dXdZ - dZdX );
}
double Grid::curl_Z( double const Y_0, double const Y_1,
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