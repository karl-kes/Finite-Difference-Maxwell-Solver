#include "grid.hpp"

// System Simulation:
void Grid::update_B() {
    // ∂B/∂t = -curl( E )
    #pragma omp parallel for collapse( 2 )
    // Start at 0; stagger for Yee-Cell grid.
    for ( std::size_t z = 0; z < Nz() - 1; ++z ) {
        for ( std::size_t y = 0; y < Ny() - 1; ++y ) {
            #pragma omp simd
            for ( std::size_t x = 0; x < Nx() - 1; ++x ) {
                // Take curl of components and apply B -= ∂B:
                // ∂B_x = ∂t * curl_x( E )
                Bx_[idx(x,y,z)] -= dt() * curl_x( Ey_[idx(x,y,z)], Ey_[idx(x,y,z+1)],
                                                  Ez_[idx(x,y,z)], Ez_[idx(x,y+1,z)] );

                // ∂B_y = ∂t * curl_y( E )
                By_[idx(x,y,z)] -= dt() * curl_y( Ex_[idx(x,y,z)], Ex_[idx(x,y,z+1)],
                                                  Ez_[idx(x,y,z)], Ez_[idx(x+1,y,z)] );

                // ∂B_z = ∂t * curl_z( E )
                Bz_[idx(x,y,z)] -= dt() * curl_z( Ey_[idx(x,y,z)], Ey_[idx(x+1,y,z)],
                                                  Ex_[idx(x,y,z)], Ex_[idx(x,y+1,z)] );
            }
        }
    }
}

void Grid::update_E() {
    // ∂E/∂t = c*c * curl(B)
    #pragma omp parallel for collapse( 2 )
    // Start at 1; stagger for Yee-Cell grid.
    for ( std::size_t z = 1; z < Nz(); ++z ) {
        for ( std::size_t y = 1; y < Ny(); ++y ) {
            #pragma omp simd
            for ( std::size_t x = 1; x < Nx(); ++x ) {
                // Curl of components and apply E += ∂E:
                // ∂E_x = ∂t * c*c * (∂E_z/∂y - ∂E_y/∂z)
                Ex_[idx(x,y,z)] += dt() * ( c_sq() * curl_x( By_[idx(x,y,z-1)], By_[idx(x,y,z)],
                                                             Bz_[idx(x,y-1,z)], Bz_[idx(x,y,z)] )
                                                    - Jx_[idx(x,y,z)] / eps() );

                // ∂E_y = ∂t * c*c * (∂Ex/∂z - ∂Ez/∂x)
                Ey_[idx(x,y,z)] += dt() * ( c_sq() * curl_y( Bx_[idx(x,y,z-1)], Bx_[idx(x,y,z)],
                                                             Bz_[idx(x-1,y,z)], Bz_[idx(x,y,z)] )
                                                    - Jy_[idx(x,y,z)] / eps() );

                // ∂E_z = ∂t * c*c * (∂Ex/∂y - ∂Ey/∂x)
                Ez_[idx(x,y,z)] += dt() * ( c_sq() * curl_z( By_[idx(x-1,y,z)], By_[idx(x,y,z)],
                                                             Bx_[idx(x,y-1,z)], Bx_[idx(x,y,z)] )
                                                    - Jz_[idx(x,y,z)] / eps() ) ;
            }
        }
    }
}

void Grid::step() {
    update_B();
    update_E();
}

void Grid::straight_wire_x( double const current,
                            std::size_t const y, std::size_t const z ) {
    for ( std::size_t x{1}; x < Nx()-1; ++x ) {
        Jx_[idx(x,y,z)] = current;
    }
}

void Grid::hard_source_inject( double const value,
                               std::size_t const x, std::size_t const y, std::size_t const z ) {
    Ex_[idx(x,y,z)] += value;
    Ey_[idx(x,y,z)] += value;
    Ez_[idx(x,y,z)] += value;
}

void Grid::soft_source_inject( double const injection,
                               std::size_t const x, std::size_t const y, std::size_t const z ) {
    
}

void Grid::dipole_antenna_inject( double const amp_one, double const amp_two,
                                  double const freq_one, double const freq_two,
                                  double const injection,
                                  std::size_t const x, std::size_t const y, std::size_t const z ) {
    Ex_[idx(x,y,z)] += amp_one * std::sin( freq_one * injection );
    Ey_[idx(x,y,z)] += amp_one * std::sin( freq_one * injection );
    Ez_[idx(x,y,z)] += amp_one * std::sin( freq_one * injection );

    Ex_[idx(Nx()-x,Ny()-y,Nz()-z)] += amp_two * std::sin( freq_two * injection );
    Ey_[idx(Nx()-x,Ny()-y,Nz()-z)] += amp_two * std::sin( freq_two * injection );
    Ez_[idx(Nx()-x,Ny()-y,Nz()-z)] += amp_two * std::sin( freq_two * injection );
}

void Grid::gaussian_pulse_inject( double const injection,
                                  std::size_t const x, std::size_t const y, std::size_t const z ) {

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
                // Average to fix the staggering for visual output.
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

// Fields
double Grid::get_field( char const field,
                        char const component,
                        std::size_t const x, std::size_t const y, std::size_t const z ) const {
    if ( field == 'e' ) {
        if ( component == 'x' ) return Ex_[idx(x,y,z)];
        else if ( component == 'y' ) return Ey_[idx(x,y,z)];
        else if ( component == 'z' ) return Ez_[idx(x,y,z)];
    } else if ( field == 'b' ) {
        if ( component == 'x' ) return Bx_[idx(x,y,z)];
        else if ( component == 'y' ) return By_[idx(x,y,z)];
        else if ( component == 'z' ) return Bz_[idx(x,y,z)];
    }
    throw std::invalid_argument{ "ERROR! CHECK PARAMETERS!" };
    return -1.0;
}

double Grid::field_mag( char const field,
                        std::size_t const x, std::size_t const y, std::size_t const z ) const {
    double Fx{ get_field( field, 'x', x, y, z ) };
    double Fy{ get_field( field, 'y', x, y, z ) };
    double Fz{ get_field( field, 'z', x, y, z ) };

    return std::sqrt( Fx*Fx + Fy*Fy + Fz*Fz );
}