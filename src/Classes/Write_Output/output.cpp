#include "output.hpp"

#include "../Grid/grid.hpp"

#include <fstream>
#include <stdexcept>
#include <cstdint>
#include <vector>
#include <cmath>

void Output::write_field( Grid const& grid, double const time_step ) const {
    // Electric:
    std::string path_E{ file_name( Field::ELECTRIC, time_step ) };
    std::ofstream file_E( path_E, std::ios::binary | std::ios::out );
    
    if ( !file_E.is_open() ) {
        throw std::runtime_error{ "Failed to open file: " + path_E };
    }

    std::size_t const nx{ grid.Nx() - 1 };
    std::size_t const ny{ grid.Ny() - 1 };
    std::size_t const nz{ grid.Nz() - 1 };

    uint64_t const dimensions[3] = {
        static_cast<uint64_t>( nx ),
        static_cast<uint64_t>( ny ),
        static_cast<uint64_t>( nz )
    };
    file_E.write( reinterpret_cast<char const*>( dimensions ), sizeof( dimensions ) );

    std::vector<double> buffer;
    buffer.reserve( nx * ny * 4 );

    double const* RESTRICT Ex{ grid.Ex_ptr() };
    double const* RESTRICT Ey{ grid.Ey_ptr() };
    double const* RESTRICT Ez{ grid.Ez_ptr() };

    double const* RESTRICT Bx{ grid.Bx_ptr() };
    double const* RESTRICT By{ grid.By_ptr() };
    double const* RESTRICT Bz{ grid.Bz_ptr() };

    for ( std::size_t z = 0; z < nz; ++z ) {
        buffer.clear();
        for ( std::size_t y = 0; y < ny; ++y ) {
            for ( std::size_t x = 0; x < nx; ++x ) {
                std::size_t const i{ grid.idx(x,y,z) };
                std::size_t const ix2{ grid.idx(x+1,y,z ) };
                std::size_t const iy2{ grid.idx(x,y+1,z) };
                std::size_t const iz2{ grid.idx(x,y,z+1) };

                double const Ex_mag{ 0.5 * ( Ex[i] + Ex[ix2] ) };
                double const Ey_mag{ 0.5 * ( Ey[i] + Ey[iy2] ) };
                double const Ez_mag{ 0.5 * ( Ez[i] + Ez[iz2] ) };

                double const E_mag{ std::sqrt( Ex_mag*Ex_mag + Ey_mag*Ey_mag + Ez_mag*Ez_mag ) };

                buffer.push_back( Ex_mag );
                buffer.push_back( Ey_mag );
                buffer.push_back( Ez_mag );
                buffer.push_back( E_mag );
            }
        }
        file_E.write( reinterpret_cast<char const*>( buffer.data() ), 
                      buffer.size() * sizeof( double ) );
    }
    file_E.close();

    // Magnetic:
    std::string path_B{ file_name( Field::MAGNETIC, time_step ) };
    std::ofstream file_B( path_B, std::ios::binary | std::ios::out );
    
    if ( !file_B.is_open() ) {
        throw std::runtime_error{ "Failed to open file: " + path_B };
    }

    file_B.write( reinterpret_cast<char const*>( dimensions ), sizeof( dimensions ) );

    for ( std::size_t z = 0; z < nz; ++z ) {
        buffer.clear();
        for ( std::size_t y = 0; y < ny; ++y ) {
            for ( std::size_t x = 0; x < nx; ++x ) {
                std::size_t const i{ grid.idx(x,y,z) };
                std::size_t const ix2{ grid.idx(x+1,y,z ) };
                std::size_t const iy2{ grid.idx(x,y+1,z) };
                std::size_t const iz2{ grid.idx(x,y,z+1) };
                
                double const Bx_mag{ 0.5 * ( Bx[i] + Bx[ix2] ) };
                double const By_mag{ 0.5 * ( By[i] + By[iy2] ) };
                double const Bz_mag{ 0.5 * ( Bz[i] + Bz[iz2] ) };

                double const B_mag{ std::sqrt( Bx_mag*Bx_mag + By_mag*By_mag + Bz_mag*Bz_mag ) };

                buffer.push_back( Bx_mag );
                buffer.push_back( By_mag );
                buffer.push_back( Bz_mag );
                buffer.push_back( B_mag );
            }
        }
        file_B.write( reinterpret_cast<char const*>( buffer.data() ), 
                      buffer.size() * sizeof( double ) );
    }
    file_B.close();
}