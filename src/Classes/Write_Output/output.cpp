#include "../Grid/grid.hpp"
#include "output.hpp"

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

    for ( std::size_t z = 0; z < nz; ++z ) {
        buffer.clear();
        for ( std::size_t y = 0; y < ny; ++y ) {
            for ( std::size_t x = 0; x < nx; ++x ) {
                // Average to cell centers for visualization; deals with Yee staggering:
                double const Ex{ 0.5 * ( grid.Ex(x,y,z) + grid.Ex(x+1,y,z) ) };
                double const Ey{ 0.5 * ( grid.Ey(x,y,z) + grid.Ey(x,y+1,z) ) };
                double const Ez{ 0.5 * ( grid.Ez(x,y,z) + grid.Ez(x,y,z+1) ) };

                double const mag{ std::sqrt( Ex*Ex + Ey*Ey + Ez*Ez ) };

                buffer.push_back( Ex );
                buffer.push_back( Ey );
                buffer.push_back( Ez );
                buffer.push_back( mag );
            }
        }
        file_E.write( reinterpret_cast<char const*>( buffer.data() ), 
                      buffer.size() * sizeof( double ) );
    }
    file_E.close();

    // Magnetic:
    std::string path_B = file_name( Field::MAGNETIC, time_step );
    std::ofstream file_B( path_B, std::ios::binary | std::ios::out );
    
    if ( !file_B.is_open() ) {
        throw std::runtime_error{ "Failed to open file: " + path_B };
    }

    file_B.write( reinterpret_cast<char const*>( dimensions ), sizeof( dimensions ) );

    for ( std::size_t z = 0; z < nz; ++z ) {
        buffer.clear();
        for ( std::size_t y = 0; y < ny; ++y ) {
            for ( std::size_t x = 0; x < nx; ++x ) {
                // Average to cell centers for visualization; deals with Yee staggering:
                double const Bx{ 0.5 * ( grid.Bx(x,y,z) + grid.Bx(x+1,y,z) ) };
                double const By{ 0.5 * ( grid.By(x,y,z) + grid.By(x,y+1,z) ) };
                double const Bz{ 0.5 * ( grid.Bz(x,y,z) + grid.Bz(x,y,z+1) ) };

                double const mag{ std::sqrt( Bx*Bx + By*By + Bz*Bz ) };

                buffer.push_back( Bx );
                buffer.push_back( By );
                buffer.push_back( Bz );
                buffer.push_back( mag );
            }
        }
        file_B.write( reinterpret_cast<char const*>( buffer.data() ), 
                      buffer.size() * sizeof( double ) );
    }
    file_B.close();
}