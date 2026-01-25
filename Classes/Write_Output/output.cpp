#include "../Grid/grid.hpp"
#include "output.hpp"

void Output::initialize() const {
    std::filesystem::remove_all( base_path_ );
    std::filesystem::create_directories( base_path_ + "/E" );
    std::filesystem::create_directories( base_path_ + "/B" );
}

std::string Output::file_name ( Field field, std::size_t time_step ) const {
    std::string prefix{ ( field == Field::ELECTRIC ) ? "/E/E" : "/B/B" };
    return base_path_ + prefix + std::to_string( time_step ) + ".bin";
}

void Output::write_field( Grid const& grid, Field field, double time_step ) const {
    std::string path{ file_name( field, time_step ) };
    std::ofstream file( path, std::ios::binary | std::ios::out );
    
    if ( !file.is_open() ) {
        throw std::runtime_error{ "Failed to open file: " + path };
    }

    std::size_t nx{ grid.Nx() - 1 };
    std::size_t ny{ grid.Ny() - 1 };
    std::size_t nz{ grid.Nz() - 1 };

    uint64_t dimensions[3] = {
        static_cast<uint64_t>( nx ),
        static_cast<uint64_t>( ny ),
        static_cast<uint64_t>( nz )
    };
    file.write( reinterpret_cast<char const*>( dimensions ), sizeof( dimensions ) );

    Field field_type{ ( field == Field::ELECTRIC ) ? Field::ELECTRIC : Field::MAGNETIC };

    std::vector<double> buffer;
    buffer.reserve( nx * ny * 4 );

    for ( std::size_t z = 0; z < nz; ++z ) {
        buffer.clear();
        
        for ( std::size_t y = 0; y < ny; ++y ) {
            for ( std::size_t x = 0; x < nx; ++x ) {
                // Average to cell centers for visualization; deals with Yee staggering:
                double Fx{ 0.5 * ( grid.field( field_type, Component::X, x, y, z ) + 
                                   grid.field( field_type, Component::X, x+1, y, z ) ) };

                double Fy{ 0.5 * ( grid.field( field_type, Component::Y, x, y, z ) + 
                                   grid.field( field_type, Component::Y, x, y+1, z ) ) };

                double Fz{ 0.5 * ( grid.field( field_type, Component::Z, x, y, z ) + 
                                   grid.field( field_type, Component::Z, x, y, z+1 ) ) };

                double mag{ std::sqrt( Fx*Fx + Fy*Fy + Fz*Fz ) };

                buffer.push_back( Fx );
                buffer.push_back( Fy );
                buffer.push_back( Fz );
                buffer.push_back( mag );
            }
        }
        
        file.write( reinterpret_cast<char const*>( buffer.data() ), 
                    buffer.size() * sizeof( double ) );
    }
    file.close();
}