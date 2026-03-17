#include "output.hpp"

#include "../Grid/grid.hpp"

#include <fstream>
#include <stdexcept>
#include <cstdint>
#include <vector>
#include <cmath>

void Output::write_field( Grid const& grid, std::size_t const time_step ) const {
    std::size_t const nx{ grid.Nx() - 1 };
    std::size_t const ny{ grid.Ny() - 1 };
    std::size_t const nz{ grid.Nz() - 1 };

    std::size_t const Sx{ 1 };
    std::size_t const Sy{ grid.Nx_padded() };
    std::size_t const Sz{ grid.Nx_padded() * grid.Ny_padded() };

    std::size_t const slab_size{ nx * ny * 4 };
    std::vector<double> buffer( slab_size );

    std::string path_E{ file_name( Field::ELECTRIC, time_step ) };
    std::ofstream file_E( path_E, std::ios::binary | std::ios::out );
    if ( !file_E.is_open() ) {
        throw std::runtime_error{ "Failed to open file: " + path_E };
    }

    uint64_t const dimensions[3] = {
        static_cast<uint64_t>( nx ),
        static_cast<uint64_t>( ny ),
        static_cast<uint64_t>( nz )
    };
    file_E.write( reinterpret_cast<char const*>( dimensions ), sizeof( dimensions ) );

    double const* RESTRICT Ex{ grid.Ex_ptr() };
    double const* RESTRICT Ey{ grid.Ey_ptr() };
    double const* RESTRICT Ez{ grid.Ez_ptr() };

    for ( std::size_t z = 0; z < nz; ++z ) {
        std::size_t buf_idx{ 0 };

        for ( std::size_t y = 0; y < ny; ++y ) {
            std::size_t const base{ y * Sy + z * Sz };

            for ( std::size_t x = 0; x < nx; ++x ) {
                std::size_t const i{ base + x };

                double const ex{ 0.5 * ( Ex[i] + Ex[i + Sx] ) };
                double const ey{ 0.5 * ( Ey[i] + Ey[i + Sy] ) };
                double const ez{ 0.5 * ( Ez[i] + Ez[i + Sz] ) };

                buffer[buf_idx    ] = ex;
                buffer[buf_idx + 1] = ey;
                buffer[buf_idx + 2] = ez;
                buffer[buf_idx + 3] = std::sqrt( ex * ex + ey * ey + ez * ez );
                buf_idx += 4;
            }
        }
        file_E.write( reinterpret_cast<char const*>( buffer.data() ),
                      static_cast<std::streamsize>( slab_size * sizeof( double ) ) );
    }
    file_E.close();

    std::string path_B{ file_name( Field::MAGNETIC, time_step ) };
    std::ofstream file_B( path_B, std::ios::binary | std::ios::out );
    if ( !file_B.is_open() ) {
        throw std::runtime_error{ "Failed to open file: " + path_B };
    }

    file_B.write( reinterpret_cast<char const*>( dimensions ), sizeof( dimensions ) );

    double const* RESTRICT Bx{ grid.Bx_ptr() };
    double const* RESTRICT By{ grid.By_ptr() };
    double const* RESTRICT Bz{ grid.Bz_ptr() };

    for ( std::size_t z = 0; z < nz; ++z ) {
        std::size_t buf_idx{ 0 };

        for ( std::size_t y = 0; y < ny; ++y ) {
            std::size_t const base{ y * Sy + z * Sz };
            
            for ( std::size_t x = 0; x < nx; ++x ) {
                std::size_t const i{ base + x };

                double const bx{ 0.5 * ( Bx[i] + Bx[i + Sx] ) };
                double const by{ 0.5 * ( By[i] + By[i + Sy] ) };
                double const bz{ 0.5 * ( Bz[i] + Bz[i + Sz] ) };

                buffer[buf_idx    ] = bx;
                buffer[buf_idx + 1] = by;
                buffer[buf_idx + 2] = bz;
                buffer[buf_idx + 3] = std::sqrt( bx * bx + by * by + bz * bz );
                buf_idx += 4;
            }
        }
        file_B.write( reinterpret_cast<char const*>( buffer.data() ),
                      static_cast<std::streamsize>( slab_size * sizeof( double ) ) );
    }
    file_B.close();
}