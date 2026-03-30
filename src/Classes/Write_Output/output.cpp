#include "output.hpp"

#include "../Grid/grid.hpp"

#include <stdexcept>

void Output::initialize( Grid const &grid ) {
    std::filesystem::remove_all( base_path_ );
    std::filesystem::create_directories( base_path_ );

    dimensions_[0] = static_cast<uint64_t>( grid.Nx() - 1 );
    dimensions_[1] = static_cast<uint64_t>( grid.Ny() - 1 );
    dimensions_[2] = static_cast<uint64_t>( grid.Nz() - 1 );

    volume_size_ = static_cast<std::size_t>( dimensions_[0] )
                 * static_cast<std::size_t>( dimensions_[1] )
                 * static_cast<std::size_t>( dimensions_[2] ) * 3;

    buffers_[0].resize( volume_size_ * 2 );
    buffers_[1].resize( volume_size_ * 2 );
    active_buf_ = 0;

    nx = static_cast<std::size_t>( dimensions_[0] );
    ny = static_cast<std::size_t>( dimensions_[1] );
    nz = static_cast<std::size_t>( dimensions_[2] );

    Sx = 1;
    Sy = grid.Nx_padded();
    Sz = grid.Nx_padded() * grid.Ny_padded();

    std::string path_E{ base_path_ + "/E.bin" };
    file_E_.open( path_E, std::ios::binary | std::ios::out );
    if ( !file_E_.is_open() ) {
        throw std::runtime_error{ "Failed to open file: " + path_E };
    }
    file_E_.write( reinterpret_cast<char const*>( dimensions_ ), sizeof( dimensions_ ) );

    std::string path_B{ base_path_ + "/B.bin" };
    file_B_.open( path_B, std::ios::binary | std::ios::out );
    if ( !file_B_.is_open() ) {
        throw std::runtime_error{ "Failed to open file: " + path_B };
    }
    file_B_.write( reinterpret_cast<char const*>( dimensions_ ), sizeof( dimensions_ ) );

    shutdown_ = false;
    has_work_ = false;
    writer_thread_ = std::thread{ &Output::writer_loop, this };
}

void Output::writer_loop() {
    while ( true ) {
        std::unique_lock<std::mutex> lock{ mtx_ };
        cv_ready_.wait( lock, [this]{ return has_work_ || shutdown_; } );

        if ( shutdown_ && !has_work_ ) return;

        int flush_buf{ 1 - active_buf_ };
        double const* e_ptr{ buffers_[flush_buf].data() };
        double const* b_ptr{ buffers_[flush_buf].data() + volume_size_ };
        std::size_t size{ volume_size_ };

        lock.unlock();

        flush_buffer( e_ptr, b_ptr, size );

        {
            std::lock_guard<std::mutex> done_lock{ mtx_ };
            has_work_ = false;
        }
        cv_done_.notify_one();
    }
}

void Output::write_field( Grid const &grid ) {
    {
        std::unique_lock<std::mutex> lock{ mtx_ };
        cv_done_.wait( lock, [this]{ return !has_work_; } );
    }

    double* e_buf{ buffers_[active_buf_].data() };
    double* b_buf{ buffers_[active_buf_].data() + volume_size_ };

    double const* RESTRICT Ex{ grid.Ex_ptr() };
    double const* RESTRICT Ey{ grid.Ey_ptr() };
    double const* RESTRICT Ez{ grid.Ez_ptr() };

    std::size_t idx{ 0 };
    for ( std::size_t z = 0; z < nz; ++z ) {
        for ( std::size_t y = 0; y < ny; ++y ) {

            std::size_t const base{ y * Sy + z * Sz };
            for ( std::size_t x = 0; x < nx; ++x ) {
                std::size_t const i{ base + x };
                e_buf[idx    ] = 0.5 * ( Ex[i] + Ex[i + Sx] );
                e_buf[idx + 1] = 0.5 * ( Ey[i] + Ey[i + Sy] );
                e_buf[idx + 2] = 0.5 * ( Ez[i] + Ez[i + Sz] );
                idx += 3;
            }
        }
    }

    double const* RESTRICT Bx{ grid.Bx_ptr() };
    double const* RESTRICT By{ grid.By_ptr() };
    double const* RESTRICT Bz{ grid.Bz_ptr() };

    idx = 0;
    for ( std::size_t z = 0; z < nz; ++z ) {
        for ( std::size_t y = 0; y < ny; ++y ) {
            
            std::size_t const base{ y * Sy + z * Sz };
            for ( std::size_t x = 0; x < nx; ++x ) {
                std::size_t const i{ base + x };
                b_buf[idx    ] = 0.5 * ( Bx[i] + Bx[i + Sx] );
                b_buf[idx + 1] = 0.5 * ( By[i] + By[i + Sy] );
                b_buf[idx + 2] = 0.5 * ( Bz[i] + Bz[i + Sz] );
                idx += 3;
            }
        }
    }

    {
        std::lock_guard<std::mutex> lock{ mtx_ };
        active_buf_ = 1 - active_buf_;
        has_work_ = true;
    }
    cv_ready_.notify_one();
}

void Output::flush_buffer( double const* e_data, double const* b_data, std::size_t size ) {
    auto const bytes{ static_cast<std::streamsize>( size * sizeof( double ) ) };
    file_E_.write( reinterpret_cast<char const*>( e_data ), bytes );
    file_B_.write( reinterpret_cast<char const*>( b_data ), bytes );
}

void Output::finalize() {
    {
        std::unique_lock<std::mutex> lock{ mtx_ };
        cv_done_.wait( lock, [this]{ return !has_work_; } );
        shutdown_ = true;
    }
    cv_ready_.notify_one();
    if ( writer_thread_.joinable() ) writer_thread_.join();
    if ( file_E_.is_open() ) file_E_.close();
    if ( file_B_.is_open() ) file_B_.close();
}