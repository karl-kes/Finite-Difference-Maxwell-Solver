#pragma once

#include "../Config/config.hpp"

#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

class Grid;

class Output {
private:
    std::size_t nx;
    std::size_t ny;
    std::size_t nz;

    std::size_t Sx;
    std::size_t Sy;
    std::size_t Sz;

    std::string base_path_;
    std::ofstream file_E_;
    std::ofstream file_H_;
    uint64_t dimensions_[3]{};

    std::size_t volume_size_{};
    std::vector<double> buffers_[2];
    int active_buf_{ 0 };

    std::thread writer_thread_;
    std::mutex mtx_;
    std::condition_variable cv_ready_;
    std::condition_variable cv_done_;
    bool has_work_{ false };
    bool shutdown_{ false };

    void writer_loop();
    void flush_buffer( double const* e_data, double const* h_data, std::size_t size );

public:
    Output( std::string const &base_path = "Output" )
    : base_path_{ base_path }
    { }

    ~Output() { finalize(); }

    Output( Output const& ) = delete;
    Output& operator=( Output const& ) = delete;

    void initialize( Grid const &grid );
    void write_field( Grid const &grid );
    void finalize();
};