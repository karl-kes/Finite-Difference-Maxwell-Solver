#pragma once

#include "../Config/config.hpp"

class Grid;

class Output {
private:
    std::string base_path_;

public:
    Output( std::string const &base_path = "Output" )
    : base_path_{ base_path }
    { }

    // Writes field data as binary:
    void write_field( Grid const &grid, double const time_step ) const;

    // Creates fresh directories by removing old data:
    void initialize() const {
        std::filesystem::remove_all( base_path_ );
        std::filesystem::create_directories( base_path_ + "/E" );
        std::filesystem::create_directories( base_path_ + "/B" );
    }

    // Generates the file for each time step:
    std::string file_name ( Field const field, std::size_t const time_step ) const {
        std::string prefix{ ( field == Field::ELECTRIC ) ? "/E/E" : "/B/B" };
        return base_path_ + prefix + std::to_string( time_step ) + ".bin";
    }
};