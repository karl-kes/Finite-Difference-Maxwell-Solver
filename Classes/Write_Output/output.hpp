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

    // Creates fresh directories by removing old data:
    void initialize() const;

    // Writes field data as binary:
    void write_field( Grid const &grid, Field field, double time_step ) const;

    // Generates the file for each time step:
    std::string file_name( Field field, std::size_t time_step ) const;
};