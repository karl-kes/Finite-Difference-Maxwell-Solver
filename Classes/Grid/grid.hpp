#pragma once

#include "../../config.hpp"

#include <memory>
#include <cmath>
#include <fstream>
#include <iostream>
#include <cctype>
#include <string>
#include <omp.h>
#include <vector>
#include <filesystem>

class Grid {
private:
    std::size_t const Nx_, Ny_, Nz_;            // Grid Size
    double const dx_, dy_, dz_;                 // Spatial Differentials
    double const eps_;                          // Epsilon
    double const mu_;                           // Mu
    double const c_;                            // Speed
    double const dt_;                           // Time Differential
    std::unique_ptr<double[]> Ex_, Ey_, Ez_;    // Electric Field
    std::unique_ptr<double[]> Bx_, By_, Bz_;    // Magnetic Field
    std::unique_ptr<double[]> Jx_, Jy_, Jz_;    // Current

public:
    // Constructor:
    Grid( std::size_t new_Nx = 12, std::size_t new_Ny = 12, std::size_t new_Nz = 12,
          double new_dx = 5.0, double new_dy = 5.0, double new_dz = 5.0,
          double new_eps = 1.0, double new_mu = 1.0 );

    // System Simulation:
    void update_B();
    void update_E();
    void step();

    void straight_wire_x( double const amp, double const freq,
                          std::size_t const time,
                          std::size_t const y, std::size_t const z );

    void hard_source_inject( double const value,
                             std::size_t const x, std::size_t const y, std::size_t const z );

    void soft_source_inject( double const injection, std::size_t const x,
                             std::size_t const y, std::size_t const z );

    void dipole_antenna_inject( double const amp_one, double const amp_two,
                                double const freq_one, double const freq_two,
                                double const injection,
                                std::size_t const x, std::size_t const y, std::size_t const z );

    void gaussian_pulse_inject( double const injection,
                                std::size_t const x, std::size_t const y, std::size_t const z );

    void vector_volume( std::string const &file_name, char const field );

    // Getters:
    // Dimensions
    std::size_t Nx() const;
    std::size_t Ny() const;
    std::size_t Nz() const;

    // Differentials
    double dx() const;
    double dy() const;
    double dz() const;

    // Wave Constants
    double c() const;
    double c_sq() const;
    double eps() const;
    double mu() const;

    // Time Step
    double dt() const;

    // Fields
    double get_field( char const field,
                      char const component,
                      std::size_t const x,
                      std::size_t const y,
                      std::size_t const z ) const;

    double field_mag( char const field,
                      std::size_t const x, std::size_t const y, std::size_t const z ) const;

    // Helpers:
    void print_progress( int curr_time, int total_time ) const;

    // Finds 3D index
    std::size_t idx( std::size_t const x, std::size_t const y, std::size_t const z ) const;

    // Curls in X, Y, Z
    double curl_x( double const Y_0, double const Y_1,
                   double const Z_0, double const Z_1 ) const;

    double curl_y( double const X_0, double const X_1,
                   double const Z_0, double const Z_1 ) const;

    double curl_z( double const Y_0, double const Y_1,
                   double const X_0, double const X_1 ) const;

    // Total energy for validation
    double total_energy() const;

    // Deletes previous data and creates new files
    void create_directories() const;
};