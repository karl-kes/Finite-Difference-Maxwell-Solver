#pragma once

#include "../../constant.hpp"

#include <memory>
#include <cmath>
#include <fstream>
#include <iostream>
#include <cctype>
#include <string>
#include <omp.h>
#include <vector>
#include <cstdint>
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

public:
    // Constructor:
    Grid( std::size_t new_Nx = 10, std::size_t new_Ny = 10, std::size_t new_Nz = 10,
          double new_dx = 5.0, double new_dy = 5.0, double new_dz = 5.0,
          double new_eps = 1.0, double new_mu = 1.0 );

    // System Simulation:
    void update_B();
    void update_E();
    void step();
    void hard_source_inject( std::size_t const x,
                             std::size_t const y,
                             std::size_t const z, 
                             double const value );
    void soft_source_inject();
    void dipole_antenna_inject();
    void gaussian_pulse_inject();
    void vector_volume( std::string const &file_name, char const field );

    // Getters:
    // Differentials
    double dx() const;
    double dy() const;
    double dz() const;
    double dt() const;
    // Speed, Mu, Epsilon
    double c() const;
    double c_sq() const;
    double eps() const;
    double mu() const;
    // Dimensions
    std::size_t Nx() const;
    std::size_t Ny() const;
    std::size_t Nz() const;
    // Fields
    double get_field( char const field,
                      char const component,
                      std::size_t const x,
                      std::size_t const y,
                      std::size_t const z ) const;
    double field_mag( char const field,
                      std::size_t const x,
                      std::size_t const y,
                      std::size_t const z ) const;

    // Helpers:
    void print_progress( int curr_time, int total_time ) const;
    // Finds 3D index
    std::size_t idx( std::size_t const x,
                     std::size_t const y,
                     std::size_t const z ) const;
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