#pragma once

#include <memory>
#include <cmath>
#include <fstream>
#include <iostream>
#include <cctype>
#include <string>
#include <omp.h>

class Grid {
private:
    std::size_t Nx_, Ny_, Nz_;                  // Grid Size
    double dx_, dy_, dz_;                       // Spatial Differentials
    double c_;                                  // Speed
    double eps_;                                // Epsilon
    double mu_;                                 // Mu
    double dt_;                                 // Time Differential
    std::unique_ptr<double[]> Ex_, Ey_, Ez_;    // Electric Field
    std::unique_ptr<double[]> Bx_, By_, Bz_;    // Magnetic Field

public:
    // Constructor:
    Grid( std::size_t new_Nx = 10, std::size_t new_Ny = 10, std::size_t new_Nz = 10,
          double new_dx = 1.0, double new_dy = 1.0, double new_dz = 1.0,
          double new_c = 1.0, double new_eps = 1.0, double new_mu = 1.0 );

    // System Simulation:
    void update_B();
    void update_E();
    void step();
    void inject_source( std::size_t const x,
                        std::size_t const y,
                        std::size_t const z, 
                        double const value );
    void component_slice( std::size_t const z,
                          std::string const &file_name,
                          char const field,
                          char const component );
    void magnitude_slice( std::size_t const z,
                          std::string const &file_name,
                          char const field );

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

    // Helpers:
    // Finds 3D index
    std::size_t idx( std::size_t const x,
                     std::size_t const y,
                     std::size_t const z ) const;
    // Curls in X, Y, Z
    double curl_X( double const Y_0, double const Y_1,
                   double const Z_0, double const Z_1 ) const;

    double curl_Y( double const X_0, double const X_1,
                   double const Z_0, double const Z_1 ) const;

    double curl_Z( double const Y_0, double const Y_1,
                   double const X_0, double const X_1 ) const;

    double total_energy() const;
};