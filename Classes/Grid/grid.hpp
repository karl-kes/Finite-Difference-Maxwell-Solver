#pragma once

#include <memory>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <omp.h>

#include "../Config/config.hpp"
#include "../Source/source.hpp"

class Grid {
private:
    // Grid Dimensions:
    std::size_t Nx_, Ny_, Nz_;

    // Spatial Cell Size:
    double dx_, dy_, dz_;

    // Material Properties:
    double eps_, mu_;

    // Derived Constants:
    double c_, dt_;

    // Field Components:
    std::unique_ptr<double[]> Ex_, Ey_, Ez_;
    std::unique_ptr<double[]> Bx_, By_, Bz_;
    std::unique_ptr<double[]> Jx_, Jy_, Jz_;

    // Sources:
    std::vector<std::unique_ptr<Source>> sources_;

    // Private Methods:
    // Finds flattened 1D index given 3D space:
    std::size_t idx( std::size_t x, std::size_t y, std::size_t z ) const;

    // Curl Calculation:
    double curl_x( double Y_0, double Y_1,
                   double Z_0, double Z_1) const;

    double curl_y( double X_0, double X_1,
                   double Z_0, double Z_1 ) const;

    double curl_z( double Y_0, double Y_1,
                   double X_0, double X_1 ) const;

    // Field Updates:
    void update_B();
    void update_E();

public:
    // Constructor:
    explicit Grid( Simulation_Config const &config );

    // Non-copyable and Non-movable:
    Grid ( Grid const& ) = delete;
    Grid &operator=( Grid const& ) = delete;
    Grid( Grid&& ) = default;
    Grid &operator=( Grid&& ) = default;

    // System Simulation:
    void step();
    void apply_sources( double time_step );
    void add_source( std::unique_ptr<Source> source );

    // Field Access:
    double field( char field, char component,
                  std::size_t x, std::size_t y, std::size_t z ) const;

    double &field( char field, char component,
                   std::size_t x, std::size_t y, std::size_t z );

    double field_magnitude( char field,
                            std::size_t x, std::size_t y, std::size_t z ) const;

    // Current Access:
    double &Jx( std::size_t x, std::size_t y, std::size_t z );
    double &Jy( std::size_t x, std::size_t y, std::size_t z );
    double &Jz( std::size_t x, std::size_t y, std::size_t z );

    // Getters:
    std::size_t Nx() const;
    std::size_t Ny() const;
    std::size_t Nz() const;

    double dx() const;
    double dy() const;
    double dz() const;

    double eps() const;
    double mu() const;

    double c() const;
    double c_sq() const;
    double dt() const;

    // Diagnostics:
    double total_energy() const;
};