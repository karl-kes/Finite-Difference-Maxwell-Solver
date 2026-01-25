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
    double c_, c_sq_, dt_;

    // Field Components:
    std::unique_ptr<double[]> Ex_, Ey_, Ez_;
    std::unique_ptr<double[]> Bx_, By_, Bz_;
    std::unique_ptr<double[]> Jx_, Jy_, Jz_;

    // Sources:
    std::vector<std::unique_ptr<Source>> sources_;

    // Private Methods:
    // Finds flattened 1D index given 3D space:
    [[nodiscard]] std::size_t idx( std::size_t x, std::size_t y, std::size_t z ) const;

    // Curl Calculation:
    // Calculation flipped to fix Plotly left-hand rule default.
    [[nodiscard]] double curl_x(
                    double Y_0, double Y_1,
                    double Z_0, double Z_1) const;

    [[nodiscard]] double curl_y(
                    double X_0, double X_1,
                    double Z_0, double Z_1 ) const;

    [[nodiscard]] double curl_z(
                    double Y_0, double Y_1,
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
    void apply_sources( double time_step = 0.0 );
    void add_source( std::unique_ptr<Source> source );

    // Field Access:
    [[nodiscard]] double field(
                    Field field, Component component,
                    std::size_t x, std::size_t y, std::size_t z ) const;

    [[nodiscard]] double &field(
                    Field field, Component component,
                    std::size_t x, std::size_t y, std::size_t z );

    [[nodiscard]] double field_magnitude(
                    Field field,
                    std::size_t x, std::size_t y, std::size_t z ) const;

    // Electric, Magnetic, Current Access:
    // Read only:
    [[nodiscard]] const double &Ex( std::size_t x, std::size_t y, std::size_t z ) const;
    [[nodiscard]] const double &Ey( std::size_t x, std::size_t y, std::size_t z ) const;
    [[nodiscard]] const double &Ez( std::size_t x, std::size_t y, std::size_t z ) const;

    [[nodiscard]] const double &Bx( std::size_t x, std::size_t y, std::size_t z ) const;
    [[nodiscard]] const double &By( std::size_t x, std::size_t y, std::size_t z ) const;
    [[nodiscard]] const double &Bz( std::size_t x, std::size_t y, std::size_t z ) const;

    [[nodiscard]] const double &Jx( std::size_t x, std::size_t y, std::size_t z ) const;
    [[nodiscard]] const double &Jy( std::size_t x, std::size_t y, std::size_t z ) const;
    [[nodiscard]] const double &Jz( std::size_t x, std::size_t y, std::size_t z ) const;

    // Writable:
    double &Ex( std::size_t x, std::size_t y, std::size_t z );
    double &Ey( std::size_t x, std::size_t y, std::size_t z );
    double &Ez( std::size_t x, std::size_t y, std::size_t z );

    double &Bx( std::size_t x, std::size_t y, std::size_t z );
    double &By( std::size_t x, std::size_t y, std::size_t z );
    double &Bz( std::size_t x, std::size_t y, std::size_t z );

    double &Jx( std::size_t x, std::size_t y, std::size_t z );
    double &Jy( std::size_t x, std::size_t y, std::size_t z );
    double &Jz( std::size_t x, std::size_t y, std::size_t z );

    // Getters:
    [[nodiscard]] std::size_t Nx() const;
    [[nodiscard]] std::size_t Ny() const;
    [[nodiscard]] std::size_t Nz() const;

    [[nodiscard]] double dx() const;
    [[nodiscard]] double dy() const;
    [[nodiscard]] double dz() const;

    [[nodiscard]] double eps() const;
    [[nodiscard]] double mu() const;

    [[nodiscard]] double c() const;
    [[nodiscard]] double c_sq() const;
    [[nodiscard]] double dt() const;

    // Diagnostics:
    [[nodiscard]] double total_energy() const;
};