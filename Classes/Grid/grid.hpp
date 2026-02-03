#pragma once

#include <memory>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <omp.h>

#include "../Config/config.hpp"
#include "../Source/source.hpp"
#include "../Write_Output/output.hpp"

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
    [[nodiscard]] std::size_t idx( std::size_t const x, std::size_t const y, std::size_t const z ) const { return x + Nx() * ( y + Ny() * z ); }

    // Curl Calculation:
    // Calculation flipped to fix Plotly left-hand rule default.
    [[nodiscard]] double curl_x(
                    double const Y_0, double const Y_1,
                    double const Z_0, double const Z_1) const
                    { return ( Z_1 - Z_0 ) / dy() - ( Y_1 - Y_0 ) / dz(); }

    [[nodiscard]] double curl_y(
                    double const X_0, double const X_1,
                    double const Z_0, double const Z_1 ) const
                    { return ( X_1 - X_0 ) / dz() - ( Z_1 - Z_0 ) / dx(); }

    [[nodiscard]] double curl_z(
                    double const Y_0, double const Y_1,
                    double const X_0, double const X_1 ) const 
                    { return ( Y_1 - Y_0 ) / dx() - ( X_1 - X_0 ) / dy(); }

    // Field Updates:
    void update_B();
    void update_E();

    void print_progress( double const current, double const total ) const {
        double const percent{ 100.0 * current / total };
        std::cout << "\rProgress: " << percent << "%" << std::flush;
    }

public:
    // Constructor:
    explicit Grid( Simulation_Config const &config );

    // Non-copyable and Non-movable:
    Grid ( Grid const& ) = delete;
    Grid &operator=( Grid const& ) = delete;
    Grid( Grid&& ) = default;
    Grid &operator=( Grid&& ) = default;

    // System Simulation:
    void step( Simulation_Config const &config, Output const &output, std::size_t const curr_time );
    void apply_sources( double const time_step = 0.0 );
    void add_source( std::unique_ptr<Source> source );

    // Field Access:
    [[nodiscard]] double field(
                    Field const field, Component const component,
                    std::size_t const x, std::size_t const y, std::size_t const z ) const;

    [[nodiscard]] double &field(
                    Field const field, Component const component,
                    std::size_t const x, std::size_t const y, std::size_t const z );

    [[nodiscard]] double field_magnitude(
                    Field const field,
                    std::size_t const x, std::size_t const y, std::size_t const z ) const;

    // Getters:
    // Dimensions:
    [[nodiscard]] std::size_t Nx() const { return Nx_; }
    [[nodiscard]] std::size_t Ny() const { return Ny_; }
    [[nodiscard]] std::size_t Nz() const { return Nz_; }

    // Grid Size:
    [[nodiscard]] double dx() const { return dx_; }
    [[nodiscard]] double dy() const { return dy_; }
    [[nodiscard]] double dz() const { return dz_; }

    // Wave Constants:
    [[nodiscard]] double eps() const { return eps_; }
    [[nodiscard]] double mu() const { return mu_; }
    [[nodiscard]] double c() const { return c_; }
    [[nodiscard]] double c_sq() const { return c_sq_; }

    // Time Step:
    [[nodiscard]] double dt() const { return dt_; }

    // Field Components:
    // Read Only:
    [[nodiscard]] const double &Ex( std::size_t const x, std::size_t const y, std::size_t const z ) const { return Ex_[idx(x,y,z)]; }
    [[nodiscard]] const double &Ey( std::size_t const x, std::size_t const y, std::size_t const z ) const { return Ey_[idx(x,y,z)]; }
    [[nodiscard]] const double &Ez( std::size_t const x, std::size_t const y, std::size_t const z ) const { return Ez_[idx(x,y,z)]; }

    [[nodiscard]] const double &Bx( std::size_t const x, std::size_t const y, std::size_t const z ) const { return Bx_[idx(x,y,z)]; }
    [[nodiscard]] const double &By( std::size_t const x, std::size_t const y, std::size_t const z ) const { return By_[idx(x,y,z)]; }
    [[nodiscard]] const double &Bz( std::size_t const x, std::size_t const y, std::size_t const z ) const { return Bz_[idx(x,y,z)]; }

    [[nodiscard]] const double &Jx( std::size_t const x, std::size_t const y, std::size_t const z ) const { return Jx_[idx(x,y,z)]; }
    [[nodiscard]] const double &Jy( std::size_t const x, std::size_t const y, std::size_t const z ) const { return Jy_[idx(x,y,z)]; }
    [[nodiscard]] const double &Jz( std::size_t const x, std::size_t const y, std::size_t const z ) const { return Jz_[idx(x,y,z)]; }

    // Writable:
    double &Ex( std::size_t const x, std::size_t const y, std::size_t const z ) { return Ex_[idx(x,y,z)]; }
    double &Ey( std::size_t const x, std::size_t const y, std::size_t const z ) { return Ey_[idx(x,y,z)]; }
    double &Ez( std::size_t const x, std::size_t const y, std::size_t const z ) { return Ez_[idx(x,y,z)]; }

    double &Bx( std::size_t const x, std::size_t const y, std::size_t const z ) { return Bx_[idx(x,y,z)]; }
    double &By( std::size_t const x, std::size_t const y, std::size_t const z ) { return By_[idx(x,y,z)]; }
    double &Bz( std::size_t const x, std::size_t const y, std::size_t const z ) { return Bz_[idx(x,y,z)]; }

    double &Jx( std::size_t const x, std::size_t const y, std::size_t const z ) { return Jx_[idx(x,y,z)]; }
    double &Jy( std::size_t const x, std::size_t const y, std::size_t const z ) { return Jy_[idx(x,y,z)]; }
    double &Jz( std::size_t const x, std::size_t const y, std::size_t const z ) { return Jz_[idx(x,y,z)]; }

    // Diagnostics:
    [[nodiscard]] double total_energy() const;
};