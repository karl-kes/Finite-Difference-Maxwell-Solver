#pragma once

#include <memory>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <omp.h>

#include "../Config/config.hpp"
#include "../Source/source.hpp"
#include "../Write_Output/output.hpp"
#include "../PML/PML.hpp"

class Grid {
private:
    // Number of field components:
    std::size_t static constexpr num_vec_components_{ 9 };

    // Grid Dimensions:
    std::size_t Nx_, Ny_, Nz_;
    
    // Total Grid Size:
    std::size_t N_;

    // Spatial Cell Size:
    double dx_, dy_, dz_;

    // Material Properties:
    double eps_, mu_;

    // Derived Constants:
    double c_, c_sq_, dt_;

    // Field Components:
    std::unique_ptr<double[]> memory_block_;

    // Sources:
    std::vector<std::unique_ptr<Source>> sources_;

    // PML Absorbing Boundary:
    PML pml_;

    // Private Methods:
    // Field Updates:
    void update_B();
    void update_E();
    
public:
    // Constructor:
    explicit Grid( Simulation_Config const &config );

    // System Simulation:
    void step();
    void apply_sources( double const time_step = 0.0 );
    void add_source( std::unique_ptr<Source> source );

    // Finds flattened 1D index given 3D space:
    [[nodiscard]] std::size_t idx(
        std::size_t const x,
        std::size_t const y,
        std::size_t const z ) const {
        return x + Nx() * ( y + Ny() * z );
    }

    // Field Access:
    [[nodiscard]] double field(
        Field const field, Component const component,
        std::size_t const x, std::size_t const y, std::size_t const z
    ) const;

    [[nodiscard]] double &field(
        Field const field, Component const component,
        std::size_t const x, std::size_t const y, std::size_t const z
    );

    [[nodiscard]] double field_magnitude(
        Field const field,
        std::size_t const x, std::size_t const y, std::size_t const z
    ) const;

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

    // Raw Pointer Field Components:
    // Mutable:
    [[nodiscard]] double* Ex_() { return memory_block_.get(); }
    [[nodiscard]] double* Ey_() { return memory_block_.get() + N_; }
    [[nodiscard]] double* Ez_() { return memory_block_.get() + 2*N_; }

    [[nodiscard]] double* Bx_() { return memory_block_.get() + 3*N_; }
    [[nodiscard]] double* By_() { return memory_block_.get() + 4*N_; }
    [[nodiscard]] double* Bz_() { return memory_block_.get() + 5*N_; }

    [[nodiscard]] double* Jx_() { return memory_block_.get() + 6*N_; }
    [[nodiscard]] double* Jy_() { return memory_block_.get() + 7*N_; }
    [[nodiscard]] double* Jz_() { return memory_block_.get() + 8*N_; }

    // Raw Pointer Field Components:
    // Immutable:
    [[nodiscard]] double const* Ex_() const { return memory_block_.get(); }
    [[nodiscard]] double const* Ey_() const { return memory_block_.get() + N_; }
    [[nodiscard]] double const* Ez_() const { return memory_block_.get() + 2*N_; }

    [[nodiscard]] double const* Bx_() const { return memory_block_.get() + 3*N_; }
    [[nodiscard]] double const* By_() const { return memory_block_.get() + 4*N_; }
    [[nodiscard]] double const* Bz_() const { return memory_block_.get() + 5*N_; }

    [[nodiscard]] double const* Jx_() const { return memory_block_.get() + 6*N_; }
    [[nodiscard]] double const* Jy_() const { return memory_block_.get() + 7*N_; }
    [[nodiscard]] double const* Jz_() const { return memory_block_.get() + 8*N_; }

    // Diagnostics:
    [[nodiscard]] double total_energy() const;
    [[nodiscard]] double source_power() const;

    [[nodiscard]] PML const &pml() const { return pml_; }
};