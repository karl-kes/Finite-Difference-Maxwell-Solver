#pragma once

#include <memory>
#include <vector>
#include <cstddef>

#include "../Config/config.hpp"
#include "../Config/memory.hpp"
#include "../PML/PML.hpp"

class Source;

class Grid {
private:
    // Number of field components:
    static constexpr std::size_t num_vec_components_{ 9 };
    static constexpr std::size_t alignment_bytes_{ SIMD_BYTES };
    static constexpr std::size_t doubles_per_alignment_{ SIMD_BYTES / sizeof( double ) };

    // Grid Dimensions:
    std::size_t Nx_, Ny_, Nz_;
    std::size_t Nx_padded_, Ny_padded_, Nz_padded_;
    
    // Total Grid Size:
    std::size_t N_;

    // Spatial Cell Size:
    double dx_, dy_, dz_;

    // Material Properties:
    double eps_, mu_;

    // Derived Constants:
    double c_, c_sq_, dt_;

    // Field Components:
    std::unique_ptr<double[], AlignedDeleter> memory_block_;
    std::size_t alignment_padding_;

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
    ~Grid();

    // System Simulation:
    void step();
    void apply_sources( double const time_step = 0.0 );
    void add_source( std::unique_ptr<Source> source );

    // Finds flattened 1D index given 3D space:
    [[nodiscard]] std::size_t idx(
        std::size_t const x,
        std::size_t const y,
        std::size_t const z ) const {
        return x + Nx_padded() * ( y + Ny_padded() * z );
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

    [[nodiscard]] std::size_t Nx_padded() const { return Nx_padded_; }
    [[nodiscard]] std::size_t Ny_padded() const { return Ny_padded_; }
    [[nodiscard]] std::size_t Nz_padded() const { return Nz_padded_; }
    [[nodiscard]] std::size_t alignment_padding() const { return alignment_padding_; }

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
    [[nodiscard]] double* Ex_ptr() { return memory_block_.get(); }
    [[nodiscard]] double* Ey_ptr() { return memory_block_.get() + alignment_padding(); }
    [[nodiscard]] double* Ez_ptr() { return memory_block_.get() + 2*alignment_padding(); }

    [[nodiscard]] double* Bx_ptr() { return memory_block_.get() + 3*alignment_padding(); }
    [[nodiscard]] double* By_ptr() { return memory_block_.get() + 4*alignment_padding(); }
    [[nodiscard]] double* Bz_ptr() { return memory_block_.get() + 5*alignment_padding(); }

    [[nodiscard]] double* Jx_ptr() { return memory_block_.get() + 6*alignment_padding(); }
    [[nodiscard]] double* Jy_ptr() { return memory_block_.get() + 7*alignment_padding(); }
    [[nodiscard]] double* Jz_ptr() { return memory_block_.get() + 8*alignment_padding(); }

    // Raw Pointer Field Components:
    // Immutable:
    [[nodiscard]] double const* Ex_ptr() const { return memory_block_.get(); }
    [[nodiscard]] double const* Ey_ptr() const { return memory_block_.get() + alignment_padding(); }
    [[nodiscard]] double const* Ez_ptr() const { return memory_block_.get() + 2*alignment_padding(); }

    [[nodiscard]] double const* Bx_ptr() const { return memory_block_.get() + 3*alignment_padding(); }
    [[nodiscard]] double const* By_ptr() const { return memory_block_.get() + 4*alignment_padding(); }
    [[nodiscard]] double const* Bz_ptr() const { return memory_block_.get() + 5*alignment_padding(); }

    [[nodiscard]] double const* Jx_ptr() const { return memory_block_.get() + 6*alignment_padding(); }
    [[nodiscard]] double const* Jy_ptr() const { return memory_block_.get() + 7*alignment_padding(); }
    [[nodiscard]] double const* Jz_ptr() const { return memory_block_.get() + 8*alignment_padding(); }

    // Diagnostics:
    [[nodiscard]] double total_energy() const;
    [[nodiscard]] double source_power() const;

    [[nodiscard]] PML const &pml() const { return pml_; }
};