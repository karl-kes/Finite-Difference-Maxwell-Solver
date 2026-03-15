#pragma once

#include <memory>
#include <vector>
#include <cstddef>

#include "../Config/config.hpp"
#include "../../Utilities/AlignedSoA.hpp"
#include "../PML/PML.hpp"

class Source;

// Monolithic field array indices:
enum GridArray : std::size_t {
    EX_,
    EY_,
    EZ_,
    BX_,
    BY_,
    BZ_,
    JX_,
    JY_,
    JZ_,
    NUM_GRID_ARRAYS_
};

class Grid {
private:
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

    // Field Components (single monolithic aligned block):
    AlignedSoA<double> fields_;

    // Sources:
    std::vector<std::unique_ptr<Source>> sources_;

    // PML Absorbing Boundary:
    PML pml_;

    // Private Methods:
    void update_B();
    void update_E();
    
public:
    // Constructor:
    explicit Grid( Simulation_Config const &config );
    ~Grid();

    // System Simulation:
    void step();
    void apply_sources( std::size_t const time_step = 0 );
    void add_source( std::unique_ptr<Source> source );

    // Finds flattened 1D index given 3D space:
    [[nodiscard]] std::size_t idx(
        std::size_t const x,
        std::size_t const y,
        std::size_t const z ) const {
        return x + Nx_padded_ * ( y + Ny_padded_ * z );
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
    [[nodiscard]] std::size_t Nx() const { return Nx_; }
    [[nodiscard]] std::size_t Ny() const { return Ny_; }
    [[nodiscard]] std::size_t Nz() const { return Nz_; }

    [[nodiscard]] std::size_t Nx_padded() const { return Nx_padded_; }
    [[nodiscard]] std::size_t Ny_padded() const { return Ny_padded_; }
    [[nodiscard]] std::size_t Nz_padded() const { return Nz_padded_; }

    [[nodiscard]] double dx() const { return dx_; }
    [[nodiscard]] double dy() const { return dy_; }
    [[nodiscard]] double dz() const { return dz_; }

    [[nodiscard]] double eps() const { return eps_; }
    [[nodiscard]] double mu() const { return mu_; }
    [[nodiscard]] double c() const { return c_; }
    [[nodiscard]] double c_sq() const { return c_sq_; }

    [[nodiscard]] double dt() const { return dt_; }

    // Raw Pointer Field Getters (mutable):
    [[nodiscard]] double* Ex_ptr() { return fields_[EX_]; }
    [[nodiscard]] double* Ey_ptr() { return fields_[EY_]; }
    [[nodiscard]] double* Ez_ptr() { return fields_[EZ_]; }

    [[nodiscard]] double* Bx_ptr() { return fields_[BX_]; }
    [[nodiscard]] double* By_ptr() { return fields_[BY_]; }
    [[nodiscard]] double* Bz_ptr() { return fields_[BZ_]; }

    [[nodiscard]] double* Jx_ptr() { return fields_[JX_]; }
    [[nodiscard]] double* Jy_ptr() { return fields_[JY_]; }
    [[nodiscard]] double* Jz_ptr() { return fields_[JZ_]; }

    // Raw Pointer Field Getters (immutable):
    [[nodiscard]] double const* Ex_ptr() const { return fields_[EX_]; }
    [[nodiscard]] double const* Ey_ptr() const { return fields_[EY_]; }
    [[nodiscard]] double const* Ez_ptr() const { return fields_[EZ_]; }

    [[nodiscard]] double const* Bx_ptr() const { return fields_[BX_]; }
    [[nodiscard]] double const* By_ptr() const { return fields_[BY_]; }
    [[nodiscard]] double const* Bz_ptr() const { return fields_[BZ_]; }

    [[nodiscard]] double const* Jx_ptr() const { return fields_[JX_]; }
    [[nodiscard]] double const* Jy_ptr() const { return fields_[JY_]; }
    [[nodiscard]] double const* Jz_ptr() const { return fields_[JZ_]; }

    // Diagnostics:
    [[nodiscard]] double total_energy() const;
    [[nodiscard]] double source_power() const;

    [[nodiscard]] PML const &pml() const { return pml_; }
};