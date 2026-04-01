#pragma once

#include <memory>
#include <vector>
#include <cstddef>

#include "../Config/config.hpp"
#include "../../Utilities/aligned_soa.hpp"
#include "../PML/pml.hpp"

class Source;

// Monolithic field array indices:
enum GridArray : std::size_t {
    EX_, EY_, EZ_,
    HX_, HY_, HZ_,
    JX_, JY_, JZ_,
    EPS_X_, EPS_Y_, EPS_Z_,
    SIG_X_, SIG_Y_, SIG_Z_,
    MU_X_, MU_Y_,   MU_Z_,
    CA_X_, CA_Y_, CA_Z_, // E damp coeffs.
    CB_X_, CB_Y_, CB_Z_, // E curl coeffs.
    DB_X_, DB_Y_, DB_Z_, // H curl coeffs.
    NUM_ARRAYS_
};

class Grid {
private:
    // Grid Dimensions:
    std::size_t Nx_, Ny_, Nz_;
    std::size_t Nx_padded_, Ny_padded_, Nz_padded_;
    
    // Spatial Cell Size:
    double dx_, dy_, dz_;

    // Derived Constants:
    double dt_;

    // Field Components (single monolithic aligned block):
    AlignedSoA<double> data_;

    // Sources:
    std::vector<std::unique_ptr<Source>> sources_;

    // PML Absorbing Boundary:
    PML pml_;

    // Private Methods:
    [[nodiscard]] PML const &pml() const { return pml_; }

    
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

    [[nodiscard]] double dt() const { return dt_; }

    // Raw Pointer Getters (mutable):
    [[nodiscard]] double* Ex_ptr() { return data_[EX_]; }
    [[nodiscard]] double* Ey_ptr() { return data_[EY_]; }
    [[nodiscard]] double* Ez_ptr() { return data_[EZ_]; }

    [[nodiscard]] double* Hx_ptr() { return data_[HX_]; }
    [[nodiscard]] double* Hy_ptr() { return data_[HY_]; }
    [[nodiscard]] double* Hz_ptr() { return data_[HZ_]; }

    [[nodiscard]] double* Jx_ptr() { return data_[JX_]; }
    [[nodiscard]] double* Jy_ptr() { return data_[JY_]; }
    [[nodiscard]] double* Jz_ptr() { return data_[JZ_]; }

    [[nodiscard]] double* eps_x_ptr() { return data_[EPS_X_]; }
    [[nodiscard]] double* eps_y_ptr() { return data_[EPS_Y_]; }
    [[nodiscard]] double* eps_z_ptr() { return data_[EPS_Z_]; }

    [[nodiscard]] double* sig_x_ptr() { return data_[SIG_X_]; }
    [[nodiscard]] double* sig_y_ptr() { return data_[SIG_Y_]; }
    [[nodiscard]] double* sig_z_ptr() { return data_[SIG_Z_]; }

    [[nodiscard]] double* mu_x_ptr() { return data_[MU_X_]; }
    [[nodiscard]] double* mu_y_ptr() { return data_[MU_Y_]; }
    [[nodiscard]] double* mu_z_ptr() { return data_[MU_Z_]; }

    [[nodiscard]] double* Ca_x_ptr() { return data_[CA_X_]; }
    [[nodiscard]] double* Ca_y_ptr() { return data_[CA_Y_]; }
    [[nodiscard]] double* Ca_z_ptr() { return data_[CA_Z_]; }

    [[nodiscard]] double* Cb_x_ptr() { return data_[CB_X_]; }
    [[nodiscard]] double* Cb_y_ptr() { return data_[CB_Y_]; }
    [[nodiscard]] double* Cb_z_ptr() { return data_[CB_Z_]; }

    [[nodiscard]] double* Db_x_ptr() { return data_[DB_X_]; }
    [[nodiscard]] double* Db_y_ptr() { return data_[DB_Y_]; }
    [[nodiscard]] double* Db_z_ptr() { return data_[DB_Z_]; }

    // Raw Pointer Getters (immutable):
    [[nodiscard]] double const* Ex_ptr() const { return data_[EX_]; }
    [[nodiscard]] double const* Ey_ptr() const { return data_[EY_]; }
    [[nodiscard]] double const* Ez_ptr() const { return data_[EZ_]; }

    [[nodiscard]] double const* Hx_ptr() const { return data_[HX_]; }
    [[nodiscard]] double const* Hy_ptr() const { return data_[HY_]; }
    [[nodiscard]] double const* Hz_ptr() const { return data_[HZ_]; }

    [[nodiscard]] double const* Jx_ptr() const { return data_[JX_]; }
    [[nodiscard]] double const* Jy_ptr() const { return data_[JY_]; }
    [[nodiscard]] double const* Jz_ptr() const { return data_[JZ_]; }
    
    [[nodiscard]] double const* eps_x_ptr() const { return data_[EPS_X_]; }
    [[nodiscard]] double const* eps_y_ptr() const { return data_[EPS_Y_]; }
    [[nodiscard]] double const* eps_z_ptr() const { return data_[EPS_Z_]; }

    [[nodiscard]] double const* sig_x_ptr() const { return data_[SIG_X_]; }
    [[nodiscard]] double const* sig_y_ptr() const { return data_[SIG_Y_]; }
    [[nodiscard]] double const* sig_z_ptr() const { return data_[SIG_Z_]; }

    [[nodiscard]] double const* mu_x_ptr() const { return data_[MU_X_]; }
    [[nodiscard]] double const* mu_y_ptr() const { return data_[MU_Y_]; }
    [[nodiscard]] double const* mu_z_ptr() const { return data_[MU_Z_]; }

    [[nodiscard]] double const* Ca_x_ptr() const { return data_[CA_X_]; }
    [[nodiscard]] double const* Ca_y_ptr() const { return data_[CA_Y_]; }
    [[nodiscard]] double const* Ca_z_ptr() const { return data_[CA_Z_]; }

    [[nodiscard]] double const* Cb_x_ptr() const { return data_[CB_X_]; }
    [[nodiscard]] double const* Cb_y_ptr() const { return data_[CB_Y_]; }
    [[nodiscard]] double const* Cb_z_ptr() const { return data_[CB_Z_]; }

    [[nodiscard]] double const* Db_x_ptr() const { return data_[DB_X_]; }
    [[nodiscard]] double const* Db_y_ptr() const { return data_[DB_Y_]; }
    [[nodiscard]] double const* Db_z_ptr() const { return data_[DB_Z_]; }

    // Diagnostics:
    void update_H();
    void update_E();
    void bake_coefficients();

    [[nodiscard]] double h_energy() const;
    [[nodiscard]] double e_energy() const;
    [[nodiscard]] double total_energy() const;
};