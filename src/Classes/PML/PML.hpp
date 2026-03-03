#pragma once

#include "../Config/config.hpp"
#include "../Config/memory.hpp"

#include <memory>
#include <cmath>
#include <cstddef>

#if defined(__GNUC__) || defined(__clang__)
    #define RESTRICT __restrict__
#elif defined(_MSC_VER)
    #define RESTRICT __restrict
#else
    #define RESTRICT
#endif

class PML {
private:
    // Alignment constant:
    static constexpr std::size_t num_coeff_arrays_{ 18 };
    static constexpr std::size_t alignment_bytes_{ SIMD_BYTES };
    static constexpr std::size_t doubles_per_alignment_{ SIMD_BYTES / sizeof(double) };

    // PML Region Size:
    std::size_t thickness_;
    std::size_t thickness_padded_;

    // Grid Dimensions:
    std::size_t Nx_, Ny_, Nz_;
    std::size_t Nx_padded_, Ny_padded_, Nz_padded_;

    // Grading Parameters:
    int order_;
    double sigma_max_;
    double kappa_max_;
    double alpha_max_;

    // CPML Coefficients:
    std::unique_ptr<double[], AlignedDeleter> coeff_memory_block_;

    // Psi (Auxiliary Convolution) Arrays:
    // E-field corrections:
    std::unique_ptr<double[]> psi_Eyx_, psi_Ezx_;  // x-faces
    std::unique_ptr<double[]> psi_Exy_, psi_Ezy_;  // y-faces
    std::unique_ptr<double[]> psi_Exz_, psi_Eyz_;  // z-faces

    // B-field corrections:
    std::unique_ptr<double[]> psi_Byx_, psi_Bzx_;  // x-faces
    std::unique_ptr<double[]> psi_Bxy_, psi_Bzy_;  // y-faces
    std::unique_ptr<double[]> psi_Bxz_, psi_Byz_;  // z-faces

    // Private Methods:
    // Polynomial grading profile (0 = interface, 1 = outer edge):
    [[nodiscard]] double sigma( double const depth_norm ) const { return sigma_max_ * std::pow( depth_norm, order_ ); }
    [[nodiscard]] double kappa( double const depth_norm ) const { return 1.0 + ( kappa_max_ - 1.0 ) * std::pow( depth_norm, order_ ); }
    [[nodiscard]] double alpha( double const depth_norm ) const { return alpha_max_ * ( 1.0 - depth_norm ); }

    // Compute b and c coefficients from grading values:
    void compute_coefficients( 
        double const sigma_val, double const kappa_val, double const alpha_val,
        double const dt, double const eps,
        double &b_out, double &c_out
    ) const {
        double const denom{ kappa_val * ( sigma_val + kappa_val * alpha_val ) };
        b_out = std::exp( -( sigma_val / kappa_val + alpha_val ) * dt / eps );
        c_out = ( denom > 1e-20 ) ? ( sigma_val / denom ) * ( b_out - 1.0 ) : 0.0;
    }

    // Flattened psi index for face arrays:
    // x-faces: pml_thickness * Ny * Nz (for each of the two x-faces)
    [[nodiscard]] std::size_t psi_idx_x( std::size_t const d, std::size_t const y, std::size_t const z ) const { return d + thickness_ * ( y + Ny_ * z ); }
    [[nodiscard]] std::size_t psi_idx_y( std::size_t const x, std::size_t const d, std::size_t const z ) const { return x + Nx_ * ( d + thickness_ * z ); }
    [[nodiscard]] std::size_t psi_idx_z( std::size_t const x, std::size_t const y, std::size_t const d ) const { return x + Nx_ * ( y + Ny_ * d ); }

    // Size of psi arrays for each face pair:
    [[nodiscard]] std::size_t psi_size_x() const { return thickness_ * Ny_ * Nz_; }
    [[nodiscard]] std::size_t psi_size_y() const { return Nx_ * thickness_ * Nz_; }
    [[nodiscard]] std::size_t psi_size_z() const { return Nx_ * Ny_ * thickness_; }

public:
    // Constructor:
    explicit PML( Simulation_Config const &config );

    // Apply PML corrections after standard update:
    void update_B_psi(
        double* RESTRICT Ex, double* RESTRICT Ey, double* RESTRICT Ez,
        double* RESTRICT Bx, double* RESTRICT By, double* RESTRICT Bz,
        double const dt, double const dx, double const dy, double const dz
    );

    void update_E_psi(
        double* RESTRICT Ex, double* RESTRICT Ey, double* RESTRICT Ez,
        double* RESTRICT Bx, double* RESTRICT By, double* RESTRICT Bz,
        double const dt, double const dx, double const dy, double const dz,
        double const c_sq
    );

    // Getters:
    [[nodiscard]] std::size_t thickness() const { return thickness_; }
    [[nodiscard]] std::size_t thickness_padded() const { return thickness_padded_; }
    [[nodiscard]] bool is_active() const { return thickness_ > 0; }

    // Raw Pointer Coefficient Access - Mutable
    [[nodiscard]] double* b_Ex_ptr() { return coeff_memory_block_.get(); }
    [[nodiscard]] double* c_Ex_ptr() { return coeff_memory_block_.get() + 1 * thickness_padded(); }
    [[nodiscard]] double* kappa_Ex_ptr() { return coeff_memory_block_.get() + 2 * thickness_padded(); }

    [[nodiscard]] double* b_Bx_ptr() { return coeff_memory_block_.get() + 3 * thickness_padded(); }
    [[nodiscard]] double* c_Bx_ptr() { return coeff_memory_block_.get() + 4 * thickness_padded(); }
    [[nodiscard]] double* kappa_Bx_ptr() { return coeff_memory_block_.get() + 5 * thickness_padded(); }

    [[nodiscard]] double* b_Ey_ptr() { return coeff_memory_block_.get() + 6 * thickness_padded(); }
    [[nodiscard]] double* c_Ey_ptr() { return coeff_memory_block_.get() + 7 * thickness_padded(); }
    [[nodiscard]] double* kappa_Ey_ptr() { return coeff_memory_block_.get() + 8 * thickness_padded(); }

    [[nodiscard]] double* b_By_ptr() { return coeff_memory_block_.get() + 9 * thickness_padded(); }
    [[nodiscard]] double* c_By_ptr() { return coeff_memory_block_.get() + 10 * thickness_padded(); }
    [[nodiscard]] double* kappa_By_ptr() { return coeff_memory_block_.get() + 11 * thickness_padded(); }

    [[nodiscard]] double* b_Ez_ptr() { return coeff_memory_block_.get() + 12 * thickness_padded(); }
    [[nodiscard]] double* c_Ez_ptr() { return coeff_memory_block_.get() + 13 * thickness_padded(); }
    [[nodiscard]] double* kappa_Ez_ptr() { return coeff_memory_block_.get() + 14 * thickness_padded(); }

    [[nodiscard]] double* b_Bz_ptr() { return coeff_memory_block_.get() + 15 * thickness_padded(); }
    [[nodiscard]] double* c_Bz_ptr() { return coeff_memory_block_.get() + 16 * thickness_padded(); }
    [[nodiscard]] double* kappa_Bz_ptr() { return coeff_memory_block_.get() + 17 * thickness_padded(); }

    // Raw Pointer Coefficient Access - Immutable:
    [[nodiscard]] double const* b_Ex_ptr() const { return coeff_memory_block_.get(); }
    [[nodiscard]] double const* c_Ex_ptr() const { return coeff_memory_block_.get() + 1 * thickness_padded(); }
    [[nodiscard]] double const* kappa_Ex_ptr() const { return coeff_memory_block_.get() + 2 * thickness_padded(); }
    
    [[nodiscard]] double const* b_Bx_ptr() const { return coeff_memory_block_.get() + 3 * thickness_padded(); }
    [[nodiscard]] double const* c_Bx_ptr() const { return coeff_memory_block_.get() + 4 * thickness_padded(); }
    [[nodiscard]] double const* kappa_Bx_ptr() const { return coeff_memory_block_.get() + 5 * thickness_padded(); }

    [[nodiscard]] double const* b_Ey_ptr() const { return coeff_memory_block_.get() + 6 * thickness_padded(); }
    [[nodiscard]] double const* c_Ey_ptr() const { return coeff_memory_block_.get() + 7 * thickness_padded(); }
    [[nodiscard]] double const* kappa_Ey_ptr() const { return coeff_memory_block_.get() + 8 * thickness_padded(); }

    [[nodiscard]] double const* b_By_ptr() const { return coeff_memory_block_.get() + 9 * thickness_padded(); }
    [[nodiscard]] double const* c_By_ptr() const { return coeff_memory_block_.get() + 10 * thickness_padded(); }
    [[nodiscard]] double const* kappa_By_ptr() const { return coeff_memory_block_.get() + 11 * thickness_padded(); }

    [[nodiscard]] double const* b_Ez_ptr() const { return coeff_memory_block_.get() + 12 * thickness_padded(); }
    [[nodiscard]] double const* c_Ez_ptr() const { return coeff_memory_block_.get() + 13 * thickness_padded(); }
    [[nodiscard]] double const* kappa_Ez_ptr() const { return coeff_memory_block_.get() + 14 * thickness_padded(); }

    [[nodiscard]] double const* b_Bz_ptr() const { return coeff_memory_block_.get() + 15 * thickness_padded(); }
    [[nodiscard]] double const* c_Bz_ptr() const { return coeff_memory_block_.get() + 16 * thickness_padded(); }
    [[nodiscard]] double const* kappa_Bz_ptr() const { return coeff_memory_block_.get() + 17 * thickness_padded(); }

    // Grid-index helpers for the PML to use the same idx as Grid:
    [[nodiscard]] std::size_t idx(
        std::size_t const x,
        std::size_t const y,
        std::size_t const z
    ) const{
        return x + Nx_padded_ * ( y + Ny_padded_ * z );
    }
};