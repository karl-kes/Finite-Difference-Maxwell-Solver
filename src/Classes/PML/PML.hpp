#pragma once

#include <memory>
#include <cmath>

#include "../Config/config.hpp"

class PML {
private:
    // PML Region Size:
    std::size_t thickness_;

    // Grid Dimensions:
    std::size_t Nx_, Ny_, Nz_;

    // Grading Parameters:
    int order_;
    double sigma_max_;
    double kappa_max_;
    double alpha_max_;

    // CPML Coefficients:
    // x-direction:
    std::unique_ptr<double[]> b_Ex_, c_Ex_;     // E-field half-grid
    std::unique_ptr<double[]> b_Bx_, c_Bx_;     // B-field half-grid
    std::unique_ptr<double[]> kappa_Ex_, kappa_Bx_;

    // y-direction:
    std::unique_ptr<double[]> b_Ey_, c_Ey_;
    std::unique_ptr<double[]> b_By_, c_By_;
    std::unique_ptr<double[]> kappa_Ey_, kappa_By_;

    // z-direction:
    std::unique_ptr<double[]> b_Ez_, c_Ez_;
    std::unique_ptr<double[]> b_Bz_, c_Bz_;
    std::unique_ptr<double[]> kappa_Ez_, kappa_Bz_;

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
    void compute_coefficients( double const sigma_val, double const kappa_val, double const alpha_val,
                               double const dt, double const eps,
                               double &b_out, double &c_out ) const {
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
    void update_B_psi( double *Ex, double *Ey, double *Ez,
                       double *Bx, double *By, double *Bz,
                       double const dt, double const dx, double const dy, double const dz );

    void update_E_psi( double *Ex, double *Ey, double *Ez,
                       double *Bx, double *By, double *Bz,
                       double const dt, double const dx, double const dy, double const dz,
                       double const c_sq );

    // Getters:
    [[nodiscard]] std::size_t thickness() const { return thickness_; }
    [[nodiscard]] bool is_active() const { return thickness_ > 0; }

    // Kappa accessors for modified curl in interior:
    [[nodiscard]] double const* kappa_Ex() const { return kappa_Ex_.get(); }
    [[nodiscard]] double const* kappa_Bx() const { return kappa_Bx_.get(); }
    [[nodiscard]] double const* kappa_Ey() const { return kappa_Ey_.get(); }
    [[nodiscard]] double const* kappa_By() const { return kappa_By_.get(); }
    [[nodiscard]] double const* kappa_Ez() const { return kappa_Ez_.get(); }
    [[nodiscard]] double const* kappa_Bz() const { return kappa_Bz_.get(); }

    // Grid-index helpers for the PML to use the same idx as Grid:
    [[nodiscard]] std::size_t idx( std::size_t const x, std::size_t const y, std::size_t const z ) const { return x + Nx_ * ( y + Ny_ * z ); }
};