#pragma once

#include "../Config/config.hpp"
#include "../../Utilities/AlignedSoA.hpp"

#include <cmath>
#include <cstddef>
#include <algorithm>

enum CoeffArray : std::size_t {
    B_EX_,  C_EX_,  KAPPA_EX_,
    B_BX_,  C_BX_,  KAPPA_BX_,
    B_EY_,  C_EY_,  KAPPA_EY_,
    B_BY_,  C_BY_,  KAPPA_BY_,
    B_EZ_,  C_EZ_,  KAPPA_EZ_,
    B_BZ_,  C_BZ_,  KAPPA_BZ_,
    NUM_COEFF_ARRAYS_
};

enum PsiArray : std::size_t {
    PSI_EYX_,  PSI_EZX_,
    PSI_EXY_,  PSI_EZY_,
    PSI_EXZ_,  PSI_EYZ_,
    PSI_BYX_,  PSI_BZX_,
    PSI_BXY_,  PSI_BZY_,
    PSI_BXZ_,  PSI_BYZ_,
    NUM_PSI_ARRAYS_
};

class PML {
private:
    std::size_t thickness_;
    std::size_t Nx_, Ny_, Nz_;
    std::size_t Nx_padded_, Ny_padded_;
    int order_;
    double sigma_max_;
    double kappa_max_;
    double alpha_max_;

    AlignedSoA<double> coeffs_;
    AlignedSoA<double> psi_;

    std::size_t psi_face_x_;
    std::size_t psi_face_y_;
    std::size_t psi_face_z_;

    [[nodiscard]] double sigma( double const depth_norm ) const { return sigma_max_ * std::pow( depth_norm, order_ ); }
    [[nodiscard]] double kappa( double const depth_norm ) const { return 1.0 + ( kappa_max_ - 1.0 ) * std::pow( depth_norm, order_ ); }
    [[nodiscard]] double alpha( double const depth_norm ) const { return alpha_max_ * ( 1.0 - depth_norm ); }

    void compute_coefficients( 
        double const sigma_val, double const kappa_val, double const alpha_val,
        double const dt, double const eps,
        double &b_out, double &c_out
    ) const {
        // Roden-Gedney CPML coefficients:
        b_out = std::exp( -( sigma_val / kappa_val + alpha_val ) * dt / eps );
        double const denom{ sigma_val + kappa_val * alpha_val };
        c_out = ( denom > 1e-20 ) ? ( sigma_val / denom ) * ( b_out - 1.0 ) : 0.0;
    }

    [[nodiscard]] std::size_t psi_idx_x( std::size_t const d, std::size_t const y, std::size_t const z ) const { return d + thickness_ * ( y + Ny_ * z ); }
    [[nodiscard]] std::size_t psi_idx_y( std::size_t const x, std::size_t const d, std::size_t const z ) const { return x + Nx_ * ( d + thickness_ * z ); }
    [[nodiscard]] std::size_t psi_idx_z( std::size_t const x, std::size_t const y, std::size_t const d ) const { return x + Nx_ * ( y + Ny_ * d ); }

public:
    explicit PML( Simulation_Config const &config );

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

    [[nodiscard]] std::size_t thickness() const { return thickness_; }
    [[nodiscard]] bool is_active() const { return thickness_ > 0; }

    // Per-layer PML coefficients (length = thickness):
    [[nodiscard]] double* b_Ex_ptr() { return coeffs_[B_EX_]; }
    [[nodiscard]] double* c_Ex_ptr() { return coeffs_[C_EX_]; }
    [[nodiscard]] double* kappa_Ex_ptr() { return coeffs_[KAPPA_EX_]; }

    [[nodiscard]] double* b_Bx_ptr() { return coeffs_[B_BX_]; }
    [[nodiscard]] double* c_Bx_ptr() { return coeffs_[C_BX_]; }
    [[nodiscard]] double* kappa_Bx_ptr() { return coeffs_[KAPPA_BX_]; }

    [[nodiscard]] double* b_Ey_ptr() { return coeffs_[B_EY_]; }
    [[nodiscard]] double* c_Ey_ptr() { return coeffs_[C_EY_]; }
    [[nodiscard]] double* kappa_Ey_ptr() { return coeffs_[KAPPA_EY_]; }

    [[nodiscard]] double* b_By_ptr() { return coeffs_[B_BY_]; }
    [[nodiscard]] double* c_By_ptr() { return coeffs_[C_BY_]; }
    [[nodiscard]] double* kappa_By_ptr() { return coeffs_[KAPPA_BY_]; }

    [[nodiscard]] double* b_Ez_ptr() { return coeffs_[B_EZ_]; }
    [[nodiscard]] double* c_Ez_ptr() { return coeffs_[C_EZ_]; }
    [[nodiscard]] double* kappa_Ez_ptr() { return coeffs_[KAPPA_EZ_]; }

    [[nodiscard]] double* b_Bz_ptr() { return coeffs_[B_BZ_]; }
    [[nodiscard]] double* c_Bz_ptr() { return coeffs_[C_BZ_]; }
    [[nodiscard]] double* kappa_Bz_ptr() { return coeffs_[KAPPA_BZ_]; }

    [[nodiscard]] double const* b_Ex_ptr() const { return coeffs_[B_EX_]; }
    [[nodiscard]] double const* c_Ex_ptr() const { return coeffs_[C_EX_]; }
    [[nodiscard]] double const* kappa_Ex_ptr() const { return coeffs_[KAPPA_EX_]; }

    [[nodiscard]] double const* b_Bx_ptr() const { return coeffs_[B_BX_]; }
    [[nodiscard]] double const* c_Bx_ptr() const { return coeffs_[C_BX_]; }
    [[nodiscard]] double const* kappa_Bx_ptr() const { return coeffs_[KAPPA_BX_]; }

    [[nodiscard]] double const* b_Ey_ptr() const { return coeffs_[B_EY_]; }
    [[nodiscard]] double const* c_Ey_ptr() const { return coeffs_[C_EY_]; }
    [[nodiscard]] double const* kappa_Ey_ptr() const { return coeffs_[KAPPA_EY_]; }

    [[nodiscard]] double const* b_By_ptr() const { return coeffs_[B_BY_]; }
    [[nodiscard]] double const* c_By_ptr() const { return coeffs_[C_BY_]; }
    [[nodiscard]] double const* kappa_By_ptr() const { return coeffs_[KAPPA_BY_]; }

    [[nodiscard]] double const* b_Ez_ptr() const { return coeffs_[B_EZ_]; }
    [[nodiscard]] double const* c_Ez_ptr() const { return coeffs_[C_EZ_]; }
    [[nodiscard]] double const* kappa_Ez_ptr() const { return coeffs_[KAPPA_EZ_]; }

    [[nodiscard]] double const* b_Bz_ptr() const { return coeffs_[B_BZ_]; }
    [[nodiscard]] double const* c_Bz_ptr() const { return coeffs_[C_BZ_]; }
    [[nodiscard]] double const* kappa_Bz_ptr() const { return coeffs_[KAPPA_BZ_]; }

    [[nodiscard]] double* psi_Eyx_ptr() { return psi_[PSI_EYX_]; }
    [[nodiscard]] double* psi_Ezx_ptr() { return psi_[PSI_EZX_]; }
    [[nodiscard]] double* psi_Exy_ptr() { return psi_[PSI_EXY_]; }
    [[nodiscard]] double* psi_Ezy_ptr() { return psi_[PSI_EZY_]; }
    [[nodiscard]] double* psi_Exz_ptr() { return psi_[PSI_EXZ_]; }
    [[nodiscard]] double* psi_Eyz_ptr() { return psi_[PSI_EYZ_]; }
    
    [[nodiscard]] double* psi_Byx_ptr() { return psi_[PSI_BYX_]; }
    [[nodiscard]] double* psi_Bzx_ptr() { return psi_[PSI_BZX_]; }
    [[nodiscard]] double* psi_Bxy_ptr() { return psi_[PSI_BXY_]; }
    [[nodiscard]] double* psi_Bzy_ptr() { return psi_[PSI_BZY_]; }
    [[nodiscard]] double* psi_Bxz_ptr() { return psi_[PSI_BXZ_]; }
    [[nodiscard]] double* psi_Byz_ptr() { return psi_[PSI_BYZ_]; }

    [[nodiscard]] std::size_t psi_face_x() const { return psi_face_x_; }
    [[nodiscard]] std::size_t psi_face_y() const { return psi_face_y_; }
    [[nodiscard]] std::size_t psi_face_z() const { return psi_face_z_; }

    [[nodiscard]] std::size_t idx(
        std::size_t const x, std::size_t const y, std::size_t const z
    ) const {
        return x + Nx_padded_ * ( y + Ny_padded_ * z );
    }
};