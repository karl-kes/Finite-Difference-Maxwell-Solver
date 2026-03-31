#pragma once

#include "../Config/config.hpp"
#include "../../Utilities/aligned_soa.hpp"

#include <cmath>
#include <cstddef>
#include <algorithm>

enum CoeffArray : std::size_t {
    B_EX_,  C_EX_,  KAPPA_EX_,
    B_HX_,  C_HX_,  KAPPA_HX_,
    B_EY_,  C_EY_,  KAPPA_EY_,
    B_HY_,  C_HY_,  KAPPA_HY_,
    B_EZ_,  C_EZ_,  KAPPA_EZ_,
    B_HZ_,  C_HZ_,  KAPPA_HZ_,
    NUM_COEFF_ARRAYS_
};

enum PsiArray : std::size_t {
    PSI_EYX_,  PSI_EZX_,
    PSI_EXY_,  PSI_EZY_,
    PSI_EXZ_,  PSI_EYZ_,
    PSI_HYX_,  PSI_HZX_,
    PSI_HXY_,  PSI_HZY_,
    PSI_HXZ_,  PSI_HYZ_,
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

    void update_H_psi(
        double* RESTRICT Ex, double* RESTRICT Ey, double* RESTRICT Ez,
        double* RESTRICT Hx, double* RESTRICT Hy, double* RESTRICT Hz,
        double const dt, double const dx, double const dy, double const dz
    );

    void update_E_psi(
        double* RESTRICT Ex, double* RESTRICT Ey, double* RESTRICT Ez,
        double* RESTRICT Hx, double* RESTRICT Hy, double* RESTRICT Hz,
        double const dt, double const dx, double const dy, double const dz
    );

    [[nodiscard]] std::size_t thickness() const { return thickness_; }
    [[nodiscard]] bool is_active() const { return thickness_ > 0; }

    // Per-layer PML coefficients (length = thickness):
    [[nodiscard]] double* b_Ex_ptr() { return coeffs_[B_EX_]; }
    [[nodiscard]] double* c_Ex_ptr() { return coeffs_[C_EX_]; }
    [[nodiscard]] double* kappa_Ex_ptr() { return coeffs_[KAPPA_EX_]; }

    [[nodiscard]] double* b_Hx_ptr() { return coeffs_[B_HX_]; }
    [[nodiscard]] double* c_Hx_ptr() { return coeffs_[C_HX_]; }
    [[nodiscard]] double* kappa_Hx_ptr() { return coeffs_[KAPPA_HX_]; }

    [[nodiscard]] double* b_Ey_ptr() { return coeffs_[B_EY_]; }
    [[nodiscard]] double* c_Ey_ptr() { return coeffs_[C_EY_]; }
    [[nodiscard]] double* kappa_Ey_ptr() { return coeffs_[KAPPA_EY_]; }

    [[nodiscard]] double* b_Hy_ptr() { return coeffs_[B_HY_]; }
    [[nodiscard]] double* c_Hy_ptr() { return coeffs_[C_HY_]; }
    [[nodiscard]] double* kappa_Hy_ptr() { return coeffs_[KAPPA_HY_]; }

    [[nodiscard]] double* b_Ez_ptr() { return coeffs_[B_EZ_]; }
    [[nodiscard]] double* c_Ez_ptr() { return coeffs_[C_EZ_]; }
    [[nodiscard]] double* kappa_Ez_ptr() { return coeffs_[KAPPA_EZ_]; }

    [[nodiscard]] double* b_Hz_ptr() { return coeffs_[B_HZ_]; }
    [[nodiscard]] double* c_Hz_ptr() { return coeffs_[C_HZ_]; }
    [[nodiscard]] double* kappa_Hz_ptr() { return coeffs_[KAPPA_HZ_]; }

    [[nodiscard]] double const* b_Ex_ptr() const { return coeffs_[B_EX_]; }
    [[nodiscard]] double const* c_Ex_ptr() const { return coeffs_[C_EX_]; }
    [[nodiscard]] double const* kappa_Ex_ptr() const { return coeffs_[KAPPA_EX_]; }

    [[nodiscard]] double const* b_Hx_ptr() const { return coeffs_[B_HX_]; }
    [[nodiscard]] double const* c_Hx_ptr() const { return coeffs_[C_HX_]; }
    [[nodiscard]] double const* kappa_Hx_ptr() const { return coeffs_[KAPPA_HX_]; }

    [[nodiscard]] double const* b_Ey_ptr() const { return coeffs_[B_EY_]; }
    [[nodiscard]] double const* c_Ey_ptr() const { return coeffs_[C_EY_]; }
    [[nodiscard]] double const* kappa_Ey_ptr() const { return coeffs_[KAPPA_EY_]; }

    [[nodiscard]] double const* b_Hy_ptr() const { return coeffs_[B_HY_]; }
    [[nodiscard]] double const* c_Hy_ptr() const { return coeffs_[C_HY_]; }
    [[nodiscard]] double const* kappa_Hy_ptr() const { return coeffs_[KAPPA_HY_]; }

    [[nodiscard]] double const* b_Ez_ptr() const { return coeffs_[B_EZ_]; }
    [[nodiscard]] double const* c_Ez_ptr() const { return coeffs_[C_EZ_]; }
    [[nodiscard]] double const* kappa_Ez_ptr() const { return coeffs_[KAPPA_EZ_]; }

    [[nodiscard]] double const* b_Hz_ptr() const { return coeffs_[B_HZ_]; }
    [[nodiscard]] double const* c_Hz_ptr() const { return coeffs_[C_HZ_]; }
    [[nodiscard]] double const* kappa_Hz_ptr() const { return coeffs_[KAPPA_HZ_]; }

    [[nodiscard]] double* psi_Eyx_ptr() { return psi_[PSI_EYX_]; }
    [[nodiscard]] double* psi_Ezx_ptr() { return psi_[PSI_EZX_]; }
    [[nodiscard]] double* psi_Exy_ptr() { return psi_[PSI_EXY_]; }
    [[nodiscard]] double* psi_Ezy_ptr() { return psi_[PSI_EZY_]; }
    [[nodiscard]] double* psi_Exz_ptr() { return psi_[PSI_EXZ_]; }
    [[nodiscard]] double* psi_Eyz_ptr() { return psi_[PSI_EYZ_]; }
    
    [[nodiscard]] double* psi_Hyx_ptr() { return psi_[PSI_HYX_]; }
    [[nodiscard]] double* psi_Hzx_ptr() { return psi_[PSI_HZX_]; }
    [[nodiscard]] double* psi_Hxy_ptr() { return psi_[PSI_HXY_]; }
    [[nodiscard]] double* psi_Hzy_ptr() { return psi_[PSI_HZY_]; }
    [[nodiscard]] double* psi_Hxz_ptr() { return psi_[PSI_HXZ_]; }
    [[nodiscard]] double* psi_Hyz_ptr() { return psi_[PSI_HYZ_]; }

    [[nodiscard]] std::size_t psi_face_x() const { return psi_face_x_; }
    [[nodiscard]] std::size_t psi_face_y() const { return psi_face_y_; }
    [[nodiscard]] std::size_t psi_face_z() const { return psi_face_z_; }

    [[nodiscard]] std::size_t idx(
        std::size_t const x, std::size_t const y, std::size_t const z
    ) const {
        return x + Nx_padded_ * ( y + Ny_padded_ * z );
    }
};