#pragma once

#include <cstddef>
#include <cmath>
#include <numbers>

class Simulation_Config {
public:
    // Grid Dimensions:
    std::size_t static constexpr Nx{ 100 };
    std::size_t static constexpr Ny{ 100 };
    std::size_t static constexpr Nz{ 100 };
    std::size_t static constexpr size{ ( Nx + 1 ) * ( Ny + 1 ) * ( Nz + 1 ) };

    // Spatial Step Size:
    double static constexpr dx{ 5.0 };
    double static constexpr dy{ 5.0 };
    double static constexpr dz{ 5.0 };

    // Material Properties:
    double static constexpr mu{ 1.0 };
    double static constexpr eps{ 1.0 };
    double const c{ 1.0 / ( std::sqrt( mu * eps ) ) };

    // Time Stepping:
    double static constexpr cfl_factor{ 0.125 };
    std::size_t static constexpr total_time{ 1000 };
    double const dt{ cfl_factor / ( c * std::sqrt( 1.0/(dx*dx) + 1.0/(dy*dy) + 1.0/(dz*dz) ) ) };

    // Validation:
    bool static constexpr run_validation{ true };

    // PML (Perfectly Matched Layer):
    bool static constexpr use_pml{ true };
    std::size_t static constexpr pml_thickness{ 8 };
    int static constexpr pml_order{ 3 };
    double const pml_sigma_max{ 0.8 * ( pml_order + 1 ) / ( dx * std::sqrt( mu / eps ) ) };
    double static constexpr pml_kappa_max{ 15.0 };
    double static constexpr pml_alpha_max{ 0.05 };

    std::size_t output_interval() const { return total_time / 100; }
};

enum class Field {
    ELECTRIC,
    MAGNETIC
};

enum class Component {
    X,
    Y,
    Z
};