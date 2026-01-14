#pragma once

#include <cstdint>

namespace config {
    static constexpr char B_field{ 'b' };
    static constexpr char E_field{ 'e' };
    static constexpr double PI{ 3.14159 };

    static constexpr std::size_t Nx{ 12 }, Ny{ 12 }, Nz{ 12 };
    static constexpr double dx{ 5.0 }, dy{ 5.0 }, dz{ 5.0 };
    static constexpr double eps{ 8.854e-12 }, mu{ 1.2566e-6 };

    static constexpr double cfl_factor{ 0.1 };
    static constexpr std::size_t total_time{ 1000 };
    static const std::size_t print_rate{ total_time / 100 };

    static constexpr double inject{ 100.0 };
    static constexpr double amp_one{ 100.0 };
    static constexpr double amp_two{ 100.0 };
    static constexpr double freq_one{ 1.0 };
    static constexpr double freq_two{ 1.0 };
}