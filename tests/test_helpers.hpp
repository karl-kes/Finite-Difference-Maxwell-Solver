#pragma once

// Helper configs for testing. The default Simulation_Config uses constexpr members
// which can't be changed, so we provide test-specific configs by specializing
// through a wrapper that builds a config and overrides what we can.
//
// For tests needing smaller grids, we use the default config (100^3) but restrict
// our iteration ranges. For PML-off tests, we rely on runtime checks.
//
// The test helper struct provides commonly-needed derived values.

#include "../src/Classes/Config/config.hpp"
#include "../src/Classes/Grid/grid.hpp"

struct TestHelper {
    // Center of the default 101^3 grid:
    static constexpr std::size_t cx{ 50 };
    static constexpr std::size_t cy{ 50 };
    static constexpr std::size_t cz{ 50 };

    // Safe interior region (well away from PML):
    static constexpr std::size_t interior_lo{ 20 };
    static constexpr std::size_t interior_hi{ 80 };

    // Check if a point is in the interior (away from PML + boundary):
    static bool is_interior( std::size_t x, std::size_t y, std::size_t z ) {
        return x >= interior_lo && x <= interior_hi
            && y >= interior_lo && y <= interior_hi
            && z >= interior_lo && z <= interior_hi;
    }

    // Max absolute value of a field across the full grid:
    static double max_abs_field( Grid const& grid, double const* ptr ) {
        double mx{ 0.0 };
        for ( std::size_t z = 0; z < grid.Nz(); ++z )
            for ( std::size_t y = 0; y < grid.Ny(); ++y )
                for ( std::size_t x = 0; x < grid.Nx(); ++x ) {
                    double v = std::abs( ptr[grid.idx(x,y,z)] );
                    if ( v > mx ) mx = v;
                }
        return mx;
    }
};