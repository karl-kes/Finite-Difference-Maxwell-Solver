#include "test_framework.hpp"
#include "../src/Utilities/aligned_soa.hpp"
#include <cstdint>

// ============================================================================
// AlignedSoA Unit Tests
// ============================================================================

TEST(AlignedSoA, RoundUpMultipleOfSimdWidth) {
    constexpr std::size_t elems_per_align = SIMD_BYTES / sizeof(double);

    ASSERT_EQ( AlignedSoA<double>::round_up( elems_per_align ), elems_per_align );
    ASSERT_EQ( AlignedSoA<double>::round_up( 2 * elems_per_align ), 2 * elems_per_align );
    ASSERT_EQ( AlignedSoA<double>::round_up( elems_per_align + 1 ), 2 * elems_per_align );
    ASSERT_EQ( AlignedSoA<double>::round_up( elems_per_align - 1 ), elems_per_align );
    ASSERT_EQ( AlignedSoA<double>::round_up( 0 ), std::size_t{0} );
    ASSERT_EQ( AlignedSoA<double>::round_up( 1 ), elems_per_align );
}

TEST(AlignedSoA, AllocationZeroInitialized) {
    AlignedSoA<double> soa{ 100, 4 };
    for ( std::size_t arr = 0; arr < 4; ++arr ) {
        for ( std::size_t i = 0; i < 100; ++i ) {
            ASSERT_EQ( soa[arr][i], 0.0 );
        }
    }
}

TEST(AlignedSoA, StrideIsAligned) {
    AlignedSoA<double> soa{ 101, 3 };
    std::size_t stride = soa.stride();
    constexpr std::size_t elems_per_align = SIMD_BYTES / sizeof(double);

    ASSERT_EQ( stride % elems_per_align, std::size_t{0} );
    ASSERT_TRUE( stride >= 101 );
}

TEST(AlignedSoA, SubArraysAreAligned) {
    AlignedSoA<double> soa{ 200, 5 };
    for ( std::size_t arr = 0; arr < 5; ++arr ) {
        auto ptr = reinterpret_cast<std::uintptr_t>( soa[arr] );
        ASSERT_EQ( ptr % SIMD_BYTES, std::uintptr_t{0} );
    }
}

TEST(AlignedSoA, SubArraysNonOverlapping) {
    AlignedSoA<double> soa{ 200, 5 };
    for ( std::size_t a = 0; a < 5; ++a ) {
        for ( std::size_t b = a + 1; b < 5; ++b ) {
            double* pa = soa[a];
            double* pb = soa[b];
            pa[0] = 1.0 + static_cast<double>(a);
            pb[0] = 1.0 + static_cast<double>(b);
            ASSERT_NEAR( pa[0], 1.0 + static_cast<double>(a), 1e-15 );
            ASSERT_NEAR( pb[0], 1.0 + static_cast<double>(b), 1e-15 );
        }
    }
}

TEST(AlignedSoA, MoveSemantics) {
    AlignedSoA<double> soa{ 50, 2 };
    soa[0][0] = 42.0;
    soa[1][0] = 99.0;

    AlignedSoA<double> moved{ std::move(soa) };
    ASSERT_NEAR( moved[0][0], 42.0, 1e-15 );
    ASSERT_NEAR( moved[1][0], 99.0, 1e-15 );

    AlignedSoA<double> assigned;
    assigned = std::move(moved);
    ASSERT_NEAR( assigned[0][0], 42.0, 1e-15 );
    ASSERT_NEAR( assigned[1][0], 99.0, 1e-15 );
}

TEST(AlignedSoA, NumElementsReported) {
    AlignedSoA<double> soa{ 137, 3 };
    ASSERT_EQ( soa.num_elements(), std::size_t{137} );
    ASSERT_TRUE( soa.stride() >= 137 );
}