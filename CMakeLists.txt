cmake_minimum_required(VERSION 3.20)
project(HPCSim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Build modes:
#
#   Debug:   -O0, warnings, debug symbols, tests.    For catching bugs.
#   Test:    -O3 -march=native etc., warnings, tests. Fast correctness check.
#   Release: -O3 -march=native etc., no warnings, no tests.  Ship/benchmark.

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(ALLOWED_BUILD_TYPES Debug Test Release)
list(FIND ALLOWED_BUILD_TYPES "${CMAKE_BUILD_TYPE}" _bt_index)
if(_bt_index EQUAL -1)
    message(FATAL_ERROR "Invalid CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}. Must be one of: Debug, Test, Release")
endif()

set(IS_DEBUG   FALSE)
set(IS_TEST    FALSE)
set(IS_RELEASE FALSE)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(IS_DEBUG TRUE)
elseif(CMAKE_BUILD_TYPE STREQUAL "Test")
    set(IS_TEST TRUE)
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(IS_RELEASE TRUE)
endif()

set(BUILD_TESTS FALSE)
if(IS_DEBUG OR IS_TEST)
    set(BUILD_TESTS TRUE)
endif()

# Sanitizer detection:
# ASan/UBSan are not available on MinGW. Only enable on Linux/macOS GCC/Clang.

set(HAS_SANITIZERS FALSE)
if(IS_DEBUG AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang" AND NOT MINGW)
    set(HAS_SANITIZERS TRUE)
endif()

# Core solver library (shared between main and tests):

file(GLOB_RECURSE LIB_SOURCES
    src/Classes/Grid/*.cpp
    src/Classes/Source/*.cpp
    src/Classes/PML/*.cpp
    src/Classes/Write_Output/*.cpp
    src/Classes/Validation/*.cpp
)

add_library(fdtd_core STATIC ${LIB_SOURCES})

target_include_directories(fdtd_core PUBLIC src)

# Compiler flags:

set(HPC_OPT_FLAGS
    -O3
    -march=native
    -mtune=native
    -funroll-loops
    -ffast-math
    -flto
    -DNDEBUG
)

set(WARN_FLAGS
    -Wall
    -Wextra
    -Wpedantic
    -Wshadow
    -Wconversion
)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")

    if(IS_DEBUG)
        target_compile_options(fdtd_core PRIVATE
            -O0 -g3
            ${WARN_FLAGS}
            -fno-omit-frame-pointer
        )
        if(HAS_SANITIZERS)
            target_compile_options(fdtd_core PRIVATE -fsanitize=address,undefined)
            target_link_options(fdtd_core PUBLIC -fsanitize=address,undefined)
        endif()

    elseif(IS_TEST)
        target_compile_options(fdtd_core PRIVATE
            ${HPC_OPT_FLAGS}
            ${WARN_FLAGS}
        )
        target_link_options(fdtd_core PUBLIC -flto)

    elseif(IS_RELEASE)
        target_compile_options(fdtd_core PRIVATE
            ${HPC_OPT_FLAGS}
        )
        target_link_options(fdtd_core PUBLIC -flto)
    endif()

elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")

    if(IS_DEBUG)
        target_compile_options(fdtd_core PRIVATE /Od /Zi /W4)
    elseif(IS_TEST)
        target_compile_options(fdtd_core PRIVATE /O2 /GL /DNDEBUG /W4)
        target_link_options(fdtd_core PUBLIC /LTCG)
    elseif(IS_RELEASE)
        target_compile_options(fdtd_core PRIVATE /O2 /GL /DNDEBUG)
        target_link_options(fdtd_core PUBLIC /LTCG)
    endif()

endif()

# OpenMP:

find_package(OpenMP REQUIRED)
target_link_libraries(fdtd_core PUBLIC OpenMP::OpenMP_CXX)

# IPO / Link-Time Optimization:

include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT ipo_output)
if(ipo_supported AND NOT IS_DEBUG)
    set_property(TARGET fdtd_core PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# Main simulation executable:

file(GLOB_RECURSE SIM_SOURCES
    src/Classes/Simulation/*.cpp
)

add_executable(main src/main.cpp ${SIM_SOURCES})
target_link_libraries(main PRIVATE fdtd_core)

if(ipo_supported AND NOT IS_DEBUG)
    set_property(TARGET main PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# Test suite:
#
#   Built in Debug and Test modes only. Excluded from Release.
#
#   Test:  -O2 -march=native + warnings (no -ffast-math so validation
#          floating-point comparisons aren't perturbed by reordering).
#
#   Debug: -O1 + warnings (+ sanitizers on Linux/macOS only).
#          -O1 not -O0 because integration tests step a 101^3 grid hundreds
#          of times — -O0 would take minutes.

enable_testing()

if(BUILD_TESTS)

    file(GLOB TEST_SOURCES tests/test_*.cpp)

    add_executable(run_tests ${TEST_SOURCES})

    target_include_directories(run_tests PRIVATE tests)
    target_link_libraries(run_tests PRIVATE fdtd_core)

    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        if(IS_DEBUG)
            target_compile_options(run_tests PRIVATE
                -O1 -g3
                ${WARN_FLAGS}
                -fno-omit-frame-pointer
            )
            if(HAS_SANITIZERS)
                target_compile_options(run_tests PRIVATE -fsanitize=address,undefined)
                target_link_options(run_tests PRIVATE -fsanitize=address,undefined)
            endif()
        elseif(IS_TEST)
            target_compile_options(run_tests PRIVATE
                -O2
                -march=native
                -DNDEBUG
                ${WARN_FLAGS}
            )
        endif()
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
        if(IS_DEBUG)
            target_compile_options(run_tests PRIVATE /Od /Zi /W4)
        elseif(IS_TEST)
            target_compile_options(run_tests PRIVATE /O2 /DNDEBUG /W4)
        endif()
    endif()

    if(ipo_supported AND IS_TEST)
        set_property(TARGET run_tests PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()

    # CTest targets — one per suite for granular reporting:
    add_test(NAME test_AlignedSoA       COMMAND run_tests AlignedSoA)
    add_test(NAME test_Grid             COMMAND run_tests Grid)
    add_test(NAME test_Source           COMMAND run_tests Source)
    add_test(NAME test_PML             COMMAND run_tests PML)
    add_test(NAME test_Integration      COMMAND run_tests Integration)
    add_test(NAME test_Output           COMMAND run_tests Output)
    add_test(NAME test_ValidationClass  COMMAND run_tests ValidationClass)
    add_test(NAME test_AllTests         COMMAND run_tests)

    if(IS_DEBUG)
        set_tests_properties(test_AlignedSoA PROPERTIES TIMEOUT 60)
        set_tests_properties(test_Grid       PROPERTIES TIMEOUT 60)
        set_tests_properties(test_Source     PROPERTIES TIMEOUT 60)
        set_tests_properties(test_Output     PROPERTIES TIMEOUT 60)
        set_tests_properties(test_PML             PROPERTIES TIMEOUT 300)
        set_tests_properties(test_Integration     PROPERTIES TIMEOUT 600)
        set_tests_properties(test_ValidationClass PROPERTIES TIMEOUT 300)
        set_tests_properties(test_AllTests        PROPERTIES TIMEOUT 600)
    else()
        set_tests_properties(test_AlignedSoA PROPERTIES TIMEOUT 30)
        set_tests_properties(test_Grid       PROPERTIES TIMEOUT 30)
        set_tests_properties(test_Source     PROPERTIES TIMEOUT 30)
        set_tests_properties(test_Output     PROPERTIES TIMEOUT 30)
        set_tests_properties(test_PML             PROPERTIES TIMEOUT 30)
        set_tests_properties(test_Integration     PROPERTIES TIMEOUT 300)
        set_tests_properties(test_ValidationClass PROPERTIES TIMEOUT 120)
        set_tests_properties(test_AllTests        PROPERTIES TIMEOUT 300)
    endif()
endif()

# Status:

message(STATUS "Build type:    ${CMAKE_BUILD_TYPE}")
message(STATUS "Compiler:      ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "IPO/LTO:       ${ipo_supported}")
message(STATUS "Tests:         ${BUILD_TESTS}")
message(STATUS "Sanitizers:    ${HAS_SANITIZERS}")

# Benchmark executable:

add_executable(benchmark benchmark/benchmark.cpp)
target_link_libraries(benchmark PRIVATE fdtd_core)

if(ipo_supported AND NOT IS_DEBUG)
    set_property(TARGET benchmark PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()