cmake_minimum_required(VERSION 3.25)
project(kestrel LANGUAGES CXX)

if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
    message(FATAL_ERROR "Kestrel supports Linux and WSL2 only.")
endif()

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Toolchain requirements:

if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(FATAL_ERROR
        "Kestrel requires GCC 15 or newer (found ${CMAKE_CXX_COMPILER_ID}). "
        "Configure with -DCMAKE_CXX_COMPILER=g++-15.")
endif()

if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 15)
    message(FATAL_ERROR
        "Kestrel requires GCC 15 or newer (found ${CMAKE_CXX_COMPILER_VERSION}). "
        "Configure with -DCMAKE_CXX_COMPILER=g++-15.")
endif()

# CUDA backend:
option(KESTREL_ENABLE_CUDA "Enable CUDA when a CUDA compiler is available" ON)

set(KESTREL_CUDA FALSE)
if(KESTREL_ENABLE_CUDA)
    include(CheckLanguage)
    check_language(CUDA)

    if(CMAKE_CUDA_COMPILER)
        if(NOT CMAKE_CUDA_HOST_COMPILER)
            set(CMAKE_CUDA_HOST_COMPILER "${CMAKE_CXX_COMPILER}")
        endif()

        # CUDA 13.3 supports C++23 even when CMake has not learned its flag.
        set(CMAKE_CUDA23_STANDARD_COMPILE_OPTION "-std=c++23")
        set(CMAKE_CUDA23_EXTENSION_COMPILE_OPTION "-std=c++23")
        enable_language(CUDA)

        if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
            set(CMAKE_CUDA_ARCHITECTURES native)
        endif()

        set(KESTREL_CUDA TRUE)
    else()
        message(STATUS "Kestrel: no CUDA compiler found; using the CPU backend")
    endif()
endif()

# xpu:
include(FetchContent)

set(XPU_ENABLE_CUDA ${KESTREL_CUDA} CACHE BOOL "" FORCE)
set(XPU_BUILD_TESTS OFF CACHE BOOL "" FORCE)

fetchcontent_declare(
    xpu
    GIT_REPOSITORY https://github.com/UWHPC/xpu.git
    GIT_TAG main
    GIT_SHALLOW TRUE
)
fetchcontent_makeavailable(xpu)

function(kestrel_configure_sources)
    if(NOT KESTREL_CUDA)
        set_source_files_properties(
            ${ARGN}
            PROPERTIES
                LANGUAGE CXX
                COMPILE_OPTIONS "-x;c++"
        )
    endif()
endfunction()

# Build modes:

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

# Sanitizers are enabled for CPU Debug builds.
set(HAS_SANITIZERS FALSE)
if(IS_DEBUG AND NOT KESTREL_CUDA)
    set(HAS_SANITIZERS TRUE)
endif()

# Core solver library (shared between main and tests):

file(GLOB_RECURSE LIB_SOURCES
    src/classes/grid/*.cu
    src/classes/source/*.cu
    src/classes/pml/*.cu
    src/classes/write-output/*.cu
    src/classes/validation/*.cu
)
kestrel_configure_sources(${LIB_SOURCES})

add_library(fdtd-core STATIC ${LIB_SOURCES})

target_include_directories(fdtd-core PUBLIC src)
target_link_libraries(fdtd-core PUBLIC xpu::xpu)

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

if(IS_DEBUG)
    if(KESTREL_CUDA)
        target_compile_options(fdtd-core PRIVATE
            -O0 -g
            -Xcompiler=-Wall,-Wextra,-Wshadow,-Wconversion,-fno-omit-frame-pointer
        )
    else()
        target_compile_options(fdtd-core PRIVATE
            -O0 -g3
            ${WARN_FLAGS}
            -fno-omit-frame-pointer
            -fsanitize=address,undefined
        )
        target_link_options(fdtd-core PUBLIC -fsanitize=address,undefined)
    endif()
elseif(IS_TEST)
    if(KESTREL_CUDA)
        target_compile_options(fdtd-core PRIVATE
            -O3
            -Xcompiler=-O3,-march=native,-mtune=native,-funroll-loops,-ffast-math,-DNDEBUG
            -Xcompiler=-Wall,-Wextra,-Wshadow,-Wconversion
        )
    else()
        target_compile_options(fdtd-core PRIVATE
            ${HPC_OPT_FLAGS}
            ${WARN_FLAGS}
        )
        target_link_options(fdtd-core PUBLIC -flto)
    endif()
elseif(IS_RELEASE)
    if(KESTREL_CUDA)
        target_compile_options(fdtd-core PRIVATE
            -O3
            -Xcompiler=-O3,-march=native,-mtune=native,-funroll-loops,-ffast-math,-DNDEBUG
        )
    else()
        target_compile_options(fdtd-core PRIVATE
            ${HPC_OPT_FLAGS}
        )
        target_link_options(fdtd-core PUBLIC -flto)
    endif()
endif()

# OpenMP:

find_package(OpenMP REQUIRED COMPONENTS CXX)
target_link_libraries(fdtd-core PUBLIC OpenMP::OpenMP_CXX)
target_compile_options(fdtd-core PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-fopenmp>
)

# IPO / Link-Time Optimization:

include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT ipo_output LANGUAGES CXX)
if(ipo_supported AND NOT IS_DEBUG AND NOT KESTREL_CUDA)
    set_property(TARGET fdtd-core PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# Main simulation executable:

file(GLOB_RECURSE SIM_SOURCES
    src/classes/simulation/*.cu
)
kestrel_configure_sources(src/main.cu ${SIM_SOURCES})

add_executable(main src/main.cu ${SIM_SOURCES})
target_link_libraries(main PRIVATE fdtd-core)

if(ipo_supported AND NOT IS_DEBUG AND NOT KESTREL_CUDA)
    set_property(TARGET main PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# Test suite:

enable_testing()

if(BUILD_TESTS)

    file(GLOB TEST_SOURCES tests/test-*.cu)
    kestrel_configure_sources(${TEST_SOURCES})

    add_executable(run-tests ${TEST_SOURCES})

    target_include_directories(run-tests PRIVATE tests)
    target_link_libraries(run-tests PRIVATE fdtd-core)

    if(IS_DEBUG)
        if(KESTREL_CUDA)
            target_compile_options(run-tests PRIVATE
                -O1 -g
                -Xcompiler=-Wall,-Wextra,-Wshadow,-Wconversion,-fno-omit-frame-pointer
            )
        else()
            target_compile_options(run-tests PRIVATE
                -O1 -g3
                ${WARN_FLAGS}
                -fno-omit-frame-pointer
                -fsanitize=address,undefined
            )
            target_link_options(run-tests PRIVATE -fsanitize=address,undefined)
        endif()
    elseif(IS_TEST)
        if(KESTREL_CUDA)
            target_compile_options(run-tests PRIVATE
                -O2
                -Xcompiler=-O2,-march=native,-DNDEBUG
                -Xcompiler=-Wall,-Wextra,-Wshadow,-Wconversion
            )
        else()
            target_compile_options(run-tests PRIVATE
                -O2
                -march=native
                -DNDEBUG
                ${WARN_FLAGS}
            )
        endif()
    endif()

    if(ipo_supported AND IS_TEST AND NOT KESTREL_CUDA)
        set_property(TARGET run-tests PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()

    # CTest targets - one per suite for granular reporting:
    add_test(NAME aligned-soa      COMMAND run-tests AlignedSoA)
    add_test(NAME grid             COMMAND run-tests Grid)
    add_test(NAME source           COMMAND run-tests Source)
    add_test(NAME pml              COMMAND run-tests PML)
    add_test(NAME integration      COMMAND run-tests Integration)
    add_test(NAME output           COMMAND run-tests Output)
    add_test(NAME validation-class COMMAND run-tests ValidationClass)
    add_test(NAME all-tests        COMMAND run-tests)

    if(IS_DEBUG)
        set_tests_properties(aligned-soa      PROPERTIES TIMEOUT 60)
        set_tests_properties(grid             PROPERTIES TIMEOUT 60)
        set_tests_properties(source           PROPERTIES TIMEOUT 60)
        set_tests_properties(output           PROPERTIES TIMEOUT 60)
        set_tests_properties(pml              PROPERTIES TIMEOUT 300)
        set_tests_properties(integration      PROPERTIES TIMEOUT 600)
        set_tests_properties(validation-class PROPERTIES TIMEOUT 300)
        set_tests_properties(all-tests        PROPERTIES TIMEOUT 600)
    else()
        set_tests_properties(aligned-soa      PROPERTIES TIMEOUT 30)
        set_tests_properties(grid             PROPERTIES TIMEOUT 30)
        set_tests_properties(source           PROPERTIES TIMEOUT 30)
        set_tests_properties(output           PROPERTIES TIMEOUT 30)
        set_tests_properties(pml              PROPERTIES TIMEOUT 120)
        set_tests_properties(integration      PROPERTIES TIMEOUT 300)
        set_tests_properties(validation-class PROPERTIES TIMEOUT 120)
        set_tests_properties(all-tests        PROPERTIES TIMEOUT 300)
    endif()
endif()

# Status:

message(STATUS "Build type:    ${CMAKE_BUILD_TYPE}")
message(STATUS "Compiler:      ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "IPO/LTO:       ${ipo_supported}")
message(STATUS "Tests:         ${BUILD_TESTS}")
message(STATUS "Sanitizers:    ${HAS_SANITIZERS}")
if(KESTREL_CUDA)
    message(STATUS "xpu backend:   CUDA")
else()
    message(STATUS "xpu backend:   CPU")
endif()

# Benchmark executable (OpenMP-enabled, default):

kestrel_configure_sources(benchmark/benchmark.cu)

add_executable(benchmark benchmark/benchmark.cu)
target_link_libraries(benchmark PRIVATE fdtd-core)

if(ipo_supported AND NOT IS_DEBUG AND NOT KESTREL_CUDA)
    set_property(TARGET benchmark PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

add_library(fdtd-core-serial STATIC ${LIB_SOURCES})
target_include_directories(fdtd-core-serial PUBLIC src)
target_link_libraries(fdtd-core-serial PUBLIC xpu::xpu)

if(IS_DEBUG)
    if(KESTREL_CUDA)
        target_compile_options(fdtd-core-serial PRIVATE
            -O0 -g -Xcompiler=-Wall,-Wextra)
    else()
        target_compile_options(fdtd-core-serial PRIVATE
            -O0 -g -Wall -Wextra -Wpedantic)
    endif()
else()
    if(KESTREL_CUDA)
        target_compile_options(fdtd-core-serial PRIVATE
            -O3 -Xcompiler=-O3,-march=native,-DNDEBUG)
    else()
        target_compile_options(fdtd-core-serial PRIVATE
            -O3 -march=native -DNDEBUG -flto)
        target_link_options(fdtd-core-serial PUBLIC -flto)
    endif()
endif()

# Deliberately NO find_package(OpenMP) / no linkage to OpenMP::OpenMP_CXX here.

if(ipo_supported AND NOT IS_DEBUG AND NOT KESTREL_CUDA)
    set_property(TARGET fdtd-core-serial PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

add_executable(benchmark-serial benchmark/benchmark.cu)
target_link_libraries(benchmark-serial PRIVATE fdtd-core-serial)

if(ipo_supported AND NOT IS_DEBUG AND NOT KESTREL_CUDA)
    set_property(TARGET benchmark-serial PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
