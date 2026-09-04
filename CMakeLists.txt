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

# Project options:

option(KESTREL_ENABLE_CUDA "Enable CUDA when a CUDA compiler is available" ON)
option(KESTREL_BUILD_EXAMPLES "Build Kestrel examples" ON)

# CUDA backend:

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

# Compile .cu files as ordinary C++ when the CPU backend is selected.

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
    message(FATAL_ERROR
        "Invalid CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}. "
        "Must be one of: Debug, Test, Release")
endif()

set(IS_DEBUG FALSE)
set(IS_TEST FALSE)
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
    enable_testing()
endif()

set(HAS_SANITIZERS FALSE)
if(IS_DEBUG AND NOT KESTREL_CUDA)
    set(HAS_SANITIZERS TRUE)
endif()

include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT ipo_output LANGUAGES CXX)

if(KESTREL_BUILD_EXAMPLES)
    add_subdirectory(examples/maxwell)
endif()

# Status:

message(STATUS "Build type:    ${CMAKE_BUILD_TYPE}")
message(STATUS "Compiler:      ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "IPO/LTO:       ${ipo_supported}")
message(STATUS "Tests:         ${BUILD_TESTS}")
message(STATUS "Examples:      ${KESTREL_BUILD_EXAMPLES}")
message(STATUS "Sanitizers:    ${HAS_SANITIZERS}")
if(KESTREL_CUDA)
    message(STATUS "xpu backend:   CUDA")
else()
    message(STATUS "xpu backend:   CPU")
endif()
