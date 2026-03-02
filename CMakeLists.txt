cmake_minimum_required(VERSION 3.20)
project(HPCSim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

file(GLOB_RECURSE SRC_FILES
    src/main.cpp
    src/Classes/Grid/*.cpp
    src/Classes/Source/*.cpp
    src/Classes/Write_Output/*.cpp
    src/Classes/PML/*.cpp
    src/Classes/Simulation/*.cpp
    src/Classes/Validation/*.cpp
)

add_executable(main ${SRC_FILES})

target_include_directories(main PRIVATE src)

# Compiler detection
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")

    # Release: aggressive optimization
    target_compile_options(main PRIVATE
        $<$<CONFIG:Release>:
            -O3
            -march=native
            -mtune=native
            -funroll-loops
            -ffast-math
            -flto
            -DNDEBUG
        >
    )

    # Debug: sanitizers and diagnostics
    target_compile_options(main PRIVATE
        $<$<CONFIG:Debug>:
            -O0
            -g3
            -Wall
            -Wextra
            -Wpedantic
            -Wshadow
            -Wconversion
            -fsanitize=address,undefined
            -fno-omit-frame-pointer
        >
    )

    target_link_options(main PRIVATE
        $<$<CONFIG:Release>:-flto>
        $<$<CONFIG:Debug>:-fsanitize=address,undefined>
    )

elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")

    target_compile_options(main PRIVATE
        $<$<CONFIG:Release>:/O2 /GL /DNDEBUG>
        $<$<CONFIG:Debug>:/Od /Zi /W4 /fsanitize=address>
    )

    target_link_options(main PRIVATE
        $<$<CONFIG:Release>:/LTCG>
    )

endif()

# OpenMP
find_package(OpenMP REQUIRED)
target_link_libraries(main PRIVATE OpenMP::OpenMP_CXX)

# IPO / Link-Time Optimization (fallback if -flto flag isn't sufficient)
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT ipo_output)
if(ipo_supported)
    set_property(TARGET main PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# Architecture reporting
message(STATUS "Build type:    ${CMAKE_BUILD_TYPE}")
message(STATUS "Compiler:      ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "IPO/LTO:       ${ipo_supported}")