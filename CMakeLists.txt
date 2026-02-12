cmake_minimum_required(VERSION 3.16)

project(HPCSim LANGUAGES CXX)

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optimization flags (equivalent to -O3 -march=native)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native")

file(GLOB_RECURSE SRC_FILES
    src/main.cpp
    src/Classes/Grid/*.cpp
    src/Classes/Source/*.cpp
    src/Classes/Write_Output/*.cpp
    src/Classes/PML/*.cpp
    src/Classes/Simulation/*.cpp
    src/Classes/Validation/*.cpp
)

# Executable
add_executable(main ${SRC_FILES})

# Include directories
target_include_directories(main PRIVATE src)

# OpenMP - equivalent to -fopenmp
find_package(OpenMP REQUIRED)

if(OpenMP_CXX_FOUND)
    target_link_libraries(main PRIVATE OpenMP::OpenMP_CXX)
endif()
