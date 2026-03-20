# Finite-Difference Maxwell Solver

A high-performance 3D Finite-Difference Time-Domain (FDTD) solver for Maxwell's equations, implemented in C++ with OpenMP parallelization, CPML absorbing boundaries, and interactive Plotly visualization.

![C++](https://img.shields.io/badge/C++-23-blue?logo=c%2B%2B)
![OpenMP](https://img.shields.io/badge/OpenMP-Parallel-green)
![Python](https://img.shields.io/badge/Python-3.x-yellow?logo=python)

## Demo

![FDTD Simulation](Assets/Fields.gif)

## Overview

This project implements the **Yee algorithm** for solving Maxwell's curl equations on a staggered grid. The FDTD method discretizes both space and time, updating electric and magnetic fields in a leapfrog scheme that naturally satisfies Faraday's and Ampère's laws.

### Maxwell's Equations

The solver implements the two curl equations with current density coupling:

$$\frac{\partial \mathbf{B}}{\partial t} = -\nabla \times \mathbf{E}$$

$$\frac{\partial \mathbf{E}}{\partial t} = c^2 \nabla \times \mathbf{B} - \frac{\mathbf{J}}{\varepsilon}$$

### CFL Stability Condition

The time step is automatically calculated to satisfy the Courant–Friedrichs–Lewy condition:

$$\Delta t \leq \frac{\alpha}{c \sqrt{\frac{1}{\Delta x^2} + \frac{1}{\Delta y^2} + \frac{1}{\Delta z^2}}}$$

where α (0 < α ≤ 1.0) is a configurable CFL factor for stability margin.

## Features

- **Yee-grid leapfrog scheme** with correct forward/backward difference staggering for E and B fields
- **CPML absorbing boundaries** (Convolutional Perfectly Matched Layer) with Roden–Gedney coefficients, polynomial-graded conductivity (σ), coordinate stretching (κ), and frequency shifting (α) for reflection-free domain truncation
- **Monolithic aligned memory layout** via `AlignedSoA<T>` — single contiguous allocation with SIMD-padded strides, `RESTRICT`-qualified pointers, and compile-time `ASSUME_ALIGNED` hints
- **Compile-time SIMD detection** (AVX-512 / AVX2 / SSE2 / scalar) for optimal alignment padding without wasting cache bandwidth
- **OpenMP parallelization** with `collapse(2)` on outer loops and `#pragma omp simd` on inner loops
- **Current density (J)** support with multiple source types: point sources, straight wire currents, and Gaussian pulses
- **Plane wave validation** with energy drift, phase correlation, and dispersion error metrics
- **Energy conservation tracking** and source power diagnostics
- **Binary I/O with slab-buffered writes** for efficient data export
- **Interactive 3D visualization** with synchronized E and B field volume rendering and vector cones via Plotly
- **External configuration file** (`config.cfg`) with runtime-parsed key-value pairs; unknown keys are warned, missing keys fall back to compiled defaults
- **Comprehensive test suite** (~47 test cases) across 7 suites: AlignedSoA, Grid, Source, PML, Integration, Output, and Validation — integrated with CTest
- **Three-tier CMake build system**: Debug (sanitizers + tests), Test (optimized + tests), Release (optimized, no tests)
- **Cross-platform**: Linux, macOS, Windows (MinGW/MSYS2); sanitizers enabled on Linux/macOS only

## Quick Start

### Prerequisites

- C++23 compiler (GCC recommended)
- CMake ≥ 3.20
- OpenMP support
- Python 3.x with NumPy and Plotly (for visualization)

### Build & Run

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j
./main

# Visualize results
python ../src/render.py
```

On Windows with MinGW:

```bash
mkdir build && cd build
cmake .. -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release
mingw32-make -j
./main.exe
python ../src/render.py
```

### Build All Modes

The `scripts/` directory contains convenience scripts that build all three modes (Debug, Test, Release), run tests in Debug and Test, and report results:

```bash
# Linux / macOS
./scripts/build.sh

# Windows (PowerShell)
.\scripts\build.ps1
```

Each script auto-detects CPU core count, produces `build_debug/`, `build_test/`, and `build_release/` directories, and runs the test suite in Debug and Test modes.

### Running Tests

Tests are built automatically in Debug and Test modes:

```bash
cmake .. -DCMAKE_BUILD_TYPE=Test
make -j
ctest --output-on-failure
```

Individual suites can be run directly:

```bash
./run_tests AlignedSoA
./run_tests PML
./run_tests Integration
```

## Configuration

Simulation parameters are set in `config.cfg` (key = value format, `#` comments). Any missing key uses the compiled default. Unknown keys produce a warning.

| Parameter | Description | Default |
|-----------|-------------|---------|
| `Nx, Ny, Nz` | Grid cell count per axis | 100 |
| `dx, dy, dz` | Spatial step size | 1.0 |
| `eps, mu` | Permittivity and permeability | 1.0 |
| `cfl_factor` | CFL stability factor (0 < α ≤ 1) | 1.0 |
| `total_steps` | Number of time steps | 1000 |
| `run_validation` | Run plane wave test on startup | true |
| `use_pml` | Enable CPML boundaries | true |
| `pml_order` | Polynomial grading order | 3 |
| `pml_kappa_max` | Maximum κ stretching | 15.0 |
| `pml_alpha_max` | Maximum α frequency shift | 0.3 |

Derived parameters computed automatically: `c` (wave speed), `dt` (CFL-limited time step), `pml_thickness` (scaled from grid volume), `pml_sigma_max` (optimal Berenger formula).

## Source Types

The simulator supports several excitation methods through a polymorphic `Source` interface:

- **`Point_Source`** — Isotropic hard source injection at a single grid point (all three J components)
- **`Straight_Wire_X`** — Sinusoidal current along the x-axis between two endpoints
- **`Gaussian_Pulse`** — Temporally Gaussian current pulse (Jz) at a single point

## CPML Absorbing Boundaries

The solver uses Convolutional Perfectly Matched Layers with Roden–Gedney coefficients to absorb outgoing waves at all six domain faces. The implementation features polynomial grading profiles for conductivity (σ), coordinate stretching (κ), and frequency shifting (α), with separate coefficients computed at E-field integer and B-field half-integer positions for proper Yee-grid alignment. Auxiliary convolution (ψ) arrays are maintained for all six face pairs across all three field directions.

The correct CPML formulation keeps the standard unmodified curl in the main E/B update; the ψ correction encodes the full stretching.

## Memory Layout

All field arrays (Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz) are stored in a single monolithic `AlignedSoA<double>` block. Each sub-array stride is rounded up to a SIMD-width boundary, ensuring every row starts at a naturally aligned address. SIMD width is detected at compile time (AVX-512 → 64 B, AVX2 → 32 B, SSE2 → 16 B, scalar → 8 B), avoiding the wasted cache bandwidth of a fixed 64-byte alignment on narrower hardware. Raw pointer access with `RESTRICT` qualification and `ASSUME_ALIGNED` hints gives the compiler full optimization freedom in the hot loops.

## Validation

The built-in plane wave validation test initializes a sinusoidal plane wave propagating in the x-direction, then tracks three metrics over a dynamically computed number of time steps:

- **Energy drift** — Total electromagnetic energy conservation (pass: <5%)
- **Phase correlation** — Pearson correlation between simulated and analytical Ey at a probe point (pass: >0.99)
- **Dispersion error** — Average phase error relative to expected propagation (pass: <10%)

Enable or disable with `run_validation` in `config.cfg`.

## Visualization

The Python renderer generates an interactive side-by-side animation showing E and B field evolution:

```bash
python src/render.py
```

Features include volume rendering for field intensity (Inferno colorscale), vector cones showing field direction and magnitude, and play/pause animation controls.

## Project Structure

```
├── config.cfg                              # Runtime simulation parameters
├── CMakeLists.txt                          # Three-tier build system
├── scripts/
│   ├── build.sh                            # Build all modes + run tests (Linux/macOS)
│   └── build.ps1                           # Build all modes + run tests (Windows/PowerShell)
├── src/
│   ├── main.cpp                            # Entry point: config → validation → simulation
│   ├── render.py                           # Plotly 3D visualization
│   ├── Utilities/
│   │   ├── AlignedSoA.hpp                  # SIMD-aligned monolithic SoA allocator
│   │   └── Macros.hpp                      # RESTRICT, ASSUME_ALIGNED portability macros
│   └── Classes/
│       ├── Config/
│       │   └── config.hpp                  # Simulation_Config with file parser & derived quantities
│       ├── Grid/
│       │   ├── grid.hpp                    # Grid class: field storage, indexing, diagnostics
│       │   └── grid.cpp                    # Yee update kernels (B, E), energy/power computation
│       ├── Source/
│       │   ├── source.hpp                  # Polymorphic Source base + Point/Wire/Gaussian
│       │   └── source.cpp                  # Source apply methods
│       ├── PML/
│       │   ├── PML.hpp                     # CPML class: coefficients, ψ arrays, indexing
│       │   └── PML.cpp                     # Roden–Gedney coefficient computation & ψ updates
│       ├── Simulation/
│       │   ├── Simulation.hpp              # Simulation driver declaration
│       │   └── Simulation.cpp              # Simulation loop, timing, progress output
│       ├── Write_Output/
│       │   ├── output.hpp                  # Binary output with Yee-averaged field export
│       │   └── output.cpp                  # Slab-buffered binary writer
│       └── Validation/
│           ├── Validation.hpp              # Plane wave test declaration
│           └── Validation.cpp              # Validation metrics & reporting
├── tests/
│   ├── test_framework.hpp                  # Zero-dependency header-only test framework
│   ├── test_helpers.hpp                    # Shared test utilities
│   ├── test_main.cpp                       # Test runner entry point
│   ├── test_aligned_soa.cpp                # AlignedSoA allocation & alignment tests
│   ├── test_grid.cpp                       # Grid construction, indexing, field access
│   ├── test_pml.cpp                        # CPML coefficient & boundary tests
│   ├── test_integration.cpp                # Full-simulation integration tests
│   └── test_output_validation.cpp          # Output I/O & validation class tests
├── output/                                 # Generated simulation data
│   ├── E/                                  # Electric field snapshots (.bin)
│   └── B/                                  # Magnetic field snapshots (.bin)
└── README.md
```

## Build Modes

| Mode | Optimization | Warnings | Tests | Sanitizers | Use Case |
|------|-------------|----------|-------|------------|----------|
| Debug | `-O0` (`-O1` for tests) | Full | Yes | ASan + UBSan (Linux/macOS) | Bug hunting |
| Test | `-O3 -march=native` | Full | Yes | No | Fast correctness verification |
| Release | `-O3 -march=native -ffast-math -flto` | No | No | No | Benchmarking / production |

## Performance Notes

- OpenMP `collapse(2)` parallelizes the two outer spatial loops; `#pragma omp simd` vectorizes the inner x-loop
- `RESTRICT`-qualified pointers on all hot-path field locals eliminate aliasing penalties
- Loop-invariant quantities (`inv_dx`, `dt_local`, stride constants) are hoisted before the parallel region
- `ASSUME_ALIGNED` macro carries zero runtime overhead — purely a compiler hint
- Compilation with `-O3 -march=native` is recommended; `-ffast-math` and LTO are enabled automatically in Release mode