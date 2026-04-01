# Finite-Difference Maxwell Solver

A high-performance 3D Finite-Difference Time-Domain (FDTD) solver for Maxwell's equations, implemented in C++ with OpenMP parallelization, CPML absorbing boundaries, and Rerun-based 3D visualization.

![C++](https://img.shields.io/badge/C++-23-blue?logo=c%2B%2B)
![OpenMP](https://img.shields.io/badge/OpenMP-Parallel-green)
![Python](https://img.shields.io/badge/Python-3.x-yellow?logo=python)

## Demo

| Combined E & H Fields | H-Field | E-Field |
|:---:|:---:|:---:|
| ![EM Fields](Assets/EM_Fields.gif) | ![H Field](Assets/H_Field.gif) | ![E Field](Assets/E_Field.gif) |

*AC current loop radiating electromagnetic waves, visualized with Rerun. E-field rendered in viridis, H-field in inferno.*

## Overview

This project implements the **Yee algorithm** for solving Maxwell's curl equations on a staggered grid. The FDTD method discretizes both space and time, updating electric and magnetic fields in a leapfrog scheme that naturally satisfies Faraday's and Ampère's laws.

### Maxwell's Equations

The solver implements the two curl equations in H-field formulation with per-cell material coefficients:

$$\frac{\partial \mathbf{H}}{\partial t} = -\frac{1}{\mu}\nabla \times \mathbf{E}$$

$$\frac{\partial \mathbf{E}}{\partial t} = \frac{1}{\varepsilon}\nabla \times \mathbf{H} - \frac{\sigma}{\varepsilon}\mathbf{E} - \frac{\mathbf{J}}{\varepsilon}$$

The E-field update uses precomputed loss-incorporating coefficients:

$$E_i^{n+1} = C_a \cdot E_i^n + C_b \cdot \left(\nabla \times \mathbf{H} - \mathbf{J}\right)_i$$

where $C_a = \frac{1 - \sigma\Delta t / 2\varepsilon}{1 + \sigma\Delta t / 2\varepsilon}$ and $C_b = \frac{1}{\varepsilon/\Delta t + \sigma/2}$.

The H-field update uses: $H_i^{n+1} = H_i^n - D_b \cdot (\nabla \times \mathbf{E})_i$ where $D_b = \Delta t / \mu$.

### CFL Stability Condition

The time step is automatically calculated to satisfy the Courant–Friedrichs–Lewy (CFL) condition:

$$\Delta t \leq \frac{\alpha}{c \sqrt{\frac{1}{\Delta x^2} + \frac{1}{\Delta y^2} + \frac{1}{\Delta z^2}}}$$

where α (0 < α ≤ 1.0) is a configurable CFL factor for stability margin.

## Features

- **Yee-grid leapfrog scheme** with correct forward/backward difference staggering for E and H fields in H-field formulation
- **Per-cell material arrays** (ε, μ, σ per component) with precomputed update coefficients (Ca, Cb for E; Db for H), baked at initialization for zero per-step overhead
- **CPML absorbing boundaries** (Convolutional Perfectly Matched Layer) with Roden–Gedney coefficients, polynomial-graded conductivity (σ), coordinate stretching (κ), and frequency shifting (α) for reflection-free domain truncation
- **Monolithic aligned memory layout** via `AlignedSoA<T>` — single contiguous allocation with SIMD-padded strides, `RESTRICT`-qualified pointers, and compile-time `ASSUME_ALIGNED` hints
- **Compile-time SIMD detection** (AVX-512 / AVX2 / SSE2 / scalar) for optimal alignment padding without wasting cache bandwidth
- **OpenMP parallelization** with `collapse(2)` on outer loops and `#pragma omp simd` on inner loops
- **Current density (J)** support with multiple source types: AC current loop, straight wire, point source, and Gaussian pulse
- **Leapfrog-synchronized energy diagnostics** with Yee-averaged `e_energy()` and `h_energy()` methods and time-interpolated total energy tracking in the renderer
- **Plane wave validation** with energy drift, phase correlation, and dispersion error metrics
- **Binary I/O with double-buffered writes** — persistent writer thread with condition-variable synchronization for overlap of computation and disk I/O
- **Rerun 3D visualization** with Points3D volume rendering and Arrows3D vector cones, separate viridis (E-field) and inferno (H-field) colormaps, and time-series energy plots
- **External configuration file** (`config.cfg`) with runtime-parsed key-value pairs; unknown keys are warned, missing keys fall back to compiled defaults
- **Comprehensive test suite** Across 7 suites: AlignedSoA, Grid, Source, PML, Integration, Output, and Validation

## Quick Start

### Prerequisites

- C++23 compiler (GCC recommended)
- CMake ≥ 3.20
- OpenMP support
- Python 3.x with NumPy and Rerun (for visualization)

### Build & Run

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j
./main

# Visualize results
python ../render.py
```

On Windows with MinGW:

```bash
mkdir build && cd build
cmake .. -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release
mingw32-make -j
./main.exe
python ../render.py
```

### Build Individual Modes

The `scripts/` directory contains per-mode PowerShell scripts for Windows (MinGW):

```powershell
.\scripts\build-debug.ps1
.\scripts\build-test.ps1
.\scripts\build-release.ps1
```

Each script creates the corresponding `build_debug/`, `build_test/`, or `build_release/` directory.

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

- **`AC_Current_Loop`** — Rectangular current loop in the xy-plane at a given z-level with sinusoidal drive, producing a magnetic dipole radiation pattern
- **`Straight_Wire_X`** — Sinusoidal current along the x-axis between two endpoints
- **`Point_Source`** — Isotropic hard source injection at a single grid point (all three J components)
- **`Gaussian_Pulse`** — Temporally Gaussian current pulse (Jz) at a single point

## CPML Absorbing Boundaries

The solver uses Convolutional Perfectly Matched Layers with Roden–Gedney coefficients to absorb outgoing waves at all six domain faces. The implementation features polynomial grading profiles for conductivity (σ), coordinate stretching (κ), and frequency shifting (α), with separate coefficients computed at E-field integer and H-field half-integer positions for proper Yee-grid alignment. Auxiliary convolution (ψ) arrays are maintained for all six face pairs across all three field directions.

The correct CPML formulation keeps the standard unmodified curl in the main E/H update; the ψ correction encodes the full `(1/s - 1)` stretching.

## Memory Layout

All field and material arrays (E(x,y,z), H(x,y,z), J(x,y,z), ε(x,y,z), μ(x,y,z), σ(x,y,z), and precomputed Ca(x,y,z), Cb(x,y,z), Db(x,y,z) coefficients) are stored in a single monolithic `AlignedSoA<double>` block. Each sub-array stride is rounded up to a SIMD-width boundary, ensuring every row starts at a naturally aligned address. SIMD width is detected at compile time (AVX-512 → 64 B, AVX2 → 32 B, SSE2 → 16 B, scalar → 8 B), avoiding the wasted cache bandwidth of a fixed 64-byte alignment on narrower hardware. Raw pointer access with `RESTRICT` qualification and `ASSUME_ALIGNED` hints gives the compiler full optimization freedom in the hot loops.

## Validation

The built-in plane wave validation test initializes a sinusoidal plane wave propagating in the x-direction with correctly matched E and H fields (Ey = sin(kx), Hz = sin(kx)/η), then tracks three metrics over a dynamically computed number of time steps:

- **Energy drift** — Total electromagnetic energy conservation (pass: <5%)
- **Phase correlation** — Pearson correlation between simulated and analytical Ey at a probe point (pass: >97%)
- **Dispersion error** — Average phase error relative to expected propagation (pass: <10%)

Enable or disable with `run_validation` in `config.cfg`.

## Visualization

The Python renderer (`render.py`) generates an interactive 3D animation using Rerun:

```bash
python render.py
```

Features include volume rendering via Points3D with magnitude-scaled radii, vector cones via Arrows3D showing field direction, separate colormaps (viridis for E-field, inferno for H-field), time-series energy plots (E, H, and leapfrog-synchronized total), and configurable thresholds and downsampling for performance.

## Project Structure

```
├── config.cfg                              # Runtime simulation parameters
├── CMakeLists.txt                          # Three-tier build system
├── render.py                               # Rerun 3D visualization
├── scripts/
│   ├── build-debug.ps1                     # Debug build (Windows/MinGW)
│   ├── build-test.ps1                      # Test build (Windows/MinGW)
│   └── build-release.ps1                   # Release build (Windows/MinGW)
├── src/
│   ├── main.cpp                            # Entry point: config → validation → simulation
│   ├── Utilities/
│   │   ├── aligned_soa.hpp                 # SIMD-aligned monolithic SoA allocator
│   │   └── Macros.hpp                      # RESTRICT, ASSUME_ALIGNED portability macros
│   └── Classes/
│       ├── Config/
│       │   └── config.hpp                  # Simulation_Config with file parser & derived quantities
│       ├── Grid/
│       │   ├── grid.hpp                    # Grid class: field storage, material arrays, indexing
│       │   └── grid.cpp                    # Yee update kernels (H, E), coefficient baking, energy
│       ├── Source/
│       │   ├── source.hpp                  # Polymorphic Source base + Loop/Wire/Point/Gaussian
│       │   └── source.cpp                  # Source apply methods
│       ├── PML/
│       │   ├── PML.hpp                     # CPML class: coefficients, ψ arrays, indexing
│       │   └── PML.cpp                     # Roden–Gedney coefficient computation & ψ updates
│       ├── Simulation/
│       │   ├── Simulation.hpp              # Simulation driver declaration
│       │   └── Simulation.cpp              # Simulation loop, timing, progress output
│       ├── Write_Output/
│       │   ├── output.hpp                  # Double-buffered binary output with writer thread
│       │   └── output.cpp                  # Yee-averaged field export, condition-variable sync
│       └── Validation/
│           ├── Validation.hpp              # Plane wave test declaration
│           └── Validation.cpp              # Validation metrics & reporting
├── tests/
│   ├── test_framework.hpp                  # Zero-dependency header-only test framework
│   ├── test_helpers.hpp                    # Shared test utilities
│   ├── test_main.cpp                       # Test runner entry point
│   ├── test_aligned_soa.cpp                # AlignedSoA allocation & alignment tests (7)
│   ├── test_grid.cpp                       # Grid, Source tests (9)
│   ├── test_pml.cpp                        # CPML coefficient & boundary tests (6)
│   ├── test_integration.cpp                # Full-simulation integration tests (16)
│   └── test_output_validation.cpp          # Output I/O, Validation, CPML reflection tests (5)
├── output/                                 # Generated simulation data (gitignored)
│   ├── E.bin                               # Electric field snapshots
│   └── H.bin                               # Magnetic field snapshots
└── README.md
```

## Build Modes

| Mode | Optimization | Warnings | Tests | Sanitizers | Use Case |
|------|-------------|----------|-------|------------|----------|
| Debug | `-O0` (`-O1` for tests) | Full | Yes | ASan + UBSan (Linux/macOS) | Bug hunting |
| Test | `-O3 -march=native` (lib), `-O2 -march=native` (tests, no `-ffast-math`) | Full | Yes | No | Fast correctness verification |
| Release | `-O3 -march=native -ffast-math -flto` | No | No | No | Benchmarking / production |

## Performance Notes

- OpenMP `collapse(2)` parallelizes the two outer spatial loops; `#pragma omp simd` vectorizes the inner x-loop
- `RESTRICT`-qualified pointers on all hot-path field locals eliminate aliasing penalties
- Loop-invariant quantities (`inv_dx`, `dt_local`, stride constants) are hoisted before the parallel region
- `ASSUME_ALIGNED` macro carries zero runtime overhead — purely a compiler hint applied to coefficient arrays (not ψ arrays, which are accessed at non-aligned offsets)
- Precomputed Ca, Cb, Db coefficient arrays eliminate per-cell material divisions from the hot loop
- Binary output uses a persistent writer thread with double-buffer swap so I/O overlaps computation
- Compilation with `-O3 -march=native` is recommended; `-ffast-math` and LTO are enabled automatically in Release mode