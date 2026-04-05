# Finite-Difference Maxwell Solver

A high-performance 3D FDTD solver for Maxwell's equations in C++ with OpenMP parallelization, CPML absorbing boundaries, and Rerun-based 3D visualization.

![C++](https://img.shields.io/badge/C++-23-blue?logo=c%2B%2B)
![OpenMP](https://img.shields.io/badge/OpenMP-Parallel-green)
![Python](https://img.shields.io/badge/Python-3.x-yellow?logo=python)

## Demo

### Combined E & H Fields
![EM Fields](Assets/EM_Fields.gif)

### H-Field
![H Field](Assets/H_Field.gif)

### E-Field
![E Field](Assets/E_Field.gif)

*AC current loop radiating electromagnetic waves, visualized with Rerun. E-field in viridis, H-field in inferno.*

## Overview

Implements the Yee algorithm on a staggered grid using leapfrog time-stepping in H-field formulation with per-cell material coefficients.

### Maxwell's Equations

$$\frac{\partial \mathbf{H}}{\partial t} = -\frac{1}{\mu}\nabla \times \mathbf{E}$$

$$\frac{\partial \mathbf{E}}{\partial t} = \frac{1}{\varepsilon}\nabla \times \mathbf{H} - \frac{\sigma}{\varepsilon}\mathbf{E} - \frac{\mathbf{J}}{\varepsilon}$$

Discretized with precomputed loss-incorporating coefficients:

$$E_i^{n+1} = C_a \cdot E_i^n + C_b \cdot \left(\nabla \times \mathbf{H} - \mathbf{J}\right)_i \qquad H_i^{n+1} = H_i^n - D_b \cdot (\nabla \times \mathbf{E})_i$$

where $C_a = \frac{1 - \sigma\Delta t / 2\varepsilon}{1 + \sigma\Delta t / 2\varepsilon}$, $C_b = \frac{1}{\varepsilon/\Delta t + \sigma/2}$, and $D_b = \Delta t / \mu$.

The time step satisfies the CFL condition: $\Delta t \leq \frac{\alpha}{c \sqrt{1/\Delta x^2 + 1/\Delta y^2 + 1/\Delta z^2}}$

## Features

- **Yee-grid leapfrog** with per-cell ε, μ, σ and precomputed Ca/Cb/Db coefficients via public `bake_coefficients()`
- **CPML absorbing boundaries** with Roden–Gedney coefficients, per-axis σ_max for anisotropic grids, and material-correct ψ application using per-cell Cb/Db
- **Monolithic `AlignedSoA<T>` memory** - 27 field/material/coefficient arrays in one SIMD-aligned contiguous block with compile-time width detection (AVX-512/AVX2/SSE2)
- **OpenMP parallelization** with `collapse(2)` outer loops, `#pragma omp simd` inner loops, `RESTRICT`/`ASSUME_ALIGNED` hints
- **Four source types**: AC current loop, AC concentric rings, straight wire, point source, Gaussian pulse
- **Double-buffered binary I/O** with persistent writer thread overlapping computation
- **Plane wave validation** tracking energy drift, phase correlation, and numerical dispersion
- **Rerun 3D visualization** with volume rendering, vector cones, and time-series energy plots
- **55 tests across 7 suites** including analytical curl checks, multi-axis CPML reflection (dB), non-unity material regression, lossy baking, independent H/E half-step verification, and anisotropic coefficient validation

## Dependencies

| Dependency | Version | Required | Purpose |
|-----------|---------|----------|---------|
| C++ compiler (GCC recommended) | C++23 | Yes | Solver |
| CMake | ≥ 3.20 | Yes | Build system |
| OpenMP | — | Yes | Parallelization |
| Python | ≥ 3.x | No | Visualization only |
| NumPy | — | No | Visualization only |
| Rerun SDK | — | No | Visualization only |

## Build & Run

```bash
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j
./main
python ../render.py  # optional visualization
```

Windows (MinGW):

```bash
cmake .. -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release
mingw32-make -j
```

### Build Scripts

The `scripts/` directory has per-mode scripts for both platforms. Each creates its build directory on first run and does incremental builds on subsequent runs.

Linux / macOS:

```bash
./scripts/build-debug.sh
./scripts/build-test.sh
./scripts/build-release.sh
```

Windows (PowerShell / MinGW):

```powershell
.\scripts\build-debug.ps1
.\scripts\build-test.ps1
.\scripts\build-release.ps1
```

### Build Modes

| Mode | Flags | Tests | Use Case |
|------|-------|-------|----------|
| Debug | `-O0 -g3`, ASan+UBSan (Linux/macOS) | Yes | Bug hunting |
| Test | `-O3 -march=native` (lib), `-O2` (tests, no `-ffast-math`) | Yes | Correctness |
| Release | `-O3 -march=native -ffast-math -flto` | No | Benchmarking |

### Running Tests

```bash
cmake .. -DCMAKE_BUILD_TYPE=Test
make -j
ctest --output-on-failure

# or run individual suites:
./run_tests PML
./run_tests Integration
```

## Configuration

Parameters in `config.cfg` (key = value, `#` comments). Missing keys use compiled defaults; unknown keys warn.

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

Derived: `dt` (CFL-limited), `pml_thickness` (from grid volume, clamped to [4, 25]), `pml_sigma_max_x/y/z` (per-axis Berenger formula).

## Project Structure

```
├── config.cfg                          # Runtime parameters
├── CMakeLists.txt                      # Three-tier build system
├── render.py                           # Rerun 3D visualization
├── scripts/
│   ├── build-debug.{sh,ps1}            # Debug build
│   ├── build-test.{sh,ps1}             # Test build
│   └── build-release.{sh,ps1}          # Release build
├── src/
│   ├── main.cpp                        # Entry point
│   ├── Utilities/
│   │   ├── aligned_soa.hpp             # SIMD-aligned SoA allocator
│   │   └── macros.hpp                  # RESTRICT, ASSUME_ALIGNED macros
│   └── Classes/
│       ├── Config/config.hpp           # Config parser & derived quantities
│       ├── Grid/grid.{hpp,cpp}         # Field storage, Yee kernels, bake_coefficients
│       ├── Source/source.{hpp,cpp}     # Source types (loop, wire, point, Gaussian)
│       ├── PML/pml.{hpp,cpp}           # CPML coefficients & ψ updates
│       ├── Simulation/simulation.{hpp,cpp}
│       ├── Write_Output/output.{hpp,cpp}
│       └── Validation/validation.{hpp,cpp}
├── tests/
│   ├── test_framework.hpp              # Header-only test framework
│   ├── test_aligned_soa.cpp            # Allocation & alignment (7)
│   ├── test_grid.cpp                   # Grid, sources, config, lossy, H/E split (19)
│   ├── test_pml.cpp                    # Coefficients, anisotropic, clamping (9)
│   ├── test_integration.cpp            # Physics: conservation, causality, dipole, superposition (14)
│   └── test_output_validation.cpp      # I/O, validation, multi-axis CPML reflection (6)
└── output/                             # Generated data (gitignored)
```