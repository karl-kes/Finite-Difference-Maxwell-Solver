# Finite-Difference Maxwell Solver

A high-performance 3D Finite-Difference Time-Domain (FDTD) solver for Maxwell's equations, implemented in C++ with OpenMP parallelization, CPML absorbing boundaries, and interactive Plotly visualization.

![C++](https://img.shields.io/badge/C++-17-blue?logo=c%2B%2B)
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

where α (0 < α ≤ 1.0) is a configurable CFL factor (default 0.125) for stability margin.

## Features

- **Yee staggering** for proper spatial discretization of E and B fields
- **CPML absorbing boundaries** (Convolutional Perfectly Matched Layer) with polynomial grading for reflection-free domain truncation
- **Current density (J)** support for realistic source modeling
- **Multiple source types**: point sources, straight wire currents, and Gaussian pulses
- **Plane wave validation test** with energy drift, phase correlation, and dispersion metrics
- **Energy conservation tracking** and source power diagnostics
- **OpenMP parallelization** with `collapse(2)` and `simd` pragmas
- **Binary I/O with buffered writes** for efficient data export
- **Interactive 3D visualization** with synchronized E and B field volume rendering and vector cones

## Quick Start

### Prerequisites

- C++17 compiler (GCC recommended)
- OpenMP support
- Python 3.x with NumPy and Plotly

### Build & Run

```bash
# Compile with OpenMP
g++ -std=c++17 -O3 -march=native -fopenmp src/main.cpp src/Classes/Grid/*.cpp src/Classes/Source/*.cpp src/Classes/Write_Output/*.cpp src/Classes/PML/*.cpp src/Classes/Simulation/*.cpp src/Classes/Validation/*.cpp -o main.exe

# Run simulation
./main.exe

# Visualize results
python src/render.py
```

### CMake Build

```bash
rm -r build
mkdir build
cd build
cmake .. -G "MinGW Makefiles"
mingw32-make -j
./main.exe
python ../src/render.py
```

## Configuration

All simulation parameters are centralized in `config.hpp`:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `Nx, Ny, Nz` | Grid dimensions | 100 |
| `dx, dy, dz` | Spatial step sizes | 5.0 |
| `eps, mu` | Permittivity and permeability | 1.0 |
| `cfl_factor` | CFL stability factor | 0.125 |
| `total_time` | Number of time steps | 1000 |
| `run_validation` | Run plane wave test on startup | true |
| `use_pml` | Enable CPML boundaries | true |
| `pml_thickness` | PML layer count | 8 |
| `pml_order` | Polynomial grading order | 3 |
| `pml_kappa_max` | Maximum kappa stretching | 15.0 |

## Source Types

The simulator supports several excitation methods through a polymorphic `Source` interface:

- **`Point_Source`** — Isotropic hard source injection at a single grid point
- **`Straight_Wire_X`** — Sinusoidal current along the x-axis between two endpoints
- **`Gaussian_Pulse`** — Temporally Gaussian current pulse at a single point

## CPML Absorbing Boundaries

The solver uses Convolutional Perfectly Matched Layers to absorb outgoing waves at domain boundaries. The CPML implementation features polynomial grading profiles for conductivity (σ), coordinate stretching (κ), and frequency shifting (α), with separate coefficients computed at E-field integer and B-field half-integer positions for proper Yee-grid alignment. Auxiliary convolution (ψ) arrays are maintained for all six face pairs across all three field directions.

## Validation

The built-in plane wave validation test initializes a sinusoidal plane wave propagating in the x-direction, then tracks three metrics over 100 time steps:

- **Energy drift** — Total electromagnetic energy conservation (pass threshold: <5%)
- **Phase correlation** — Pearson correlation between simulated and analytical Ey at a probe point (pass threshold: >0.99)
- **Dispersion error** — Average phase error relative to expected propagation (pass threshold: <10%)

Enable or disable with `run_validation` in `config.hpp`.

## Visualization

The Python renderer generates an interactive side-by-side animation showing E and B field evolution:

```bash
python src/render.py
```

Features include volume rendering for field intensity (Inferno colorscale), vector cones showing field direction and magnitude, and play/pause animation controls.

## Project Structure

```
├── src/
│   ├── main.cpp                        # Simulation driver
│   ├── render.py                       # Plotly visualization
│   └── Classes/
│       ├── Config/
│       │   └── config.hpp              # Centralized configuration & enums
│       ├── Grid/
│       │   ├── grid.hpp                # Grid class declaration
│       │   └── grid.cpp                # Field updates, curl, energy diagnostics
│       ├── Source/
│       │   ├── source.hpp              # Source base class & implementations
│       │   └── source.cpp              # Source apply methods
│       ├── Write_Output/
│       │   ├── output.hpp              # Output class declaration
│       │   └── output.cpp              # Binary field writer with Yee averaging
│       ├── PML/
│       │   ├── PML.hpp                 # CPML class declaration
│       │   └── PML.cpp                 # CPML coefficient computation & psi updates
│       ├── Simulation/
│       │   ├── Simulation.hpp          # Simulation class declaration
│       │   └── Simulation.cpp          # Simulation loop & timing
│       └── Validation/
│           ├── Validation.hpp          # Plane wave test declaration
│           └── Validation.cpp          # Validation metrics & reporting
├── output/                             # Generated simulation data
│   ├── E/                              # Electric field snapshots (.bin)
│   └── B/                              # Magnetic field snapshots (.bin)
└── README.md
```

## Output

The simulation prints progress during execution, followed by total duration (ms) and physical time simulated (s). When validation is enabled, a detailed report shows energy drift, phase correlation, and dispersion error with pass/fail status.

## Performance Notes

OpenMP parallelization with `collapse(2)` optimizes the nested spatial loops, while `simd` pragmas enable vectorization of the innermost loop. The `__restrict__` qualifier on field pointers helps the compiler avoid aliasing penalties. Compilation with `-O3 -march=native` is recommended for best performance.