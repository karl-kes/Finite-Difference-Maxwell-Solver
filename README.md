# 3D Finite-Difference Time-Domain Maxwell Solver

A high-performance Finite-Difference Time-Domain (FDTD) solver for Maxwell's equations in 3D, implemented in C++ with OpenMP parallelization and interactive Plotly visualization.

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

The time step is automatically calculated to satisfy the Courant-Friedrichs-Lewy condition:

$$\Delta t \leq \frac{\alpha}{c \sqrt{\frac{1}{\Delta x^2} + \frac{1}{\Delta y^2} + \frac{1}{\Delta z^2}}}$$

where α (0 < α <= 1.0) is a configurable CFL factor (default 0.1) for stability margin.

## Features

- **Yee staggering** for proper spatial discretization of E and B fields
- **Current density (J)** support for realistic source modeling
- **Multiple source types**: hard sources, dipole antennas, straight wire currents
- **Energy conservation tracking** for simulation validation
- **OpenMP parallelization** with `collapse(2)` and `simd` pragmas
- **Binary I/O with buffering** for efficient data export
- **Interactive 3D visualization** with synchronized E and B field animations

## Quick Start

### Prerequisites

- C++17 compiler (GCC recommended)
- OpenMP support
- Python 3.x with NumPy and Plotly

### Build & Run

```bash
# Compile with OpenMP (recommended)
g++ -std=c++17 main.cpp
               Classes/Grid/grid_constructor.cpp
               Classes/Grid/grid_getters.cpp
               Classes/Grid/grid_simulation.cpp
               Classes/Grid/grid_helpers.cpp
               -o main.exe -fopenmp

# Run simulation
./main.exe

# Visualize results
python render.py
```

### Without OpenMP

```bash
g++ -std=c++17 main.cpp
               Classes/Grid/grid_constructor.cpp
               Classes/Grid/grid_getters.cpp
               Classes/Grid/grid_simulation.cpp
               Classes/Grid/grid_helpers.cpp
               -o main.exe
```

## Configuration

All simulation parameters are centralized in `config.hpp`:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `Nx, Ny, Nz` | Grid dimensions | 15 |
| `dx, dy, dz` | Spatial step sizes | 5.0 |
| `eps, mu` | Permittivity and permeability | 1.0 |
| `cfl_factor` | CFL stability factor | 0.1 |
| `total_time` | Number of time steps | 1000 |
| `amp_one, amp_two` | Source amplitudes | 100.0 |
| `freq_one, freq_two` | Source frequencies | 1.0 |

## Source Types

The simulator supports several excitation methods:

- **`straight_wire_x()`** - Sinusoidal current along x-axis
- **`hard_source_inject()`** - Direct field injection at a point
- **`dipole_antenna_inject()`** - Dual-point oscillating source

## Visualization

The Python renderer generates an interactive side-by-side animation showing E and B field evolution with volume rendering and vector cones:

```bash
python render.py
```

Features:
- Synchronized E and B field display
- Volume rendering for field intensity
- Vector cones showing field direction
- Play/pause animation controls

## Project Structure

```
├── main.cpp                           # Simulation driver
├── config.hpp                         # Centralized configuration
├── Classes/Grid/
│   ├── grid.hpp                       # Grid class declaration
│   ├── grid_constructor.cpp           # Initialization
│   ├── grid_getters.cpp               # Accessor methods
│   ├── grid_simulation.cpp            # Field updates & sources
│   └── grid_helpers.cpp               # Utilities & energy calc
├── render.py                          # Plotly visualization
├── output/                            # Generated data
│   ├── E/                             # Electric field snapshots
│   └── B/                             # Magnetic field snapshots
└── README.md
```

## Output

The simulation prints:
- Progress percentage during execution
- Total simulation duration (ms)
- Physical time simulated (s)
- Maximum energy drift (%) for validation

## Performance Notes

The modular file structure enables better compiler optimization through improved instruction cache locality. OpenMP parallelization with `collapse(2)` optimizes the nested spatial loops, while `simd` pragmas enable vectorization of the innermost loop.