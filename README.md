# Kestrel

A C++23 framework in development for GPU-accelerated three-dimensional finite-difference solvers. Kestrel uses XPU for dual CPU/CUDA backends; HIP support is planned.

The existing Maxwell FDTD solver now lives in `examples/maxwell` and serves as the first reference application while the generic Grid and Field APIs are developed.

![C++](https://img.shields.io/badge/C++-23-blue?logo=c%2B%2B)
![OpenMP](https://img.shields.io/badge/OpenMP-Parallel-green)
![Python](https://img.shields.io/badge/Python-3.x-yellow?logo=python)

## Demo

### Combined E & H Fields
![EM Fields](examples/maxwell/assets/em-fields.gif)

### H-Field
![H Field](examples/maxwell/assets/h-field.gif)

### E-Field
![E Field](examples/maxwell/assets/e-field.gif)

*AC current loop radiating electromagnetic waves, visualized with Rerun. E-field in viridis, H-field in inferno.*

## Maxwell Example

Implements the Yee algorithm on a staggered grid using leapfrog time-stepping in H-field formulation with per-cell material coefficients.

### Maxwell's Equations

$$\frac{\partial \mathbf{H}}{\partial t} = -\frac{1}{\mu}\nabla \times \mathbf{E}$$

$$\frac{\partial \mathbf{E}}{\partial t} = \frac{1}{\varepsilon}\nabla \times \mathbf{H} - \frac{\sigma}{\varepsilon}\mathbf{E} - \frac{\mathbf{J}}{\varepsilon}$$

Discretized with precomputed loss-incorporating coefficients:

$$E_i^{n+1} = C_a \cdot E_i^n + C_b \cdot \left(\nabla \times \mathbf{H} - \mathbf{J}\right)_i \qquad H_i^{n+1} = H_i^n - D_b \cdot (\nabla \times \mathbf{E})_i$$

where $C_a = \frac{1 - \sigma\Delta t / 2\varepsilon}{1 + \sigma\Delta t / 2\varepsilon}$, $C_b = \frac{1}{\varepsilon/\Delta t + \sigma/2}$, and $D_b = \Delta t / \mu$.

The time step satisfies the CFL condition: $\Delta t \leq \frac{\alpha}{c \sqrt{1/\Delta x^2 + 1/\Delta y^2 + 1/\Delta z^2}}$

### Features

- **Yee-grid leapfrog** with per-cell ε, μ, σ and precomputed Ca/Cb/Db coefficients via public `bake_coefficients()`
- **Roden-Gedney PML absorbing boundaries** (σ/α grading with κ = 1; no coordinate stretching), per-axis σ_max for anisotropic grids, and material-correct ψ application using per-cell Cb/Db
- **Monolithic `AlignedSoA<T>` memory** - 27 field/material/coefficient arrays in one SIMD-aligned contiguous block with compile-time width detection (AVX-512/AVX2/SSE2)
- **OpenMP parallelization** with `collapse(2)` outer loops, `#pragma omp simd` inner loops, `RESTRICT`/`ASSUME_ALIGNED` hints
- **Four source types**: AC current loop, AC concentric rings, straight wire, point source, Gaussian pulse
- **Double-buffered binary I/O** with persistent writer thread overlapping computation
- **Plane wave validation** tracking energy drift, phase correlation, and numerical dispersion
- **Rerun 3D visualization** with volume rendering, vector cones, and time-series energy plots
- **Tests across 7 suites** including analytical curl checks, multi-axis PML reflection (dB), non-unity material regression, lossy baking, independent H/E half-step verification, and anisotropic coefficient validation

## Dependencies

| Dependency | Version | Required | Purpose |
|-----------|---------|----------|---------|
| Linux or WSL2 | — | Yes | Supported platform |
| GCC | ≥ 15 | Yes | C++23 compiler |
| CMake | ≥ 3.25 | Yes | Build system |
| xpu | `main` | Yes | CPU/GPU portability layer |
| OpenMP | — | Yes | Parallelization |
| Python | ≥ 3.x | No | Visualization only |
| NumPy | — | No | Visualization only |
| Rerun SDK | — | No | Visualization only |

## Build and Run

```bash
cmake -S . -B build-release \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_COMPILER=g++-15
cmake --build build-release -j

cd build-release/examples/maxwell
./kestrel-maxwell
python ../../../examples/maxwell/render.py  # optional visualization
```

Set `KESTREL_BUILD_EXAMPLES=OFF` to configure Kestrel without its examples.

### Build Scripts

The Linux/WSL2 scripts create their build directories on first run and perform
incremental builds on subsequent runs.

Linux / WSL2:

```bash
./scripts/build-debug.sh
./scripts/build-test.sh
./scripts/build-release.sh
```

### Build Modes

| Mode | Flags | Tests | Use Case |
|------|-------|-------|----------|
| Debug | `-O0 -g3`, ASan+UBSan | Yes | Bug hunting |
| Test | `-O3 -march=native` (lib), `-O2` (tests, no `-ffast-math`) | Yes | Correctness |
| Release | `-O3 -march=native -ffast-math -flto` | No | Benchmarking |

### Running Tests

```bash
cmake -S . -B build-test \
    -DCMAKE_BUILD_TYPE=Test \
    -DCMAKE_CXX_COMPILER=g++-15
cmake --build build-test --target kestrel-maxwell-tests -j
ctest --test-dir build-test --output-on-failure

# or run individual Maxwell suites:
./build-test/examples/maxwell/kestrel-maxwell-tests PML
./build-test/examples/maxwell/kestrel-maxwell-tests Integration
```

## Configuration

Parameters in `examples/maxwell/config.cfg` use `key = value` syntax. Missing keys use compiled defaults; unknown keys warn.

| Parameter | Description | Default |
|-----------|-------------|---------|
| `Nx, Ny, Nz` | Grid cell count per axis | 100 |
| `dx, dy, dz` | Spatial step size | 1.0 |
| `eps, mu` | Permittivity and permeability | 1.0 |
| `cfl_factor` | CFL stability factor (0 < α ≤ 1) | 1.0 |
| `total_steps` | Number of time steps | 1000 |
| `run_validation` | Run plane wave test on startup | true |
| `use_pml` | Enable PML boundaries | true |
| `pml_order` | Polynomial grading order | 3 |
| `pml_alpha_max` | Maximum α frequency shift | 0.3 |

Derived: `dt` (CFL-limited), `pml_thickness` (from grid volume, clamped to [4, 25]), `pml_sigma_max_x/y/z` (per-axis Berenger formula).

## Project Structure

```text
├── CMakeLists.txt                         # Shared toolchain, XPU, and backend setup
├── examples/
│   └── maxwell/
│       ├── CMakeLists.txt                 # Maxwell targets
│       ├── config.cfg                     # Runtime parameters
│       ├── render.py                      # Rerun visualization
│       ├── assets/                        # Demo media
│       ├── benchmark/                     # Maxwell benchmark
│       ├── src/                           # Existing Maxwell implementation
│       │   ├── main.cu
│       │   ├── utilities/
│       │   └── classes/
│       └── tests/                         # Existing Maxwell test suite
└── scripts/                               # Linux and WSL2 workflows
```

The future framework API will live under `include/kestrel`; framework tests will be introduced alongside the Grid and Field abstractions.
