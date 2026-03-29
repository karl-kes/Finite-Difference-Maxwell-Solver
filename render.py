import numpy as np
import rerun as rr
import glob
import os

# Colour Maps:
def inferno(t):
    """Inferno-like colormap, t in [0,1] -> (R,G,B) in [0,255]."""
    t = np.clip(t, 0.0, 1.0)
    r = np.clip(255 * (1.0122 * t**3 - 0.3378 * t**2 + 0.7446 * t - 0.0013), 0, 255)
    g = np.clip(255 * (-0.8916 * t**3 + 1.5906 * t**2 - 0.2520 * t + 0.0026), 0, 255)
    b = np.clip(255 * (2.3137 * t**3 - 4.8890 * t**2 + 2.6561 * t + 0.0892), 0, 255)
    return np.stack([r, g, b], axis=-1).astype(np.uint8)

def cyan(t):
    """Cyan/teal colormap, t in [0,1] -> (R,G,B) in [0,255]."""
    t = np.clip(t, 0.0, 1.0)
    r = np.clip(255 * (0.1 * t), 0, 255)
    g = np.clip(255 * (0.4 + 0.6 * t), 0, 255)
    b = np.clip(255 * (0.5 + 0.5 * t), 0, 255)
    return np.stack([r, g, b], axis=-1).astype(np.uint8)

# Binary Loader:
def load_binary_volume(file_path):
    with open(file_path, 'rb') as f:
        header = np.fromfile(f, dtype=np.uint64, count=3)
        if len(header) < 3:
            return None, None, (0, 0, 0)
        nx, ny, nz = header
        data = np.fromfile(f, dtype=np.float64, count=int(nx * ny * nz * 4))

    grid = data.reshape((int(nz), int(ny), int(nx), 4))
    grid = np.transpose(grid, (2, 1, 0, 3))  # -> (x, y, z, 4)

    mag = grid[:, :, :, 3]
    vectors = grid[:, :, :, :3]
    return mag, vectors, (int(nx), int(ny), int(nz))

# Config:
vol_step    = 2      # Downsample for volume point cloud
arrow_step  = 1      # Downsample for vector arrows
arrow_scale = 1.5    # Arrow length scaling

# Per-field thresholds and radii:
e_mag_thresh   = 0.15
e_arrow_thresh = 0.08
e_radius_min   = 0.10
e_radius_max   = 0.35

b_mag_thresh   = 0.10
b_arrow_thresh = 0.20
b_radius_min   = 0.10
b_radius_max   = 0.25

# Load Files:
e_files = sorted(
    glob.glob("output/E/E*.bin"),
    key=lambda x: int(''.join(filter(str.isdigit, os.path.basename(x))))
)
b_files = sorted(
    glob.glob("output/B/B*.bin"),
    key=lambda x: int(''.join(filter(str.isdigit, os.path.basename(x))))
)

assert len(e_files) > 0, "No E-field files found in output/E/"
assert len(b_files) > 0, "No B-field files found in output/B/"

num_frames = min(len(e_files), len(b_files))
print(f"Found {num_frames} frames")

_, _, shape = load_binary_volume(e_files[0])
nx, ny, nz = shape
print(f"Grid: {nx} x {ny} x {nz}")

# Compute Grid Coords:
vx, vy, vz = np.mgrid[0:nx:vol_step, 0:ny:vol_step, 0:nz:vol_step]
vol_positions = np.column_stack([
    vx.ravel().astype(np.float32),
    vy.ravel().astype(np.float32),
    vz.ravel().astype(np.float32),
])

ax, ay, az = np.mgrid[0:nx:arrow_step, 0:ny:arrow_step, 0:nz:arrow_step]
arrow_origins_full = np.column_stack([
    ax.ravel().astype(np.float32),
    ay.ravel().astype(np.float32),
    az.ravel().astype(np.float32),
])

# Global Extrema:
print("Scanning for global max values...")
e_mag_max = 0.0
b_mag_max = 0.0
e_vec_max = 0.0
b_vec_max = 0.0

for i in range(num_frames):
    e_mag, e_vec, _ = load_binary_volume(e_files[i])
    b_mag, b_vec, _ = load_binary_volume(b_files[i])

    e_mag_max = max(e_mag_max, np.sqrt(np.maximum(e_mag, 0)).max())
    b_mag_max = max(b_mag_max, np.sqrt(np.maximum(b_mag, 0)).max())
    e_vec_max = max(e_vec_max, np.linalg.norm(e_vec, axis=-1).max())
    b_vec_max = max(b_vec_max, np.linalg.norm(b_vec, axis=-1).max())

if e_mag_max < 1e-10: e_mag_max = 1.0
if b_mag_max < 1e-10: b_mag_max = 1.0
if e_vec_max < 1e-10: e_vec_max = 1.0
if b_vec_max < 1e-10: b_vec_max = 1.0

print(f"E mag max: {e_mag_max:.4f}, B mag max: {b_mag_max:.4f}")
print(f"E vec max: {e_vec_max:.4f}, B vec max: {b_vec_max:.4f}")

# Initialize:
rr.init("FDTD_EM_Solver", spawn=True)

rr.log("world", rr.ViewCoordinates.RIGHT_HAND_Z_UP, static=True)

rr.log("world/bounds", rr.Boxes3D(
    centers=[[nx / 2, ny / 2, nz / 2]],
    half_sizes=[[nx / 2, ny / 2, nz / 2]],
    colors=[[80, 80, 80, 40]],
), static=True)

# Log Frames:
print("Logging frames to Rerun...")

for frame_idx in range(num_frames):
    rr.set_time("timestep", sequence=frame_idx)

    e_mag_raw, e_vec, _ = load_binary_volume(e_files[frame_idx])
    b_mag_raw, b_vec, _ = load_binary_volume(b_files[frame_idx])

    e_mag = np.sqrt(np.maximum(e_mag_raw, 0))
    b_mag = np.sqrt(np.maximum(b_mag_raw, 0))

    # E-Field Volume:
    e_mag_sub = e_mag[::vol_step, ::vol_step, ::vol_step].ravel()
    e_norm = e_mag_sub / e_mag_max
    e_mask = e_norm > e_mag_thresh

    if e_mask.any():
        e_colors = cyan(e_norm[e_mask])
        e_radii = (e_radius_min + (e_radius_max - e_radius_min) * e_norm[e_mask]).astype(np.float32)

        rr.log("world/E_field/volume", rr.Points3D(
            positions=vol_positions[e_mask],
            colors=e_colors,
            radii=e_radii,
        ))
    else:
        rr.log("world/E_field/volume", rr.Clear(recursive=False))

    # E-Field Vectors:
    e_vx = e_vec[::arrow_step, ::arrow_step, ::arrow_step, 0].ravel()
    e_vy = e_vec[::arrow_step, ::arrow_step, ::arrow_step, 1].ravel()
    e_vz = e_vec[::arrow_step, ::arrow_step, ::arrow_step, 2].ravel()
    e_vmag = np.sqrt(e_vx**2 + e_vy**2 + e_vz**2)
    e_amask = e_vmag > e_vec_max * e_arrow_thresh

    if e_amask.any():
        scale = arrow_scale / e_vec_max
        e_vectors = np.column_stack([
            e_vx[e_amask] * scale,
            e_vy[e_amask] * scale,
            e_vz[e_amask] * scale,
        ]).astype(np.float32)

        e_anorm = e_vmag[e_amask] / e_vec_max
        e_acolors = cyan(e_anorm)

        rr.log("world/E_field/vectors", rr.Arrows3D(
            origins=arrow_origins_full[e_amask],
            vectors=e_vectors,
            colors=e_acolors,
        ))
    else:
        rr.log("world/E_field/vectors", rr.Clear(recursive=False))

    # B-Field Volume:
    b_mag_sub = b_mag[::vol_step, ::vol_step, ::vol_step].ravel()
    b_norm = b_mag_sub / b_mag_max
    b_mask = b_norm > b_mag_thresh

    if b_mask.any():
        b_colors = inferno(b_norm[b_mask])
        b_radii = (b_radius_min + (b_radius_max - b_radius_min) * b_norm[b_mask]).astype(np.float32)

        rr.log("world/B_field/volume", rr.Points3D(
            positions=vol_positions[b_mask],
            colors=b_colors,
            radii=b_radii,
        ))
    else:
        rr.log("world/B_field/volume", rr.Clear(recursive=False))

    # B-Field Vectors:
    b_vx = b_vec[::arrow_step, ::arrow_step, ::arrow_step, 0].ravel()
    b_vy = b_vec[::arrow_step, ::arrow_step, ::arrow_step, 1].ravel()
    b_vz = b_vec[::arrow_step, ::arrow_step, ::arrow_step, 2].ravel()
    b_vmag = np.sqrt(b_vx**2 + b_vy**2 + b_vz**2)
    b_amask = b_vmag > b_vec_max * b_arrow_thresh

    if b_amask.any():
        scale = arrow_scale / b_vec_max
        b_vectors = np.column_stack([
            b_vx[b_amask] * scale,
            b_vy[b_amask] * scale,
            b_vz[b_amask] * scale,
        ]).astype(np.float32)

        b_anorm = b_vmag[b_amask] / b_vec_max
        b_acolors = inferno(b_anorm)

        rr.log("world/B_field/vectors", rr.Arrows3D(
            origins=arrow_origins_full[b_amask],
            vectors=b_vectors,
            colors=b_acolors,
        ))
    else:
        rr.log("world/B_field/vectors", rr.Clear(recursive=False))

    # Energy Scalars:
    e_total = float(np.sum(e_mag**2))
    b_total = float(np.sum(b_mag**2))
    rr.log("plots/E_energy", rr.Scalars(e_total))
    rr.log("plots/B_energy", rr.Scalars(b_total))
    rr.log("plots/total_energy", rr.Scalars(e_total + b_total))

    if (frame_idx + 1) % 10 == 0 or frame_idx == 0:
        print(f"  Frame {frame_idx + 1}/{num_frames}")

print("Done! Use the timeline slider in Rerun to scrub through timesteps.")