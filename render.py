import numpy as np
import rerun as rr
import os

# Colour Maps:
def inferno(t):
    # Inferno-like colormap, t in [0,1] -> (R,G,B) in [0,255].
    t = np.clip(t, 0.0, 1.0)
    r = np.clip(255 * (1.0122 * t**3 - 0.3378 * t**2 + 0.7446 * t - 0.0013), 0, 255)
    g = np.clip(255 * (-0.8916 * t**3 + 1.5906 * t**2 - 0.2520 * t + 0.0026), 0, 255)
    b = np.clip(255 * (2.3137 * t**3 - 4.8890 * t**2 + 2.6561 * t + 0.0892), 0, 255)
    return np.stack([r, g, b], axis=-1).astype(np.uint8)

def viridis(t):
    t = np.clip(t, 0.0, 1.0)
    r = np.clip(255 * ( t * t * 1.0 ), 0, 255)
    g = np.clip(255 * ( 0.15 + t * 0.75 ), 0, 255)
    b = np.clip(255 * ( 0.55 - t * 0.35 ), 0, 255)
    return np.stack([r, g, b], axis=-1).astype(np.uint8)

# Binary Loader:
def load_binary_volume(f, nx, ny, nz):
    # Read one frame from an already-open file handle.
    count = int(nx * ny * nz * 3)
    data = np.fromfile(f, dtype=np.float64, count=count)
    if len(data) < count:
        return None, None
    grid = data.reshape((int(nz), int(ny), int(nx), 3))
    grid = np.transpose(grid, (2, 1, 0, 3))  # -> (x, y, z, 3)
    mag = np.linalg.norm(grid, axis=-1)
    return mag, grid

# Config:
vol_step    = 2      # Downsample for volume point cloud
arrow_step  = 1      # Downsample for vector arrows
arrow_scale = 1.5    # Arrow length scaling

# Per-field thresholds and radii:
e_mag_thresh   = 0.15
e_arrow_thresh = 0.25
e_radius_min   = 0.10
e_radius_max   = 0.35

h_mag_thresh   = 0.10
h_arrow_thresh = 0.20
h_radius_min   = 0.10
h_radius_max   = 0.25

# File paths:
e_path = "output/E.bin"
h_path = "output/H.bin"

assert os.path.exists(e_path), "E.bin not found in output/"
assert os.path.exists(h_path), "H.bin not found in output/"

# Read header:
with open(e_path, 'rb') as f:
    header = np.fromfile(f, dtype=np.uint64, count=3)
    assert len(header) == 3, "Invalid header in E.bin"

nx, ny, nz = int(header[0]), int(header[1]), int(header[2])
print(f"Grid: {nx} x {ny} x {nz}")

header_bytes = 3 * 8
frame_bytes = nx * ny * nz * 3 * 8
file_size = os.path.getsize(e_path)
num_frames = (file_size - header_bytes) // frame_bytes
print(f"Found {num_frames} frames")

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
h_mag_max = 0.0
e_vec_max = 0.0
h_vec_max = 0.0

with open(e_path, 'rb') as fe, open(h_path, 'rb') as fb:
    fe.seek(header_bytes)
    fb.seek(header_bytes)
    for i in range(num_frames):
        e_mag, e_vec = load_binary_volume(fe, nx, ny, nz)
        h_mag, h_vec = load_binary_volume(fb, nx, ny, nz)

        e_mag_max = max(e_mag_max, e_mag.max())
        h_mag_max = max(h_mag_max, h_mag.max())
        e_vec_max = max(e_vec_max, np.linalg.norm(e_vec, axis=-1).max())
        h_vec_max = max(h_vec_max, np.linalg.norm(h_vec, axis=-1).max())

if e_mag_max < 1e-10: e_mag_max = 1.0
if h_mag_max < 1e-10: h_mag_max = 1.0
if e_vec_max < 1e-10: e_vec_max = 1.0
if h_vec_max < 1e-10: h_vec_max = 1.0

print(f"E mag max: {e_mag_max:.4f}, H mag max: {h_mag_max:.4f}")
print(f"E vec max: {e_vec_max:.4f}, H vec max: {h_vec_max:.4f}")

# Initialize:
rr.init("FDTD_EM_Solver", spawn=True)

rr.log("world", rr.ViewCoordinates.RIGHT_HAND_Z_UP, static=True)

rr.log("world/bounds", rr.Boxes3D(
    centers=[[nx / 2, ny / 2, nz / 2]],
    half_sizes=[[nx / 2, ny / 2, nz / 2]],
    colors=[[80, 80, 80, 40]],
), static=True)

# Initialize caching variable for time-synchronization
prev_h_total = None

# Log Frames:
print("Logging frames to Rerun...")

with open(e_path, 'rb') as fe, open(h_path, 'rb') as fb:
    fe.seek(header_bytes)
    fb.seek(header_bytes)

    for frame_idx in range(num_frames):
        rr.set_time("timestep", sequence=frame_idx)

        e_mag, e_vec = load_binary_volume(fe, nx, ny, nz)
        h_mag, h_vec = load_binary_volume(fb, nx, ny, nz)

        # E-Field Volume:
        e_mag_sub = e_mag[::vol_step, ::vol_step, ::vol_step].ravel()
        e_norm = e_mag_sub / e_mag_max
        e_mask = e_norm > e_mag_thresh

        if e_mask.any():
            e_colors = viridis(e_norm[e_mask])
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
            e_acolors = viridis(e_anorm)

            rr.log("world/E_field/vectors", rr.Arrows3D(
                origins=arrow_origins_full[e_amask],
                vectors=e_vectors,
                colors=e_acolors,
            ))
        else:
            rr.log("world/E_field/vectors", rr.Clear(recursive=False))

        # H-Field Volume:
        h_mag_sub = h_mag[::vol_step, ::vol_step, ::vol_step].ravel()
        h_norm = h_mag_sub / h_mag_max
        h_mask = h_norm > h_mag_thresh

        if h_mask.any():
            h_colors = inferno(h_norm[h_mask])
            h_radii = (h_radius_min + (h_radius_max - h_radius_min) * h_norm[h_mask]).astype(np.float32)

            rr.log("world/H_field/volume", rr.Points3D(
                positions=vol_positions[h_mask],
                colors=h_colors,
                radii=h_radii,
            ))
        else:
            rr.log("world/H_field/volume", rr.Clear(recursive=False))

        # H-field Vectors:
        h_vx = h_vec[::arrow_step, ::arrow_step, ::arrow_step, 0].ravel()
        h_vy = h_vec[::arrow_step, ::arrow_step, ::arrow_step, 1].ravel()
        h_vz = h_vec[::arrow_step, ::arrow_step, ::arrow_step, 2].ravel()
        h_vmag = np.sqrt(h_vx**2 + h_vy**2 + h_vz**2)
        h_amask = h_vmag > h_vec_max * h_arrow_thresh

        if h_amask.any():
            scale = arrow_scale / h_vec_max
            h_vectors = np.column_stack([
                h_vx[h_amask] * scale,
                h_vy[h_amask] * scale,
                h_vz[h_amask] * scale,
            ]).astype(np.float32)

            h_anorm = h_vmag[h_amask] / h_vec_max
            h_acolors = inferno(h_anorm)

            rr.log("world/H_field/vectors", rr.Arrows3D(
                origins=arrow_origins_full[h_amask],
                vectors=h_vectors,
                colors=h_acolors,
            ))
        else:
            rr.log("world/H_field/vectors", rr.Clear(recursive=False))

        # Energy Scalars:
        e_total = float(np.sum(e_mag**2))
        h_total = float(np.sum(h_mag**2))

        rr.log("plots/E_energy", rr.Scalars(e_total))
        rr.log("plots/H_energy", rr.Scalars(h_total))

        # Time-synchronized Total Energy (Interpolating H-field to match E-field time t)
        if prev_h_total is not None:
            sync_h_total = (prev_h_total + h_total) / 2.0
            rr.log("plots/total_energy", rr.Scalars(e_total + sync_h_total))
        else:
            rr.log("plots/total_energy", rr.Scalars(e_total + h_total))

        prev_h_total = h_total

        if (frame_idx + 1) % 10 == 0 or frame_idx == 0:
            print(f"  Frame {frame_idx + 1}/{num_frames}")

print("Done! Use the timeline slider in Rerun to scrub through timesteps.")