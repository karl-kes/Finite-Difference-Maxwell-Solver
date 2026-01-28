import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import glob
import os

def load_binary_volume_full(file_path):
    with open(file_path, 'rb') as f:
        header = np.fromfile(f, dtype=np.uint64, count=3)
        if len(header) < 3: return None, None, (0,0,0)
        nx, ny, nz = header
        expected_size = int(nx * ny * nz * 4)
        data = np.fromfile(f, dtype=np.float64, count=expected_size)
    
    grid = data.reshape((nz, ny, nx, 4))
    mag = grid[:, :, :, 3]
    vectors = grid[:, :, :, :3]
    return mag, vectors, (nx, ny, nz)

# 1. Load E and B files
e_files = sorted(glob.glob("output/E/E*.bin"), 
                 key=lambda x: int(''.join(filter(str.isdigit, os.path.basename(x)))))
b_files = sorted(glob.glob("output/B/B*.bin"), 
                 key=lambda x: int(''.join(filter(str.isdigit, os.path.basename(x)))))

e_volume_list, e_vector_list = [], []
b_volume_list, b_vector_list = [], []
shape = (0, 0, 0)

for ef, bf in zip(e_files, b_files):
    e_mag, e_vec, shape = load_binary_volume_full(ef)
    b_mag, b_vec, _ = load_binary_volume_full(bf)
    
    if e_mag is not None:
        e_volume_list.append(np.sqrt(np.nan_to_num(np.maximum(e_mag, 0))))
        e_vector_list.append(e_vec)
    if b_mag is not None:
        b_volume_list.append(np.sqrt(np.nan_to_num(np.maximum(b_mag, 0))))
        b_vector_list.append(b_vec)

nx, ny, nz = shape
z, y, x = np.mgrid[0:nz, 0:ny, 0:nx]

step = 1
x_sub = x[::step, ::step, ::step]
y_sub = y[::step, ::step, ::step]
z_sub = z[::step, ::step, ::step]

# 2. Scaling
e_vmax = max([v.max() for v in e_volume_list]) if e_volume_list else 1.0
b_vmax = max([v.max() for v in b_volume_list]) if b_volume_list else 1.0
if e_vmax < 1e-10: e_vmax = 1.0
if b_vmax < 1e-10: b_vmax = 1.0

e_vec_max = max([np.sqrt((v[:,:,:,0]**2 + v[:,:,:,1]**2 + v[:,:,:,2]**2)).max() for v in e_vector_list]) if e_vector_list else 1.0
b_vec_max = max([np.sqrt((v[:,:,:,0]**2 + v[:,:,:,1]**2 + v[:,:,:,2]**2)).max() for v in b_vector_list]) if b_vector_list else 1.0
if e_vec_max < 1e-10: e_vec_max = 1.0
if b_vec_max < 1e-10: b_vec_max = 1.0

# 3. Build Frames
frames = []
for i in range(len(e_volume_list)):
    e_vol, e_vec = e_volume_list[i], e_vector_list[i]
    b_vol, b_vec = b_volume_list[i], b_vector_list[i]
    
    e_Fx = e_vec[::step, ::step, ::step, 0].flatten()
    e_Fy = e_vec[::step, ::step, ::step, 1].flatten()
    e_Fz = e_vec[::step, ::step, ::step, 2].flatten()
    e_mag_sub = np.sqrt(e_Fx**2 + e_Fy**2 + e_Fz**2)
    e_mask = (e_mag_sub > e_vec_max * 0.1).flatten()
    
    b_Fx = b_vec[::step, ::step, ::step, 0].flatten()
    b_Fy = b_vec[::step, ::step, ::step, 1].flatten()
    b_Fz = b_vec[::step, ::step, ::step, 2].flatten()
    b_mag_sub = np.sqrt(b_Fx**2 + b_Fy**2 + b_Fz**2)
    b_mask = (b_mag_sub > b_vec_max * 0.1).flatten()
    
    frames.append(go.Frame(
        data=[
            go.Volume(
                x=x.flatten(), y=y.flatten(), z=z.flatten(),
                value=e_vol.flatten(),
                isomin=e_vmax * 0.01, isomax=e_vmax,
                opacity=0.08, surface_count=15,
                colorscale='Inferno',
                scene='scene'
            ),
            go.Cone(
                x=x_sub.flatten()[e_mask], y=y_sub.flatten()[e_mask], z=z_sub.flatten()[e_mask],
                u=e_Fx[e_mask], v=e_Fy[e_mask], w=e_Fz[e_mask],
                colorscale='Reds', reversescale=False,
                sizemode='scaled', sizeref=0.5, anchor='center',
                showscale=False,
                scene='scene'
            ),
            go.Volume(
                x=x.flatten(), y=y.flatten(), z=z.flatten(),
                value=b_vol.flatten(),
                isomin=b_vmax * 0.01, isomax=b_vmax,
                opacity=0.08, surface_count=15,
                colorscale='Inferno',
                scene='scene2'
            ),
            go.Cone(
                x=x_sub.flatten()[b_mask], y=y_sub.flatten()[b_mask], z=z_sub.flatten()[b_mask],
                u=b_Fx[b_mask], v=b_Fy[b_mask], w=b_Fz[b_mask],
                colorscale='Reds', reversescale=False,
                sizemode='scaled', sizeref=0.5, anchor='center',
                showscale=False,
                scene='scene2'
            )
        ],
        name=str(i)
    ))

# 4. Initial frame
start_idx = next((i for i, v in enumerate(e_volume_list) if v.max() > 0.001), 0)

e_vol, e_vec = e_volume_list[start_idx], e_vector_list[start_idx]
b_vol, b_vec = b_volume_list[start_idx], b_vector_list[start_idx]

e_Fx = e_vec[::step, ::step, ::step, 0].flatten()
e_Fy = e_vec[::step, ::step, ::step, 1].flatten()
e_Fz = e_vec[::step, ::step, ::step, 2].flatten()
e_mag_sub = np.sqrt(e_Fx**2 + e_Fy**2 + e_Fz**2)
e_mask = e_mag_sub > e_vec_max * 0.1

b_Fx = b_vec[::step, ::step, ::step, 0].flatten()
b_Fy = b_vec[::step, ::step, ::step, 1].flatten()
b_Fz = b_vec[::step, ::step, ::step, 2].flatten()
b_mag_sub = np.sqrt(b_Fx**2 + b_Fy**2 + b_Fz**2)
b_mask = b_mag_sub > b_vec_max * 0.1

# 5. Simulate
fig = go.Figure(
    data=[
        go.Volume(
            x=x.flatten(), y=y.flatten(), z=z.flatten(),
            value=e_vol.flatten(),
            isomin=e_vmax * 0.05, isomax=e_vmax,
            opacity=0.08, surface_count=15,
            colorscale='Inferno',
            colorbar=dict(title="E Field", x=0.45),
            scene='scene'
        ),
        go.Cone(
            x=x_sub.flatten()[e_mask], y=y_sub.flatten()[e_mask], z=z_sub.flatten()[e_mask],
            u=e_Fx[e_mask], v=e_Fy[e_mask], w=e_Fz[e_mask],
            colorscale='Reds', reversescale=False,
            sizemode='scaled', sizeref=0.5, anchor='center',
            showscale=False,
            scene='scene'
        ),
        go.Volume(
            x=x.flatten(), y=y.flatten(), z=z.flatten(),
            value=b_vol.flatten(),
            isomin=b_vmax * 0.05, isomax=b_vmax,
            opacity=0.08, surface_count=15,
            colorscale='Inferno',
            colorbar=dict(title="B Field", x=1.0),
            scene='scene2'
        ),
        go.Cone(
            x=x_sub.flatten()[b_mask], y=y_sub.flatten()[b_mask], z=z_sub.flatten()[b_mask],
            u=b_Fx[b_mask], v=b_Fy[b_mask], w=b_Fz[b_mask],
            colorscale='Reds', reversescale=False,
            sizemode='scaled', sizeref=0.5, anchor='center',
            showscale=False,
            scene='scene2'
        )
    ],
    layout=go.Layout(
        template="plotly_dark",
        title=f"EM Vector Field Propagation: ( {nx} x {ny} x {nz} )",
        scene=dict(
            aspectmode='cube',
            domain=dict(x=[0, 0.45], y=[0, 1])
        ),
        scene2=dict(
            aspectmode='cube',
            domain=dict(x=[0.55, 1], y=[0, 1])
        ),
        updatemenus=[dict(
            type="buttons",
            x=0.5, xanchor='center',
            buttons=[
                dict(label="Play", method="animate", 
                     args=[None, {"frame": {"duration": 50, "redraw": True}, "fromcurrent": True}]),
                dict(label="Pause", method="animate", args=[[None], {"mode": "immediate"}])
            ]
        )],
        annotations=[
            dict(text="Electric Field (E)", x=0.22, y=1.05, xref="paper", yref="paper", showarrow=False, font=dict(size=14)),
            dict(text="Magnetic Field (B)", x=0.78, y=1.05, xref="paper", yref="paper", showarrow=False, font=dict(size=14))
        ]
    ),
    frames=frames
)

fig.show()