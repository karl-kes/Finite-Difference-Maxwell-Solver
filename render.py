import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import glob
import os

e_files = sorted(glob.glob("output/E*.csv"), 
                 key=lambda x: int(os.path.basename(x).replace('E', '').replace('.csv', '')))
b_files = sorted(glob.glob("output/B*.csv"), 
                 key=lambda x: int(os.path.basename(x).replace('B', '').replace('.csv', '')))

e_data = [np.loadtxt(f, delimiter=',') for f in e_files]
b_data = [np.loadtxt(f, delimiter=',') for f in b_files]

e_limit = max(d.max() for d in e_data)
b_limit = max(d.max() for d in b_data)

n_frames = len(e_data)

frames = []
for i in range(n_frames):
    frames.append(go.Frame(data=[
        go.Surface(z=e_data[i], cmin=0, cmax=e_limit, colorscale='Inferno', showscale=True),
        go.Surface(z=b_data[i], cmin=0, cmax=b_limit, colorscale='Viridis', showscale=True)
    ], name=str(i)))

fig = make_subplots(
    rows=1, cols=2,
    specs=[[{'type': 'surface'}, {'type': 'surface'}]],
    subplot_titles=['|E| Electric Field', '|B| Magnetic Field']
)

fig.add_trace(go.Surface(z=e_data[0], cmin=0, cmax=e_limit, colorscale='Inferno'), row=1, col=1)
fig.add_trace(go.Surface(z=b_data[0], cmin=0, cmax=b_limit, colorscale='Viridis'), row=1, col=2)

fig.update_layout(
    updatemenus=[dict(
        type="buttons",
        buttons=[dict(label="Play", method="animate", 
                      args=[None, {"frame": {"duration": 100}, "transition": {"duration": 0}}])]
    )]
)

for scene in ['scene', 'scene2']:
    limit = e_limit if scene == 'scene' else b_limit
    fig.update_layout(**{
        scene: dict(
            xaxis=dict(autorange=False, range=[0, e_data[0].shape[1]]),
            yaxis=dict(autorange=False, range=[0, e_data[0].shape[0]]),
            zaxis=dict(autorange=False, range=[0, limit]),
            aspectmode='manual',
            aspectratio=dict(x=1, y=1, z=0.5)
        )
    })

fig.frames = frames
fig.show()