import numpy as np
import plotly.graph_objects as go
import glob
import os

files = sorted(glob.glob("output/output*.csv"), key=lambda x: int(os.path.basename(x).replace('output', '').replace('.csv', '')))

all_data = [np.loadtxt(f, delimiter=',') for f in files]
vmin = min(d.min() for d in all_data)
vmax = max(d.max() for d in all_data)
vlimit =  0.8 * max(abs(vmin), abs(vmax))

frames = []
for i, data in enumerate(all_data):
    frames.append(go.Frame(data=[go.Surface(z=data, cmin=-vlimit, cmax=vlimit, colorscale='RdBu')], name=str(i)))

fig = go.Figure(
    data=[go.Surface(z=all_data[0], cmin=-vlimit, cmax=vlimit, colorscale='RdBu')],
    frames=frames
)

fig.update_layout(
    scene=dict(
        xaxis=dict(range=[0, all_data[0].shape[1]], autorange=False),
        yaxis=dict(range=[0, all_data[0].shape[0]], autorange=False),
        zaxis=dict(range=[-vlimit, vlimit], autorange=False),
        aspectmode='manual',
        aspectratio=dict(x=1, y=1, z=0.5)
    ),
    updatemenus=[dict(
        type="buttons",
        buttons=[dict(label="Play", method="animate", args=[None, {"frame": {"duration": 100}, "transition": {"duration": 0}}])]
    )]
)

fig.show()