#!/usr/bin/env python3
"""Generate circular skyrmion visualization from LAMMPS dump."""

import numpy as np
import plotly.graph_objects as go

def read_dump(filename):
    frames = []
    current = {}
    with open(filename) as f:
        lines = f.readlines()
    i = 0
    natoms = 0
    while i < len(lines):
        line = lines[i].strip()
        if line == "ITEM: TIMESTEP":
            current['ts'] = int(lines[i+1].strip()); i += 2
        elif line == "ITEM: NUMBER OF ATOMS":
            natoms = int(lines[i+1].strip()); i += 2
        elif line.startswith("ITEM: BOX BOUNDS"):
            i += 4
        elif line.startswith("ITEM: ATOMS"):
            atoms = []
            for j in range(natoms):
                p = lines[i+j+1].split()
                atoms.append([float(x) for x in p])
            current['atoms'] = np.array(atoms)
            frames.append(dict(current)); current = {}
            i += natoms + 1
        else:
            i += 1
    return frames

print("Reading dump.skyrmion.lammpstrj ...")
frames = read_dump("dump.skyrmion.lammpstrj")
print(f"Loaded {len(frames)} frames")

frame = frames[-1]
atoms = frame['atoms']
x, y, z = atoms[:, 1], atoms[:, 2], atoms[:, 3]
spx, spy, spz = atoms[:, 4], atoms[:, 5], atoms[:, 6]

# Top-view colored by Sz (out-of-plane component)
fig = go.Figure()
fig.add_trace(go.Scatter(
    x=x, y=y,
    mode='markers',
    marker=dict(
        size=8,
        color=spz,
        colorscale='RdBu',
        cmin=-1, cmax=1,
        showscale=True,
        colorbar=dict(title="Sz (↑ out-of-plane ↓"),
        opacity=0.9,
        line=dict(width=0.2, color='gray')
    ),
    text=[f"S=({spx[i]:.2f}, {spy[i]:.2f}, {spz[i]:.2f})" for i in range(len(atoms))],
    hoverinfo='text+x+y',
    name='Spins'
))

# Add in-plane arrows every 8th atom
scale = 1.5
for i in range(0, len(atoms), 8):
    sx_n = spx[i] * scale
    sy_n = spy[i] * scale
    if abs(sx_n) + abs(sy_n) > 0.05:
        fig.add_annotation(
            x=x[i] + sx_n, y=y[i] + sy_n,
            ax=x[i], ay=y[i],
            xref='x', yref='y', axref='x', ayref='y',
            showarrow=True,
            arrowhead=2, arrowsize=1.0, arrowwidth=0.8,
            arrowcolor='rgba(20,20,20,0.3)',
        )

fig.update_layout(
    title=dict(
        text=f"Circular Skyrmion - Top View ({len(atoms)} atoms)<br>"
             f"<sub>Color = Sz (red=up, blue=down) | Arrows = in-plane spin direction</sub>",
        font=dict(size=16)
    ),
    xaxis=dict(title="X (Å)", scaleanchor="y", scaleratio=1),
    yaxis=dict(title="Y (Å)"),
    width=900, height=850,
    plot_bgcolor='white',
    paper_bgcolor='white'
)

fig.write_html("skyrmion_topview.html")
print("✓ skyrmion_topview.html")

# 3D spin landscape (height = Sz)
fig3d = go.Figure()
fig3d.add_trace(go.Scatter3d(
    x=x, y=y, z=spz,
    mode='markers',
    marker=dict(
        size=4,
        color=spz,
        colorscale='RdBu',
        cmin=-1, cmax=1,
        showscale=True,
        colorbar=dict(title="Sz"),
        opacity=0.9
    ),
    text=[f"({x[i]:.0f},{y[i]:.0f}): S=({spx[i]:.2f},{spy[i]:.2f},{spz[i]:.2f})" 
          for i in range(len(atoms))],
    hoverinfo='text'
))
fig3d.update_layout(
    title="Skyrmion Spin Landscape (x,y = atom position, z = Sz)",
    scene=dict(
        xaxis_title="X (Å)", yaxis_title="Y (Å)", zaxis_title="Sz",
        aspectmode='data',
        camera=dict(eye=dict(x=1.3, y=-1.3, z=1.1))
    ),
    width=1000, height=700
)
fig3d.write_html("skyrmion_3d_landscape.html")
print("✓ skyrmion_3d_landscape.html")
print("\nDone! Open .html files in your browser.")
