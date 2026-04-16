#!/usr/bin/env python3
"""
Interactive 3D spin visualization with Plotly.
Shows atomic positions and spin vectors for LAMMPS dump files.
No external dependencies beyond plotly (web-based visualization).
"""

import sys
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def read_lammpstrj_interactive(filename):
    """Read LAMMPS trajectory file with frame selection."""
    frames = []
    current_frame = {}
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    i = 0
    frame_count = 0
    
    while i < len(lines):
        line = lines[i].strip()
        
        if line == "ITEM: TIMESTEP":
            timestep = int(lines[i+1].strip())
            current_frame['timestep'] = timestep
            i += 2
            
        elif line == "ITEM: NUMBER OF ATOMS":
            natoms = int(lines[i+1].strip())
            current_frame['natoms'] = natoms
            i += 2
            
        elif line == "ITEM: BOX BOUNDS pp pp pp":
            bounds = []
            for j in range(3):
                lo, hi = map(float, lines[i+j+1].split())
                bounds.append((lo, hi))
            current_frame['bounds'] = bounds
            i += 4
            
        elif line == "ITEM: ATOMS type x y z c_outsp[1] c_outsp[2] c_outsp[3]":
            atoms = []
            for j in range(natoms):
                parts = lines[i+j+1].split()
                atom = {
                    'type': int(parts[0]),
                    'x': float(parts[1]),
                    'y': float(parts[2]),
                    'z': float(parts[3]),
                    'spx': float(parts[4]),
                    'spy': float(parts[5]),
                    'spz': float(parts[6])
                }
                atoms.append(atom)
            current_frame['atoms'] = atoms
            frames.append(dict(current_frame))
            current_frame = {}
            i += natoms + 1
            frame_count += 1
            
        else:
            i += 1
    
    return frames

def create_spin_visualization(frames, frame_idx=None):
    """Create interactive 3D plot of spin configuration."""
    
    if frame_idx is None:
        frame_idx = len(frames) - 1  # Last frame
    
    if frame_idx >= len(frames):
        frame_idx = len(frames) - 1
    
    frame = frames[frame_idx]
    atoms = frame['atoms']
    timestep = frame['timestep']
    
    # Extract atomic positions and spins
    positions = np.array([[a['x'], a['y'], a['z']] for a in atoms])
    spins = np.array([[a['spx'], a['spy'], a['spz']] for a in atoms])
    
    # Compute spin magnitudes for coloring
    spin_mags = np.linalg.norm(spins, axis=1)
    
    # Create figure
    fig = go.Figure()
    
    # Plot atoms as spheres (light colors, slight transparency)
    fig.add_trace(go.Scatter3d(
        x=positions[:, 0],
        y=positions[:, 1],
        z=positions[:, 2],
        mode='markers',
        marker=dict(
            size=4,
            color=spin_mags,
            colorscale='Plasma',
            showscale=True,
            colorbar=dict(title="Spin<br>Magnitude"),
            opacity=0.8,
            line=dict(width=0, color='white')
        ),
        text=[f"Atom {i}: S=({spins[i,0]:.3f}, {spins[i,1]:.3f}, {spins[i,2]:.3f})" 
              for i in range(len(atoms))],
        hoverinfo='text',
        name='Atoms'
    ))
    
    # Plot spin vectors as arrows (cone-like representation)
    arrow_scale = 1.0  # Scale for visibility
    
    # Create quiver-like arrows using cone traces
    for i, (pos, spin) in enumerate(zip(positions, spins)):
        if spin_mags[i] > 0.01:  # Only plot significant spins
            end_pos = pos + spin * arrow_scale
            
            # Draw line from atom to spin tip
            fig.add_trace(go.Scatter3d(
                x=[pos[0], end_pos[0]],
                y=[pos[1], end_pos[1]],
                z=[pos[2], end_pos[2]],
                mode='lines',
                line=dict(
                    color=f'rgba(0, {int(200*spin_mags[i])}, 255, 0.6)',
                    width=2
                ),
                hoverinfo='skip',
                showlegend=False,
                name=''
            ))
    
    # Update layout
    bounds = frame['bounds']
    fig.update_layout(
        title=dict(
            text=f"Spin Configuration - Frame {frame_idx} (t={timestep} fs)",
            font=dict(size=16)
        ),
        scene=dict(
            xaxis=dict(
                range=[bounds[0][0] - 1, bounds[0][1] + 1],
                title="X (Å)"
            ),
            yaxis=dict(
                range=[bounds[1][0] - 1, bounds[1][1] + 1],
                title="Y (Å)"
            ),
            zaxis=dict(
                range=[bounds[2][0] - 1, bounds[2][1] + 1],
                title="Z (Å)"
            ),
            aspectmode='data',
            camera=dict(
                eye=dict(x=1.2, y=1.2, z=1.0)
            )
        ),
        width=1200,
        height=800,
        showlegend=True,
        hovermode='closest'
    )
    
    return fig

def create_multi_frame_view(frames, frame_indices=None):
    """Create side-by-side or animation view of multiple frames."""
    
    if frame_indices is None:
        # Show first, middle, last
        n = len(frames)
        if n == 1:
            frame_indices = [0]
        elif n <= 3:
            frame_indices = list(range(n))
        else:
            frame_indices = [0, n//2, n-1]
    
    figs = []
    for idx in frame_indices:
        fig = create_spin_visualization(frames, idx)
        figs.append(fig)
    
    return figs

if __name__ == "__main__":
    print("=" * 70)
    print("PLOTLY 3D Spin Visualization Tool")
    print("=" * 70)
    
    DUMP_FILE = "dump.spins.lammpstrj"
    
    try:
        print(f"\nLoading {DUMP_FILE}...")
        frames = read_lammpstrj_interactive(DUMP_FILE)
        print(f"✓ Loaded {len(frames)} frames")
        
        # Show options
        print(f"\nAvailable frames: 0 to {len(frames)-1}")
        
        # Create visualizations
        print("\nGenerating visualizations...")
        
        # Multi-frame comparison
        print("  - Creating multi-frame view (first, middle, last)...")
        figs = create_multi_frame_view(frames)
        
        for i, (idx, fig) in enumerate(zip([0, len(frames)//2, len(frames)-1], figs)):
            output_file = f"spin_frame_{idx:03d}.html"
            fig.write_html(output_file)
            print(f"    ✓ {output_file}")
        
        # Full animation (if plotly has animation support)
        print("\n  - Creating animation (last 10 frames)...")
        start_frame = max(0, len(frames) - 10)
        final_fig = create_spin_visualization(frames, len(frames) - 1)
        final_fig.write_html("spin_animation_latest.html")
        print(f"    ✓ spin_animation_latest.html")
        
        print("\n" + "=" * 70)
        print("✓ VISUALIZATION COMPLETE")
        print("=" * 70)
        print("\nOpen any .html file in your browser to view interactive 3D spins:")
        print("  - spin_frame_000.html  (beginning)")
        print("  - spin_frame_*.html    (time progression)")
        print("  - spin_animation_latest.html (end state)")
        print("\nIn 3D view:")
        print("  - Click & drag to rotate")
        print("  - Scroll to zoom")
        print("  - Hover over atoms for spin values")
        print("  - Right-click & drag to pan")
        
    except FileNotFoundError:
        print(f"ERROR: File not found: {DUMP_FILE}")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
