#!/usr/bin/env python3
"""
Analyze skyrmion trajectory from LAMMPS dump file.
Computes center of mass of spin texture and velocity.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import center_of_mass

def read_lammpstrj(filename):
    """Read LAMMPS trajectory dump file and return timesteps and atom data."""
    timesteps = []
    boxes = []
    atoms_data = []  # list of arrays: [type, x, y, z, spx, spy, spz, sp]
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        if 'ITEM: TIMESTEP' in lines[i]:
            timestep = int(lines[i+1].strip())
            timesteps.append(timestep)
            i += 2
            
            # Skip ITEM: NUMBER OF ATOMS
            i += 1
            natoms = int(lines[i+1].strip())
            i += 2
            
            # Skip ITEM: BOX BOUNDS
            i += 1
            box_bounds = []
            for _ in range(3):
                bounds = list(map(float, lines[i].strip().split()))
                box_bounds.append(bounds)
                i += 1
            boxes.append(box_bounds)
            
            # Skip ITEM: ATOMS
            i += 1
            
            # Read atom data
            atom_array = []
            for _ in range(natoms):
                data = list(map(float, lines[i].strip().split()))
                # data: [type, x, y, z, spx, spy, spz, sp]
                # type is int but we keep as float for uniform array
                atom_array.append(data)
                i += 1
            
            atoms_data.append(np.array(atom_array))
        else:
            i += 1
    
    return timesteps, boxes, atoms_data

def compute_skyrmion_center(spin_config, r_threshold=0.5):
    """
    Compute skyrmion center as center of mass of (1 - Sz) texture.
    
    Atoms with Sz < r_threshold (inverted spins) mark skyrmion.
    r_threshold: threshold on Sz to identify skyrmion region
    """
    # Extract Sz (column 6) and x,y positions (columns 1,2)
    x = spin_config[:, 1]
    y = spin_config[:, 2]
    spz = spin_config[:, 6]
    
    # Skyrmion texture: inverted spins have Sz ~ -1, so (1-Sz) is large
    intensity = np.maximum(0, 1.0 - spz)  # 0 for aligned, ~2 for anti-aligned
    
    # Normalize
    if np.sum(intensity) > 0:
        cx = np.sum(x * intensity) / np.sum(intensity)
        cy = np.sum(y * intensity) / np.sum(intensity)
    else:
        cx, cy = np.nan, np.nan
    
    return cx, cy

def analyze_trajectory(filename):
    """Read dump and compute skyrmion trajectory."""
    print(f"Reading {filename}...")
    timesteps, boxes, atoms_data = read_lammpstrj(filename)
    
    print(f"Found {len(timesteps)} timesteps")
    
    # Compute skyrmion center for each timestep
    centers = []
    times = []
    
    for t, atoms in zip(timesteps, atoms_data):
        cx, cy = compute_skyrmion_center(atoms)
        if not np.isnan(cx):
            centers.append([cx, cy])
            times.append(t)
    
    centers = np.array(centers)
    times = np.array(times)
    
    # Compute velocity and displacement
    if len(centers) > 1:
        dt = (times[1] - times[0]) * 0.0001  # timestep in ps (metal units)
        displacements = np.diff(centers, axis=0)
        distances = np.linalg.norm(displacements, axis=1)
        velocities = distances / dt  # velocity in Angstrom/ps
        
        avg_velocity = np.mean(velocities)
        total_displacement = np.linalg.norm(centers[-1] - centers[0])
        
        print("\n=== Skyrmion Motion Analysis ===")
        print(f"Initial center: ({centers[0,0]:.2f}, {centers[0,1]:.2f})")
        print(f"Final center:   ({centers[-1,0]:.2f}, {centers[-1,1]:.2f})")
        print(f"Total displacement: {total_displacement:.3f} Angstrom")
        print(f"Average velocity: {avg_velocity:.4f} Angstrom/ps")
        print(f"Average velocity: {avg_velocity * 100:.2f} m/s")  # convert to m/s (rough)
        
        # Direction of motion
        if total_displacement > 0:
            direction = (centers[-1] - centers[0]) / total_displacement
            print(f"Direction of motion: ({direction[0]:.3f}, {direction[1]:.3f})")
        
        # Plot trajectory
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Trajectory plot
        ax1.plot(centers[:, 0], centers[:, 1], 'b-o', linewidth=2, markersize=3)
        ax1.plot(centers[0, 0], centers[0, 1], 'go', markersize=8, label='Start')
        ax1.plot(centers[-1, 0], centers[-1, 1], 'r*', markersize=15, label='End')
        ax1.set_xlabel('X (Angstrom)')
        ax1.set_ylabel('Y (Angstrom)')
        ax1.set_title('Skyrmion Trajectory (SOT-driven)')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        ax1.axis('equal')
        
        # Velocity plot
        ax2.plot(times[1:], velocities, 'b-o', linewidth=2, markersize=3)
        ax2.axhline(y=avg_velocity, color='r', linestyle='--', label=f'Avg: {avg_velocity:.4f} Å/ps')
        ax2.set_xlabel('Timestep')
        ax2.set_ylabel('Velocity (Angstrom/ps)')
        ax2.set_title('Skyrmion Velocity vs Time')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        
        plt.tight_layout()
        plt.savefig('skyrmion_analysis.png', dpi=100)
        print("\nPlot saved to: skyrmion_analysis.png")
        plt.show()
        
        return centers, times, velocities
    else:
        print("ERROR: Not enough timesteps to compute velocity")
        return None, None, None

if __name__ == '__main__':
    import sys
    
    dumpfile = 'dump.skyrmion.lammpstrj'
    if len(sys.argv) > 1:
        dumpfile = sys.argv[1]
    
    centers, times, velocities = analyze_trajectory(dumpfile)
