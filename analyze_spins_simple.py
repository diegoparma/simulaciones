#!/usr/bin/env python3
"""
Analyze spin dynamics from LAMMPS dump file.
Computes average spin magnetization over time to detect SOT effect.
"""

import numpy as np

def read_lammpstrj(filename):
    """Read LAMMPS trajectory dump file."""
    timesteps = []
    atoms_data = []
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        if 'ITEM: TIMESTEP' in lines[i]:
            timestep = int(lines[i+1].strip())
            timesteps.append(timestep)
            i += 2
            
            # Read NUMBER OF ATOMS
            if 'ITEM: NUMBER OF ATOMS' in lines[i]:
                i += 1
                natoms = int(lines[i].strip())
                i += 1
            
            # Skip BOX BOUNDS (3 lines)
            if 'ITEM: BOX BOUNDS' in lines[i]:
                i += 1
                for _ in range(3):
                    i += 1
            
            # Skip ITEM: ATOMS line
            if 'ITEM: ATOMS' in lines[i]:
                i += 1
            
            # Read atom data: [type, x, y, z, spx, spy, spz]
            atom_array = []
            for _ in range(natoms):
                if i >= len(lines):
                    break
                data = list(map(float, lines[i].strip().split()))
                atom_array.append(data)
                i += 1
            
            atoms_data.append(np.array(atom_array))
        else:
            i += 1
    
    return timesteps, atoms_data

def analyze_spins(atoms_data):
    """Compute average spin magnetization for each frame."""
    results = []
    
    for frame_idx, atoms in enumerate(atoms_data):
        # Extract spin components (columns 4, 5, 6 = spx, spy, spz)
        spx = atoms[:, 4]
        spy = atoms[:, 5]
        spz = atoms[:, 6]
        
        # Average spin
        avg_spx = np.mean(spx)
        avg_spy = np.mean(spy)
        avg_spz = np.mean(spz)
        avg_mag = np.sqrt(avg_spx**2 + avg_spy**2 + avg_spz**2)
        
        # Spin variance (indicator of spin disorder)
        var_sp = np.var(np.sqrt(spx**2 + spy**2 + spz**2))
        
        results.append({
            'frame': frame_idx,
            'avg_spx': avg_spx,
            'avg_spy': avg_spy,
            'avg_spz': avg_spz,
            'avg_mag': avg_mag,
            'var_sp': var_sp
        })
    
    return results

def main():
    import sys
    
    dumpfile = 'dump.spins.lammpstrj'
    if len(sys.argv) > 1:
        dumpfile = sys.argv[1]
    
    print(f"Reading {dumpfile}...")
    timesteps, atoms_data = read_lammpstrj(dumpfile)
    
    print(f"Found {len(timesteps)} frames\n")
    
    # Analyze spin dynamics
    results = analyze_spins(atoms_data)
    
    print("Frame Index | Time (fs) | Avg Sx   | Avg Sy   | Avg Sz   | Avg |S| | Variance")
    print("-" * 85)
    
    for i, res in enumerate(results):
        time_fs = timesteps[i] * 0.0001 * 1000  # metal timestep -> fs
        print(f"{i:5d}       | {time_fs:8.2f}  | {res['avg_spx']:+.5f}  | {res['avg_spy']:+.5f}  | {res['avg_spz']:+.5f}  | {res['avg_mag']:.5f}  | {res['var_sp']:.5f}")
    
    # Detect changes
    print("\n=== ANALYSIS ===")
    
    first_half = results[:len(results)//2]
    second_half = results[len(results)//2:]
    
    avg_mag_first = np.mean([r['avg_mag'] for r in first_half])
    avg_mag_second = np.mean([r['avg_mag'] for r in second_half])
    
    avg_sx_first = np.mean([r['avg_spx'] for r in first_half])
    avg_sx_second = np.mean([r['avg_spx'] for r in second_half])
    
    avg_sy_first = np.mean([r['avg_spy'] for r in first_half])
    avg_sy_second = np.mean([r['avg_spy'] for r in second_half])
    
    print(f"First half avg magnetization:  |S| = {avg_mag_first:.5f}")
    print(f"Second half avg magnetization: |S| = {avg_mag_second:.5f}")
    print(f"Change in |S|: {(avg_mag_second - avg_mag_first):.6f}")
    
    print(f"\nFirst half avg Sx: {avg_sx_first:+.5f}")
    print(f"Second half avg Sx: {avg_sx_second:+.5f}")
    print(f"Change in Sx: {(avg_sx_second - avg_sx_first):+.6f}")
    
    print(f"\nFirst half avg Sy: {avg_sy_first:+.5f}")
    print(f"Second half avg Sy: {avg_sy_second:+.5f}")
    print(f"Change in Sy: {(avg_sy_second - avg_sy_first):+.6f}")
    
    # Detect SOT effect
    if abs(avg_sx_second - avg_sx_first) > 0.01 or abs(avg_sy_second - avg_sy_first) > 0.01:
        print("\n✓ SOT EFFECT DETECTED: Significant change in spin polarization")
    else:
        print("\n⚠ Minimal spin polarization change (spins remain randomized)")

if __name__ == '__main__':
    main()
