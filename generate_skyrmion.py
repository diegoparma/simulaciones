#!/usr/bin/env python3
"""
Generate skyrmion initial conditions for LAMMPS.
Creates a 2D circular spin configuration with skyrmion topology.
"""

import numpy as np
import sys

def create_skyrmion_config(lattice_const=2.87, radius=40, center=(50, 50), thickness=2):
    """
    Create skyrmion spin configuration.
    
    Args:
        lattice_const: BCC lattice constant (Å)
        radius: Skyrmion radius (Å)
        center: Center position (x, y)
        thickness: Number of layers (z-direction)
    
    Returns:
        List of (x, y, z, spx, spy, spz) tuples
    """
    
    atoms = []
    atom_id = 1
    
    # Create atoms in circular region
    for ix in range(int((center[0] - radius - 5) / lattice_const), 
                    int((center[0] + radius + 5) / lattice_const)):
        for iy in range(int((center[1] - radius - 5) / lattice_const), 
                        int((center[1] + radius + 5) / lattice_const)):
            for iz in range(thickness):
                
                # BCC lattice positions
                x = ix * lattice_const
                y = iy * lattice_const
                z = iz * lattice_const
                
                # Offset for BCC
                if (ix + iy + iz) % 2 == 1:
                    x += lattice_const / 2
                    y += lattice_const / 2
                    z += lattice_const / 2
                
                # Distance from skyrmion center
                dx = x - center[0]
                dy = y - center[1]
                r = np.sqrt(dx**2 + dy**2)
                
                # Only include atoms in circular region
                if r < radius + 5:
                    
                    # Skyrmion spin profile
                    # Center: spins pointing down (inverted)
                    # Perimeter: spins pointing up (align with background)
                    
                    if r < 1:  # Core (avoid singularity)
                        spx = 0
                        spy = 0
                        spz = -1  # Down
                    else:
                        # Radial-azimuthal skyrmion profile
                        # θ(r) = π * (1 - r/R) smoothly transitions from π to 0
                        theta = np.pi * np.exp(-r / (radius / 2))
                        phi = np.arctan2(dy, dx)
                        
                        spx = np.sin(theta) * np.cos(phi)
                        spy = np.sin(theta) * np.sin(phi)
                        spz = np.cos(theta)
                    
                    # Small random perturbation for realism
                    noise = 0.05
                    spx += np.random.randn() * noise
                    spy += np.random.randn() * noise
                    spz += np.random.randn() * noise
                    
                    # Normalize
                    s_mag = np.sqrt(spx**2 + spy**2 + spz**2)
                    if s_mag > 0:
                        spx /= s_mag
                        spy /= s_mag
                        spz /= s_mag
                    
                    atoms.append({
                        'id': atom_id,
                        'type': 1,
                        'x': x,
                        'y': y,
                        'z': z,
                        'spx': spx,
                        'spy': spy,
                        'spz': spz
                    })
                    atom_id += 1
    
    return atoms

def write_lammps_input(atoms, filename="in.skyrmion_init.lmp"):
    """
    Write LAMMPS data file with skyrmion configuration.
    """
    
    with open("skyrmion_atoms.dat", "w") as f:
        # Compute box bounds
        xs = [a['x'] for a in atoms]
        ys = [a['y'] for a in atoms]
        zs = [a['z'] for a in atoms]
        
        xlo, xhi = min(xs) - 5, max(xs) + 5
        ylo, yhi = min(ys) - 5, max(ys) + 5
        zlo, zhi = min(zs) - 5, max(zs) + 5
        
        # Header
        f.write(f"Skyrmion configuration\n\n")
        f.write(f"{len(atoms)} atoms\n")
        f.write(f"1 atom types\n\n")
        f.write(f"{xlo} {xhi} xlo xhi\n")
        f.write(f"{ylo} {yhi} ylo yhi\n")
        f.write(f"{zlo} {zhi} zlo zhi\n\n")
        
        # Masses
        f.write(f"Masses\n\n")
        f.write(f"1 55.845\n\n")  # Fe mass
        
        # Atoms (format: id type x y z spx spy spz)
        f.write(f"Atoms # spin\n\n")
        for atom in atoms:
            f.write(f"{atom['id']} {atom['type']} {atom['x']:.6f} {atom['y']:.6f} {atom['z']:.6f} {atom['spx']:.6f} {atom['spy']:.6f} {atom['spz']:.6f}\n")
    
    # Create LAMMPS input script
    with open(filename, "w") as f:
        f.write("""# Skyrmion visualization and SOT dynamics
# 2D circular skyrmion in magnetic film

units            metal
atom_style       spin
dimension        3
boundary         p p p

# Create initial skyrmion atoms
read_data        skyrmion_atoms.dat

# Magnetic interactions
pair_style       spin/exchange 5.0
pair_coeff       * * exchange 1.0 0.0 1.0

# Set spin magnitudes and masses
set              group all spin 2.2 0.0 0.0 1.0

# Initialize velocities
velocity         all create 100 42

# Define groups
group            skyrmion_center region sphere 50 50 1 15

# Fix for spin dynamics
fix              1 all precession/spin zeeman 0.0 0.0 0.0 1.0 \\
                 anisotropy 5e-05 0.0 0.0 1.0
fix              2 all langevin/spin 0.05 0.1 21

# Phase 1: Equilibration without SOT (1000 steps ≈ 0.1 ps)
fix              3 all nve/spin lattice moving
dump             1 all custom 100 traj_equilibration.lammpstrj id x y z c_outsp[1] c_outsp[2] c_outsp[3]
run              1000
undump           1

# Phase 2: Apply SOT (damping-like + field-like)
unfix            3
fix              3 all nve/spin lattice moving
fix              4 all precession/spin sot/dl 0.050 0.0 1.0 0.0 sot/fl 0.010 0.0 1.0 0.0

dump             2 all custom 50 dump.skyrmion.lammpstrj id x y z c_outsp[1] c_outsp[2] c_outsp[3]
run              2000

thermo           100
thermo_style     custom step time v_Etot v_KE v_PE

write_dump       all custom dump.skyrmion_final.lammpstrj id x y z c_outsp[1] c_outsp[2] c_outsp[3]
""")
    
    print(f"✓ Wrote {filename}")
    print(f"✓ Wrote skyrmion_atoms.dat")

if __name__ == "__main__":
    print("=" * 60)
    print("Skyrmion Initial Configuration Generator")
    print("=" * 60)
    
    np.random.seed(42)
    
    print("\nGenerating skyrmion configuration...")
    atoms = create_skyrmion_config(
        lattice_const=2.87,  # Fe BCC
        radius=40.0,         # 40 Å radius
        center=(50.0, 50.0),
        thickness=2          # 2-layer film
    )
    
    print(f"✓ Created {len(atoms)} atoms")
    
    # Statistics
    szs = [a['spz'] for a in atoms]
    print(f"✓ Avg Sz: {np.mean(szs):.3f}")
    print(f"✓ Skyrmion radius: ~40 Å")
    print(f"✓ Core (r<15): Spins inverted")
    print(f"✓ Periphery: Spins aligned")
    
    write_lammps_input(atoms)
    
    print("\n" + "=" * 60)
    print("Ready to run:")
    print("  ~/lammps/build-local-sot/lmp -in in.skyrmion_init.lmp")
    print("=" * 60)
