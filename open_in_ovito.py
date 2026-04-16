#!/usr/bin/env python3
"""
Script to prepare LAMMPS dump for OVITO visualization and launch OVITO.
Uses OVITO's native capabilities to visualize spin vectors.
"""

import subprocess
import os
import sys

def open_in_ovito(dump_file):
    """Open LAMMPS dump file in OVITO GUI application."""
    
    if not os.path.exists(dump_file):
        print(f"ERROR: File not found: {dump_file}")
        return False
    
    # Path to OVITO application
    ovito_app = "/Applications/Ovito.app/Contents/MacOS/ovito"
    
    if not os.path.exists(ovito_app):
        print("ERROR: OVITO not found at /Applications/Ovito.app")
        print("Please install OVITO with: brew install ovito")
        return False
    
    print(f"Opening {dump_file} in OVITO...")
    print("(This will launch the OVITO GUI window)")
    print()
    
    try:
        # Launch OVITO with the dump file
        subprocess.Popen([ovito_app, dump_file])
        print("✓ OVITO launched!")
        print()
        print("VISUALIZATION TIPS:")
        print("-" * 60)
        print("1. In the PIPELINE panel (left), you'll see:")
        print("   - Data Source: LAMMPS Dump File")
        print("   - Particle properties: type, x, y, z, c_outsp[1], c_outsp[2], c_outsp[3]")
        print()
        print("2. To visualize spin vectors:")
        print("   a) Go to 'Add Modifier' → search 'Vector'")
        print("   b) Add 'Vector Display' modifier")
        print("   c) In Vector Display settings:")
        print("      - Property: c_outsp[1] or create custom vector")
        print("      - Vector mapping: [c_outsp[1], c_outsp[2], c_outsp[3]]")
        print()
        print("3. For better visibility:")
        print("   - Increase vector length/width")
        print("   - Color by spin magnitude or component")
        print("   - Particle size: affect rendering quality")
        print()
        print("4. Playback:")
        print("   - Use timeline at bottom to step through frames")
        print("   - Scrubber shows time progression")
        print()
        print("5. Camera control:")
        print("   - Rotate: Left mouse button + drag")
        print("   - Zoom: Scroll wheel")
        print("   - Pan: Middle mouse button (or Cmd+Click on Mac) + drag")
        print()
        print("📌 For advanced visualization, use the Python console or")
        print("   the ovito_visualize_spins.py script with ovitos command")
        print("-" * 60)
        
        return True
        
    except Exception as e:
        print(f"ERROR: Failed to launch OVITO: {e}")
        return False

if __name__ == "__main__":
    
    dump_file = "dump.spins.lammpstrj"
    
    if len(sys.argv) > 1:
        dump_file = sys.argv[1]
    
    print("=" * 60)
    print("OVITO Launcher for Spin Configuration Visualization")
    print("=" * 60)
    print()
    
    success = open_in_ovito(dump_file)
    
    if not success:
        sys.exit(1)
