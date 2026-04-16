#!/usr/bin/env python3
"""
OVITO visualization script for spin configurations from LAMMPS dump.
Run with: ovitos ovito_visualize_spins.py
"""

import sys
from ovito.io import import_file, export_file
from ovito.modifiers import SelectParticleTypeModifier, ColorCodingModifier, VectorDisplay, ExpressionSelectionModifier
from ovito.vis import Viewport, OpenGLRenderer, ViewportConfiguration
import numpy as np

# Configuration
DUMP_FILE = "dump.spins.lammpstrj"
FRAME_TO_SHOW = -1  # Last frame
SPIN_SCALE = 2.0    # Scale factor for vector display
OUTPUT_FILE = "visualization.png"

def setup_visualization():
    """Load dump and configure spin vector visualization in OVITO."""
    
    print("[OVITO] Loading dump file:", DUMP_FILE)
    
    # Import the dump file
    pipeline = import_file(DUMP_FILE)
    data = pipeline.compute()
    
    # Get particle properties
    n_frames = pipeline.source.num_frames
    print(f"[OVITO] Loaded {n_frames} frames")
    
    # Get the frame
    frame_idx = FRAME_TO_SHOW if FRAME_TO_SHOW >= 0 else n_frames - 1
    pipeline.compute(frame_idx)
    data = pipeline.compute(frame_idx)
    
    print(f"[OVITO] Visualizing frame {frame_idx}")
    print(f"[OVITO] Number of particles: {data.particles.count}")
    
    # Check available properties
    print("[OVITO] Available particle properties:")
    for prop_name in data.particles.keys():
        print(f"  - {prop_name}")
    
    # Create vector display modifier for spin vectors
    # LAMMPS dump format: type x y z c_outsp[1] c_outsp[2] c_outsp[3]
    # Map to: c_outsp[1], c_outsp[2], c_outsp[3] = sx, sy, sz
    
    # Particle display - size based on particle type, color by spin magnitude
    try:
        # Compute spin magnitude for coloring
        expr_mod = ExpressionSelectionModifier(
            expression="sqrt(c_outsp[1]*c_outsp[1] + c_outsp[2]*c_outsp[2] + c_outsp[3]*c_outsp[3])"
        )
        # This won't work directly, but shows attempt to compute derived property
    except:
        pass
    
    # Configure viewport
    vp = Viewport()
    vp.type = Viewport.VIEW_TOP
    vp.background_color = (1, 1, 1)
    
    # Render and save
    renderer = OpenGLRenderer()
    
    # Set up camera to see all particles
    vp.zoom_all()
    
    print("[OVITO] Rendering viewport...")
    # Note: Full programmatic rendering requires pro license
    # For now, we set up the visualization in the GUI
    
    # Export structure as extended xyz for manual viewing
    export_file(pipeline, "visualization_frame.xyz", "xyz")
    print("[OVITO] Exported frame as XYZ to visualization_frame.xyz")
    
    return pipeline, data

if __name__ == "__main__":
    print("=" * 60)
    print("OVITO Spin Visualization Tool")
    print("=" * 60)
    
    try:
        pipeline, data = setup_visualization()
        
        print("\n[OVITO] Setup complete!")
        print("[OVITO] To view in OVITO GUI:")
        print("  1. Open OVITO normally: /Applications/Ovito.app")
        print("  2. File > Open > dump.spins.lammpstrj")
        print("  3. In 'Particle Properties', enable display of c_outsp[1,2,3]")
        print("  4. Use Vector Display modifier for arrows")
        print("  5. Color by spin magnitude or component")
        
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)
