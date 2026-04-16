#!/bin/bash
# Quick launcher for spin visualizations

cd /Users/diegoemilioparma/Documents/simulaciones

echo "=========================================================="
echo "3D Spin Visualization Launcher"
echo "=========================================================="
echo ""
echo "Available visualizations:"
echo ""
echo "1. spin_frame_000.html    - Frame 0 (equilibrium)"
echo "2. spin_frame_025.html    - Frame 25 (middle)"
echo "3. spin_frame_050.html    - Frame 50 (final state)"
echo "4. spin_animation_latest.html - Final with more details"
echo "5. Open OVITO GUI with dump file"
echo "6. View analysis report"
echo ""
echo "Enter choice (1-6) or visualization name: "
read choice

case $choice in
    1) 
        echo "Opening spin_frame_000.html..."
        open spin_frame_000.html
        ;;
    2)
        echo "Opening spin_frame_025.html..."
        open spin_frame_025.html
        ;;
    3)
        echo "Opening spin_frame_050.html..."
        open spin_frame_050.html
        ;;
    4)
        echo "Opening spin_animation_latest.html..."
        open spin_animation_latest.html
        ;;
    5)
        echo "Launching OVITO..."
        python3 open_in_ovito.py
        ;;
    6)
        echo "Viewing analysis report..."
        open VALIDATION_REPORT.md
        ;;
    *)
        if [ -f "$choice" ]; then
            echo "Opening $choice..."
            open "$choice"
        else
            echo "Invalid choice"
        fi
        ;;
esac
