# Skyrmion SOT Playground

A practical simulation workspace for spin dynamics and skyrmion-oriented experiments using a local LAMMPS fork with Spin-Orbit Torque (SOT) extensions.

[![LAMMPS](https://img.shields.io/badge/LAMMPS-local%20source-1f6feb)](lammps)
[![Language](https://img.shields.io/badge/Python-analysis%20%2B%20viz-0b7a75)](.)
[![Status](https://img.shields.io/badge/status-public%20repo-2ea043)](.)

## Why this repo exists

This project brings together three things in one place:

- A local LAMMPS codebase with SOT-related changes in the SPIN package.
- Input decks to run minimal and skyrmion-focused scenarios.
- Analysis and visualization scripts to validate and communicate results.

## Highlights

- Local LAMMPS source is included directly under `lammps/`.
- SOT implementation lives in:
  - `lammps/src/SPIN/fix_precession_spin.cpp`
  - `lammps/src/SPIN/fix_precession_spin.h`
  - `lammps/doc/src/fix_precession_spin.rst`
- Ready-to-run simulation inputs for quick validation and skyrmion experiments.
- Interactive HTML outputs and OVITO workflows.

## Repository map

```text
simulaciones/
|-- lammps/                        # Vendored LAMMPS source + local modifications
|-- in.*.lmp                       # LAMMPS input scripts
|-- analyze_*.py                   # Numeric analysis helpers
|-- generate_*.py                  # Data generation / setup helpers
|-- visualize_*.py                 # Plotly visualization helpers
|-- ovito_*.py                     # OVITO pipeline scripts
|-- dump*.lammpstrj                # Sample trajectories
|-- spin_*.html                    # Interactive 3D visualization snapshots
|-- VALIDATION_REPORT.md           # Functional validation of SOT behavior
|-- VISUALIZATION_SETUP.md         # Setup and implementation notes
|-- VISUALIZATION_GUIDE.md         # Usage guide for visualization tools
```

## Quick start

### 1) Prepare Python environment

```bash
cd /Users/diegoemilioparma/Documents/simulaciones
python3 -m venv .venv
source .venv/bin/activate
pip install -U pip numpy plotly
```

### 2) Build local LAMMPS with SPIN package

```bash
cd lammps
cmake -S cmake -B build-local-sot -D PKG_SPIN=on -D BUILD_MPI=off
cmake --build build-local-sot -j4
```

### 3) Run a reference SOT test

```bash
cd ..
./lammps/build-local-sot/lmp -in in.sot_simple.lmp > log.sot_simple.txt
```

### 4) Analyze and visualize

```bash
python3 analyze_spins_simple.py dump.spins.lammpstrj
python3 visualize_spins_plotly.py
```

Open generated files such as `spin_frame_025.html` in your browser.

## Main scripts

- `generate_skyrmion.py`: build initial atom/spin configurations.
- `analyze_skyrmion.py`: skyrmion-centric post-processing.
- `analyze_spins_simple.py`: quick spin statistics over trajectories.
- `visualize_spins_plotly.py`: interactive HTML 3D views.
- `open_in_ovito.py`: direct launch into OVITO.
- `ovito_visualize_spins.py`: programmable OVITO workflows.

## Main inputs

- `in.sot_simple.lmp`: compact SOT validation scenario.
- `in.skyrmion_init.lmp`: initialization flow.
- `in.skyrmion_sot.lmp`: skyrmion + SOT workflow.
- `in.skyrmion_test.lmp`: test-oriented skyrmion run.
- `in.test_sot_minimal.lmp`: parser/minimal behavior check.
- `in.test_sot_parsing.lmp`: SOT argument parsing check.

## Documentation index

- `VALIDATION_REPORT.md`
- `VISUALIZATION_SETUP.md`
- `VISUALIZATION_GUIDE.md`
- `docs/README.md`

## Notes on data and size

This repository includes sample dumps and generated HTML views for reproducibility.
Large local build and virtual environment folders are excluded through `.gitignore`.

## Versioning and tags

Release tags are used to mark important milestones:

- `v0.1.0`: initial public project state.
- `v0.2.0`: LAMMPS source vendored with local SOT modifications.
- `v0.3.0`: repository cleanup + polished root documentation.

## Citation and credit

- LAMMPS upstream project: https://github.com/lammps/lammps
- Local changes in this repository are focused on SOT support and workflow integration.
