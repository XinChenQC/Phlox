# Phlox Quick Start Guide

Get started with Phlox in 5 minutes!

---

## Installation

```bash
cd /home/xchen/Work/ReaxANA/Phlox

# Using uv (recommended)
uv venv
uv pip install -e ".[dev]"

# Or using pip
pip install -e ".[dev]"
```

---

## Your First Analysis

### Step 1: Prepare Your Trajectory

Phlox supports:
- **XYZ format** - Standard XYZ trajectory files
- **LAMMPS dump** - LAMMPS dump files

Example XYZ format:
```
57
Lattice="18.0 0.0 0.0 0.0 18.0 0.0 0.0 0.0 25.0"
C  1.234  5.678  9.012
H  2.345  6.789  0.123
...
```

### Step 2: Create Analysis Script

```python
from phlox import Analyzer, Config

# Configure analysis
config = Config(
    n_cores=1,                    # Number of CPU cores (use 1 for now)
    traj_type="xyz",              # Or "lammps"
    traj_files=["trajectory.xyz"], # List of trajectory files
    pbc_box=((0, 0, 0), (18, 18, 25)), # Periodic boundary box
    output_dir="./output",        # Where to save results
    stability_lag=40,             # Lifetime threshold (frames)
)

# Run analysis
analyzer = Analyzer(config)
results = analyzer.run()

print(f"Analysis complete!")
print(f"Nodes: {results.network.number_of_nodes()}")
print(f"Edges: {results.network.number_of_edges()}")
print(f"Reactions: {len(results.reactions)}")
```

### Step 3: Run It!

```bash
python my_analysis.py
```

### Step 4: Check Results

```bash
ls output/
# reaction_network.dot  - Network graph file
# stats.txt             - Statistics summary
```

---

## Summary

```python
# That's it! Just 3 steps:

# 1. Configure
config = Config(traj_type="xyz", traj_files=["data.xyz"], pbc_box=((0,0,0), (50,50,50)))

# 2. Analyze
analyzer = Analyzer(config)
results = analyzer.run()

# 3. Check results
ls output/
```

**Happy analyzing!** 🚀
