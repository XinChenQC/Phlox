# Phlox

**Reaction mechanism discovery from reactive molecular dynamics trajectories**

Phlox (formerly ReaxANA/CPX-MechGen) analyzes reactive molecular dynamics simulations to automatically discover reaction mechanisms by tracking molecular fragmentation, building reaction networks, and classifying reaction pathways.

## Development Status

**Current Implementation Stage:**
- ✅ **Fragment Analysis**: Complete - Molecular fragmentation detection and tracking fully implemented
- ⚠️ **Network Filters**: In Progress - Oscillation removal and node contraction currently disabled
- 🚧 **Reaction Classification**: Partial - Core classification framework in place, some methods pending

## Features

- 🔬 **Trajectory Analysis**: Process XYZ and LAMMPS trajectory files
- 🧬 **Molecular Tracking**: Automatic fragmentation detection with PBC support
- 📊 **Reaction Networks**: Build directed graphs of chemical reactions
- 🚀 **Parallel Processing**: Multi-core support for large trajectories
- ⚗️ **Catalyst Support**: Track catalytic species separately
- 📈 **Time Series**: Species abundance evolution over simulation time
- 🎨 **Visualization**: Generate molecular structures and reaction graphs

## Installation

### From PyPI (when published)
```bash
pip install phlox
```

### From source
```bash
git clone https://github.com/chenxin/phlox.git
cd phlox
pip install -e .
```

### With development dependencies
```bash
pip install -e ".[dev]"
```

### Using Docker
```bash
docker build -t phlox:latest -f docker/Dockerfile .
docker run -v $(pwd)/data:/data phlox:latest phlox analyze /data/input.conf
```

## Quick Start

### 1. Prepare input configuration

Create `input.conf`:
```
ncores          4
trajtype        lmp
trajfiles       trajectory.trj
molelifttime    1 ps
timeperframe    0.1 fs
pbcmole         True
pbcbox          0 100   0 100   0 100
lmpelement      C H O Ni
```

### 2. Run analysis

```bash
phlox analyze input.conf
```

### 3. Check results

- `result.log`: Species list, time series, and reaction events
- `data.json`: Structured data for downstream analysis
- `specRec/`: XYZ and SVG files for each species
- `*.dot`: GraphViz reaction network files

## Python API

```python
from phlox import Analyzer
from phlox.io import read_config

# Load configuration
config = read_config("input.conf")

# Create analyzer
analyzer = Analyzer(config)

# Process trajectory
results = analyzer.run()

# Access results
print(f"Found {len(results.species)} unique species")
print(f"Detected {len(results.reactions)} reaction types")

# Export results
results.to_json("output.json")
results.plot_network("network.png")
```

## Documentation

Full documentation is available at [https://phlox.readthedocs.io](https://phlox.readthedocs.io)

## Citation

If you use Phlox in your research, please cite:

```bibtex
@article{phlox2024,
  title={Phlox: Automated Reaction Mechanism Discovery from Reactive MD},
  author={Chen, Xin},
  journal={Journal of Chemical Theory and Computation},
  year={2024}
}
```

## License

MIT License - see LICENSE file for details

## Author

**Chen Xin**
- Email: chenxin199261@gmail.com
- GitHub: [@chenxin](https://github.com/chenxin)

## Acknowledgments

Originally developed as ReaxANA/CPX-MechGen for complex reaction network analysis.
