# Migration Guide: ReaxANA to Phlox

This guide helps you migrate from the old ReaxANA codebase to the new Phlox package.

## Overview of Changes

### 1. Package Structure

**Old (ReaxANA):**
```
ReaxANA/
├── Main.py
├── tool.py
├── global_var.py
├── AnalysFlow_parallel.py
├── ...
```

**New (Phlox):**
```
Phlox/
├── src/phlox/
│   ├── base/         # Core data structures and config
│   ├── frag/         # Fragmentation and distance calculations
│   ├── filter/       # Network cleaning
│   ├── io/           # File I/O
│   ├── network/      # Reaction network analysis
│   ├── parallel/     # Parallel processing
│   └── utils/        # Utilities
├── tests/            # Unit tests
├── docker/           # Docker configuration
└── examples/         # Example scripts
```

### 2. Configuration

**Old input format** (`inputfile`):
```
ncores 4
trajtype lmp
trajfiles traj.trj
...
```

**New format** (same, but validated):
```
# Phlox configuration
ncores 4
trajtype lmp
trajfiles traj.trj
timeperframe 0.1 fs
molelifttime 1 ps
...
```

Generate example config:
```bash
phlox config my_input.conf
```

### 3. Command Line Interface

**Old:**
```bash
python Main.py  # Uses hardcoded "inputfile"
```

**New:**
```bash
phlox analyze input.conf                    # Basic usage
phlox analyze input.conf --cores 8          # Override cores
phlox analyze input.conf --output results/  # Custom output
phlox analyze input.conf --verbose          # Verbose logging
```

### 4. Python API

**Old (not available):**
No programmatic API, only command-line script.

**New:**
```python
from phlox import Analyzer
from phlox.io import read_config

# Load configuration
config = read_config("input.conf")

# Create analyzer
analyzer = Analyzer(config)

# Run analysis
results = analyzer.run()

# Access results
print(f"Species: {len(results.species)}")
print(f"Reactions: {len(results.reactions)}")
```

### 5. Data Structures

**Old (global variables):**
```python
import global_var as gvar
gvar.atomList = [...]
gvar.blockList = [...]
gvar.GR = nx.MultiDiGraph()
```

**New (proper classes):**
```python
from phlox.base import Atom, Molecule, Config

atom = Atom(element="C", position=[1,2,3], index=0)
molecule = Molecule(mol_id=0, atom_indices=[0,1,2])
config = Config(n_cores=4, traj_type="lmp", ...)
```

### 6. Module Mapping

| Old File | New Module | Purpose |
|----------|------------|---------|
| `Main.py` | `phlox.cli` | Command-line interface |
| `global_var.py` | `phlox.base.constants` | Constants and defaults |
| `inputio.py` | `phlox.io.config_reader` | Configuration I/O |
| `tool.py` | `phlox.frag.distance`, `phlox.frag.fragmenter` | Distance and fragmentation |
| `union_find.py` | `phlox.frag.union_find` | UnionFind data structure |
| `AnalysFlow.py` | `phlox.base.analyzer` (serial) | Serial analysis |
| `AnalysFlow_parallel.py` | `phlox.parallel.worker` | Parallel analysis |
| `tool_reactCount.py` | `phlox.network.classifier` | Reaction classification |
| `tool_removeOscill.py` | `phlox.filter.oscillation` | Remove oscillations |
| `tool_contract.py` | `phlox.filter.contraction` | Contract nodes |
| `result_log.py` | `phlox.io.result_writer` | Output generation |

### 7. Installation

**Old:**
```bash
# No installation, just run scripts
python Main.py
```

**New:**
```bash
# Install as package
pip install phlox

# Or from source
git clone https://github.com/chenxin/phlox
cd phlox
pip install -e .

# With Docker
docker build -t phlox:latest -f docker/Dockerfile .
docker run -v $(pwd)/data:/data phlox analyze /data/input.conf
```

### 8. Testing

**Old:**
No test suite.

**New:**
```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=phlox --cov-report=html

# Run specific test
pytest tests/test_config.py
```

### 9. Output Files

**Output format remains the same:**
- `result.log` - Species and reactions
- `data.json` - Structured data
- `specRec/` - XYZ and SVG files
- `*.dot` - Graph files

**New features:**
- Configurable output directory
- Structured logging
- Progress tracking

### 10. Docker Support

**Old:**
Not available.

**New:**
```bash
# Build image
docker build -t phlox:latest -f docker/Dockerfile .

# Run analysis
docker run -v $(pwd)/data:/data -v $(pwd)/output:/output \
  phlox:latest analyze /data/input.conf --output /output

# Using docker-compose
cd docker
docker-compose up
```

## Migration Steps

### Step 1: Install Phlox

```bash
cd /home/xchen/Work/ReaxANA/Phlox
pip install -e .
```

### Step 2: Convert Input File

Your existing `inputfile` should work as-is. Optionally validate:

```bash
phlox config --help
```

### Step 3: Run Analysis

```bash
# Old way (still in ReaxANA directory)
cd /home/xchen/Work/ReaxANA/ReaxANA
python Main.py

# New way (from anywhere)
phlox analyze /path/to/inputfile --output /path/to/output
```

### Step 4: Update Scripts (if any)

If you have custom analysis scripts, update imports:

```python
# Old
import global_var as gvar
from tool import buildNeigh_AtomicBased

# New
from phlox.base import Config
from phlox.frag import Fragmenter
```

## Feature Parity Status

### ✅ Completed
- [x] Package structure
- [x] Configuration management
- [x] Base data structures (Atom, Molecule, Config)
- [x] CLI interface
- [x] Docker support
- [x] Unit tests
- [x] Documentation

### 🚧 In Progress (TODO items in code)
- [ ] Serial analysis implementation
- [ ] Parallel analysis implementation
- [ ] Network cleaning (oscillation removal, contraction)
- [ ] Reaction classification
- [ ] Output writing (result.log, data.json)
- [ ] SMILES generation
- [ ] Molecular hashing

### 📋 Planned
- [ ] Performance optimizations
- [ ] Additional trajectory formats
- [ ] Interactive visualization
- [ ] Web dashboard

## Backward Compatibility

The old ReaxANA code will continue to work in the `ReaxANA/` directory. The new Phlox package is in a separate `Phlox/` directory.

To run old code:
```bash
cd /home/xchen/Work/ReaxANA/ReaxANA
python Main.py
```

To run new code:
```bash
phlox analyze inputfile
```

## Getting Help

- Documentation: See `docs/` directory
- Issues: Report at GitHub issues page
- Email: chenxin199261@gmail.com

## Next Steps

1. **Complete Core Implementation**: Finish implementing the TODO items in analyzer.py
2. **Port Algorithms**: Migrate fragmentation, network cleaning, and reaction classification
3. **Add Examples**: Create example workflows in `examples/`
4. **Benchmark**: Compare performance with old code
5. **Publish**: Release to PyPI when ready
