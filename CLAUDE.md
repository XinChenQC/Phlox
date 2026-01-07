# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Phlox (formerly ReaxANA/CPX-MechGen) is a Python package for automated reaction mechanism discovery from reactive molecular dynamics trajectories. It analyzes MD simulations to detect molecular fragmentation, build reaction networks, and classify reaction pathways.

**Current Development Status:**
- Fragment Analysis: Complete
- Network Filters: In progress (oscillation removal and node contraction currently disabled)
- Reaction Classification: Partial

## Commands

### Development Setup

```bash
# Install with development dependencies (recommended: use uv)
uv venv
uv pip install -e ".[dev]"

# Or with pip
pip install -e ".[dev]"
```

### Running Tests

```bash
# Run all tests with coverage
make test
# Or directly:
pytest tests/ -v --cov=phlox --cov-report=html --cov-report=term

# Run a single test file
pytest tests/test_atom.py -v

# Run a specific test
pytest tests/test_config.py::TestConfig::test_validate -v
```

### Code Quality

```bash
# Format code (black + isort)
make format

# Run linters
make lint

# Individual tools
black src/phlox tests
isort src/phlox tests
flake8 src/phlox tests
mypy src/phlox
```

### Running Analysis

```bash
# CLI usage (once implemented)
phlox analyze input.conf

# Python API
python test_real_data.py  # Example test script
```

### Building

```bash
# Build distribution packages
make build  # or: python -m build

# Build Docker image
make docker
```

## Architecture Overview

### Package Structure

```
src/phlox/
├── base/           # Core data structures and orchestration
│   ├── analyzer.py       # Main analysis orchestration (serial/parallel)
│   ├── config.py         # Configuration management
│   ├── atom.py           # Atom representation
│   ├── molecule.py       # Molecule representation
│   ├── results.py        # Results container
│   └── frame_processor.py # Per-frame processing logic
├── frag/           # Molecular fragmentation
│   ├── fragmenter.py     # Main fragmentation logic
│   ├── distance.py       # PBC-aware distance calculations
│   └── union_find.py     # Union-Find for connectivity
├── network/        # Reaction network construction
│   ├── builder.py        # Network builder
│   └── classifier.py     # Reaction classification
├── filter/         # Network cleaning (CURRENTLY DISABLED)
│   ├── oscillation.py    # Remove spurious oscillations
│   └── contraction.py    # Contract unstable nodes
├── io/             # Input/output operations
│   ├── config_reader.py  # Configuration file parsing
│   └── trajectory_reader.py # XYZ/LAMMPS trajectory readers
├── utils/          # Chemistry utilities
│   ├── smiles.py         # SMILES generation
│   ├── xyz2mol.py        # XYZ to molecular graph conversion
│   └── batch_smiles.py   # Batch SMILES computation
└── cli.py          # Command-line interface
```

### Execution Flow

1. **Config Loading** (`phlox.io.config_reader`) → Parse input configuration
2. **Analyzer Init** (`phlox.base.analyzer`) → Validate and setup analysis
3. **Analysis Execution**:
   - **Serial** (`n_cores=1`): Single process, frame-by-frame
   - **Parallel** (`n_cores>1`): Split trajectory across workers, merge results
4. **Per-Frame Processing** (`phlox.base.frame_processor`):
   - Read atoms from trajectory
   - Identify molecular fragments (`phlox.frag.fragmenter`)
   - Build reaction network incrementally (`phlox.network.builder`)
   - Track species over time
5. **Post-Processing**:
   - SMILES generation for all species (`phlox.utils.batch_smiles`)
   - Network cleaning - DISABLED (`phlox.filter`)
   - Reaction classification (`phlox.network.classifier`)
6. **Output** (`phlox.base.results`):
   - `reaction_network.dot` - GraphViz network file
   - `stats.txt` - Summary statistics
   - `species/*.xyz` - Individual species structures
   - `species/species_list.txt` - Species list with SMILES

### Key Data Structures

- **Config**: Immutable configuration object, validated on creation
- **Atoms** → **Molecules**: Raw atoms grouped into molecules by fragmenter
- **struct_dict**: Maps species ID (molecular hash) → structure info
- **molecule_info**: Maps molecule ID → atom indices
- **time_series**: List of Counters tracking species abundance per frame
- **graph**: NetworkX MultiDiGraph (nodes=species, edges=reactions)
- **Results**: Container aggregating all output data

### Parallel Processing

When `n_cores > 1`:
- Trajectory is split into chunks (first worker gets 120% of base frames for warm-up, last gets 80%)
- Each worker processes its chunk independently with local data structures
- Main process merges all worker results
- Cross-boundary reactions are detected by comparing last frame of worker N with first frame of worker N+1

## Important Implementation Details

### Molecular Fragmentation
- Uses Union-Find algorithm for molecular connectivity
- Distance calculations respect periodic boundary conditions (PBC)
- Neighbor lists updated every 20 frames (configurable via `neighbor_update_freq`)
- Initial frame uses full distance matrix; subsequent frames use incremental updates

### Network Filtering
**CRITICAL**: Network filtering is currently DISABLED to preserve all reaction events during development. The oscillation removal and node contraction modules exist but are not called from the main analysis flow. Do not re-enable without careful review.

### SMILES Generation
- Performed in batch after trajectory processing to avoid redundant calculations
- Uses RDKit via `xyz2mol` conversion
- SMILES strings stored in `struct_dict` for each unique species

### Configuration Parameters
- `stability_lag`: Number of frames a molecule must exist to be considered stable
- `bond_multiplier`: Multiplier on covalent radii for bond detection (default 1.5)
- `neighbor_cutoff`: Cutoff distance for neighbor search in Angstroms
- `pbc_molecule`: Whether molecules can be split across PBC boundaries

### Trajectory Format Support
- **XYZ**: Requires explicit `pbc_box` parameter
- **LAMMPS**: Reads box from dump file; requires `lammps_atom_types` mapping
- **ExtXYZ**: Extended XYZ with lattice information in comment line

## Testing Patterns

- Tests use pytest with fixtures defined in `tests/conftest.py`
- Coverage reports generated in `htmlcov/`
- Example test data paths point to `/home/xchen/Work/ReaxANA/ReaxANA/`
- Test real data with `test_real_data.py` (adjust paths as needed)

## Code Style

- Line length: 100 characters (black/isort configured)
- Python 3.8+ required
- Type hints optional (mypy configured but `disallow_untyped_defs = false`)
- Imports organized by isort with black profile

## Dependencies

Core:
- numpy, scipy: Numerical operations
- networkx: Reaction network graphs
- rdkit: Chemistry utilities and SMILES generation
- pydot: GraphViz output

Dev:
- pytest, pytest-cov: Testing
- black, isort, flake8, mypy: Code quality
