#!/usr/bin/env python
"""Test Phlox with real ReaxANA trajectory data."""

from phlox import Analyzer, Config
import logging

# Enable INFO logging to see fragmentation summaries
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Create config matching the ReaxANA inputfile
config = Config(
    n_cores=2,  # Use single core for testing
    traj_type="xyz",
#    traj_type="lmp",
#    traj_files=["/home/xchen/Work/ReaxANA/ReaxANA/reax2.trj"],
    traj_files=["/home/xchen/Work/ReaxANA/ReaxANA/trajectory_DFT2SP.xyz"],
    pbc_box=((0, 0, 0), (18, 18, 25)),  # From Lattice in XYZ file
    output_dir="./phlox_test_output",
#    lammps_atom_types={1: "C", 2: "H", 3:"N", 4:"O"},
    stability_lag=40,  # From molelifttime 20 ps / 0.5 ps per frame
)

print("=" * 80)
print("Testing Phlox with Real ReaxANA Trajectory")
print("=" * 80)
print(f"Input file: {config.traj_files[0]}")
print(f"PBC box: {config.pbc_box}")
print(f"Stability lag: {config.stability_lag} frames")
print()

# Run analysis
analyzer = Analyzer(config)
results = analyzer.run()

print()
print("=" * 80)
print("Analysis Complete!")
print("=" * 80)
print(f"Total frames processed: {results.total_frames}")
print(f"Network nodes: {results.graph.number_of_nodes()}")
print(f"Network edges: {results.graph.number_of_edges()}")

if results.reactions:
    print(f"Reaction types found: {len(results.reactions)}")
    total_events = sum(len(v) for v in results.reactions.values())
    print(f"Total reaction events: {total_events}")

print()
print("=" * 80)
print("[species]")
print("=" * 80)
print("#" + "-" * 124)
print(f"#{'SMILES':<50} {'filename':<25} {'Formula':<15} {'Abundance':>10}")
print("#" + "-" * 124)

# Build time series by speciesID (like ReaxANA's SpeciesCount)
# results.time_series is List[Counter({speciesID: count})]
time_series_by_species = {}

# Collect all species IDs from struct_dict
all_species_ids = set(results.struct_dict.keys())

# Initialize time series for each species
for species_id in all_species_ids:
    time_series_by_species[species_id] = []

# Build time series for each species across all frames
if results.time_series:
    for frame_counter in results.time_series:
        # For each species, get count from Counter (0 if not present)
        for species_id in all_species_ids:
            count = frame_counter.get(species_id, 0)
            time_series_by_species[species_id].append(count)

# Group by SMILES (like ReaxANA's smileRev)
smiles_time_series = {}
for species_id, time_series in time_series_by_species.items():
    species_data = results.struct_dict[species_id]
    smiles = species_data.get('smiles', '')
    if smiles:
        if smiles not in smiles_time_series:
            smiles_time_series[smiles] = time_series[:]
        else:
            # Combine time series for same SMILES
            smiles_time_series[smiles] = [
                smiles_time_series[smiles][i] + time_series[i]
                for i in range(len(time_series))
            ]

# Calculate abundance as MAXIMUM count across all timesteps (like ReaxANA line 80)
smiles_info = {}
max_smiles_length = 0
for species_id, species_data in results.struct_dict.items():
    smiles = species_data.get('smiles', '')
    formula = species_data.get('formula', 'Unknown')

    if smiles and smiles not in smiles_info:
        if len(smiles) > max_smiles_length:
            max_smiles_length = len(smiles)

        # Abundance = max count in any single timestep
        time_series = smiles_time_series.get(smiles, [])
        abundance = max(time_series) if time_series else 0
        # Remove "H" prefix from speciesID for filename (if present)
        species_id_short = species_id[1:8] if species_id.startswith('H') else species_id[:7]
        filename = f"S_{formula}_H{species_id_short}"
        smiles_info[smiles] = {
            'abundance': abundance,
            'filename': filename,
            'formula': formula,
            'time_series': smiles_time_series.get(smiles, [])
        }

# Sort by abundance (descending)
sorted_smiles = sorted(smiles_info.items(), key=lambda x: x[1]['abundance'], reverse=True)

# Print species table
for smiles, info in sorted_smiles:
    print(f"{smiles:<50} {info['filename']:<25} {info['formula']:<15} {info['abundance']:>10}")

print()
print(f"Total unique species: {len(sorted_smiles)}")
print("[species_end]")
print()

# Print species evolution (like ReaxANA's [speciesrev])
print("[speciesrev]")
total_frames = len(smiles_time_series[list(smiles_time_series.keys())[0]]) if smiles_time_series else 0
print(f"{len(smiles_info):5}{total_frames:10}")
for smiles, info in sorted_smiles:
    print(f"{smiles:<50}")
    time_series = info['time_series']
    for i, count in enumerate(time_series, 1):
        if i % 17 == 0:
            print(f"{count:5}")
        else:
            print(f"{count:5}", end='')
    print()
print("[speciesrev_end]")
print()

# Print formula evolution (like ReaxANA's [formularev])
print("[formularev]")
# Group by formula
formula_time_series = {}
for smiles, info in smiles_info.items():
    formula = info['formula']
    if formula not in formula_time_series:
        formula_time_series[formula] = [0] * total_frames
    # Add this species' time series to the formula
    for i, count in enumerate(info['time_series']):
        formula_time_series[formula][i] += count

for formula, time_series in sorted(formula_time_series.items(), key=lambda x: max(x[1]) if x[1] else 0, reverse=True):
    print(f"{formula:<50}")
    for i, count in enumerate(time_series, 1):
        if i % 17 == 0:
            print(f"{count:5}")
        else:
            print(f"{count:5}", end='')
    print()
print("[formularev_end]")
print()
'''
# Print reaction classification if available
if results.reactions:
    print()
    print("=" * 80)
    print("[reactions]")
    print("=" * 80)

    reaction_count = 0
    event_number = 0
    for reaction_type, events in results.reactions.items():
        for event in events:
            event_number += 1

            # Get SMILES for reactants and products
            reactant_smiles = event.get('reactant_smiles', [''])
            product_smiles = event.get('product_smiles', [''])
            reactant_atoms = event.get('reactant_atoms', [[]])
            product_atoms = event.get('product_atoms', [[]])
            timestep = event.get('timestep', 0)

            # Format reactants and products
            react_str = '.'.join([s if s else '?' for s in reactant_smiles])
            prod_str = '.'.join([s if s else '?' for s in product_smiles])

            # Format atom indices
            react_atom_str = '   -->   '.join([str(atoms) for atoms in reactant_atoms])
            prod_atom_str = '   -->   '.join([str(atoms) for atoms in product_atoms])

            # Print reaction in ReaxANA format
            print(f" #{event_number} {react_str}  -->  {prod_str}  # of events: 1")
            print(f"   {timestep}:    {react_atom_str}   -->   {prod_atom_str}")

        reaction_count += 1

    print(f"\n[Reactionrev_end]")
    print(f"\nTotal reaction types: {reaction_count}")
    print(f"Total reaction events: {event_number}")
'''
print()
print(f"Output files in: {config.output_dir}/")
print("  - reaction_network.dot")
print("  - stats.txt")
print()
