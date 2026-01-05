"""
Example: Basic usage of Phlox for trajectory analysis.

This example shows how to:
1. Create a configuration
2. Run analysis
3. Access results
4. Generate custom outputs
"""

from pathlib import Path
from phlox import Analyzer, Config
from phlox.io import read_config, write_config


def example_1_create_config():
    """Create configuration programmatically."""
    print("Example 1: Creating configuration programmatically")
    print("=" * 60)

    config = Config(
        n_cores=4,
        traj_type="lmp",
        traj_files=["data/trajectory.trj"],
        time_per_frame=(0.1, "fs"),
        pbc_molecule=True,
        stability_time=(1.0, "ps"),
        lammps_atom_types={1: "C", 2: "H", 3: "O", 4: "Ni"},
        output_dir=Path("output/example1"),
    )

    # Save configuration
    write_config(config, "output/example1_config.conf")
    print(f"Configuration saved to output/example1_config.conf")
    print(f"Stability lag: {config.stability_lag} frames\n")


def example_2_load_and_analyze():
    """Load configuration from file and run analysis."""
    print("Example 2: Load configuration and run analysis")
    print("=" * 60)

    # Load configuration
    config = read_config("input.conf")
    print(f"Loaded config: {config.n_cores} cores, {len(config.traj_files)} files")

    # Create analyzer
    analyzer = Analyzer(config)

    # Run analysis (this will take time for real trajectories)
    # results = analyzer.run()

    # Print summary
    # print(results.summary())


def example_3_custom_analysis():
    """Perform custom analysis on results."""
    print("Example 3: Custom analysis")
    print("=" * 60)

    from phlox.io import read_config

    config = read_config("input.conf")

    # Custom workflow
    print("Custom analysis workflow:")
    print("1. Read trajectory")
    print("2. Identify fragments")
    print("3. Build reaction network")
    print("4. Filter and classify")
    print("5. Export results")


def example_4_catalyst_tracking():
    """Example with catalyst tracking enabled."""
    print("Example 4: Catalyst tracking")
    print("=" * 60)

    config = Config(
        n_cores=4,
        traj_type="lmp",
        traj_files=["data/catalytic_reaction.trj"],
        time_per_frame=(0.1, "fs"),
        pbc_molecule=True,
        stability_time=(1.0, "ps"),
        lammps_atom_types={1: "C", 2: "H", 3: "O", 4: "Ni"},
        enable_catalyst=True,
        catalyst_selection_type=2,
        catalyst_indices=list(range(0, 10)),  # First 10 atoms are catalyst
        output_dir=Path("output/catalyst"),
    )

    print(f"Catalyst tracking enabled")
    print(f"Catalyst atoms: {config.catalyst_indices}")


def example_5_parallel_processing():
    """Example with parallel processing."""
    print("Example 5: Parallel processing")
    print("=" * 60)

    config = Config(
        n_cores=8,  # Use 8 cores
        traj_type="xyz",
        traj_files=["data/traj1.xyz", "data/traj2.xyz"],
        time_per_frame=(1.0, "fs"),
        pbc_molecule=True,
        pbc_box=((0, 0, 0), (50, 50, 50)),
        stability_time=(10.0, "fs"),
        output_dir=Path("output/parallel"),
    )

    print(f"Parallel processing with {config.n_cores} cores")
    print(f"Processing {len(config.traj_files)} trajectory files")


if __name__ == "__main__":
    # Create output directory
    Path("output").mkdir(exist_ok=True)

    # Run examples
    try:
        example_1_create_config()
    except Exception as e:
        print(f"Example 1 failed: {e}\n")

    try:
        example_2_load_and_analyze()
    except Exception as e:
        print(f"Example 2 skipped (no input file): {e}\n")

    try:
        example_3_custom_analysis()
    except Exception as e:
        print(f"Example 3 failed: {e}\n")

    try:
        example_4_catalyst_tracking()
    except Exception as e:
        print(f"Example 4 failed: {e}\n")

    try:
        example_5_parallel_processing()
    except Exception as e:
        print(f"Example 5 failed: {e}\n")

    print("\n" + "=" * 60)
    print("Examples complete!")
    print("Note: Full analysis requires actual trajectory files")
