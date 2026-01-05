"""
Configuration file reader and writer.
"""

from pathlib import Path
from typing import Union
from phlox.base.config import Config


def read_config(filepath: Union[str, Path]) -> Config:
    """
    Read configuration from input file.

    Args:
        filepath: Path to configuration file

    Returns:
        Config object

    Raises:
        FileNotFoundError: If config file doesn't exist
        ValueError: If config file has invalid format
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Configuration file not found: {filepath}")

    config_dict = {}

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()

            # Skip empty lines and comments
            if not line or line.startswith('#'):
                continue

            parts = line.split()
            if len(parts) < 2:
                continue

            key = parts[0].lower()
            values = parts[1:]

            # Parse different configuration options
            if key == 'ncores':
                config_dict['n_cores'] = int(values[0])

            elif key == 'trajtype':
                config_dict['traj_type'] = values[0].lower()

            elif key == 'trajfiles':
                config_dict['traj_files'] = values

            elif key == 'molelifttime':
                config_dict['stability_time'] = (float(values[0]), values[1].lower())

            elif key == 'timeperframe':
                config_dict['time_per_frame'] = (float(values[0]), values[1].lower())

            elif key == 'pbcmole':
                config_dict['pbc_molecule'] = values[0].lower() in ['true', 't', '1', 'yes']

            elif key == 'pbcbox':
                # Format: x0 x1 y0 y1 z0 z1
                config_dict['pbc_box'] = (
                    (float(values[0]), float(values[2]), float(values[4])),
                    (float(values[1]), float(values[3]), float(values[5]))
                )

            elif key == 'lmpelement':
                # Map atom types to elements: 1->C, 2->H, etc.
                config_dict['lammps_atom_types'] = {
                    i+1: elem for i, elem in enumerate(values)
                }

            elif key == 'catalysis':
                config_dict['enable_catalyst'] = values[0].lower() in ['true', 't', '1', 'yes']

            elif key == 'selectype':
                config_dict['catalyst_selection_type'] = int(values[0])

            elif key == 'cataselection':
                selection = values[0]
                if '-' in selection:
                    # Range: "1-10"
                    start, end = selection.split('-')
                    config_dict['catalyst_indices'] = list(range(int(start)-1, int(end)))
                elif ',' in selection:
                    # List: "1,5,9,12"
                    config_dict['catalyst_indices'] = [int(x)-1 for x in selection.split(',')]
                else:
                    # Single element name
                    config_dict['catalyst_element'] = selection

            elif key == 'bondmultiplier':
                config_dict['bond_multiplier'] = float(values[0])

            elif key == 'outputdir':
                config_dict['output_dir'] = Path(values[0])

            elif key == 'loglevel':
                config_dict['log_level'] = values[0].upper()

    # Create Config object
    config = Config(**config_dict)
    config.calculate_stability_lag()

    return config


def write_config(config: Config, filepath: Union[str, Path]):
    """
    Write configuration to file.

    Args:
        config: Config object to write
        filepath: Output file path
    """
    filepath = Path(filepath)

    with open(filepath, 'w') as f:
        f.write("# Phlox Configuration File\n")
        f.write("# Generated automatically\n\n")

        f.write("# Computation Settings\n")
        f.write(f"ncores          {config.n_cores}\n\n")

        f.write("# Trajectory Settings\n")
        f.write(f"trajtype        {config.traj_type}\n")
        f.write(f"trajfiles       {' '.join(config.traj_files)}\n")
        f.write(f"timeperframe    {config.time_per_frame[0]} {config.time_per_frame[1]}\n\n")

        f.write("# PBC Settings\n")
        f.write(f"pbcmole         {config.pbc_molecule}\n")
        if config.pbc_box:
            f.write(f"pbcbox          {config.pbc_box[0][0]} {config.pbc_box[1][0]}   "
                   f"{config.pbc_box[0][1]} {config.pbc_box[1][1]}   "
                   f"{config.pbc_box[0][2]} {config.pbc_box[1][2]}\n\n")

        f.write("# Analysis Parameters\n")
        f.write(f"bondmultiplier  {config.bond_multiplier}\n")
        f.write(f"molelifttime    {config.stability_time[0]} {config.stability_time[1]}\n\n")

        if config.lammps_atom_types:
            f.write("# LAMMPS Settings\n")
            elements = [config.lammps_atom_types[i+1] for i in range(len(config.lammps_atom_types))]
            f.write(f"lmpelement      {' '.join(elements)}\n\n")

        if config.enable_catalyst:
            f.write("# Catalyst Settings\n")
            f.write(f"catalysis       True\n")
            f.write(f"selectype       {config.catalyst_selection_type}\n")
            if config.catalyst_indices:
                # Write as range or list
                indices_str = f"{config.catalyst_indices[0]+1}-{config.catalyst_indices[-1]+1}"
                f.write(f"cataselection   {indices_str}\n")
            elif config.catalyst_element:
                f.write(f"cataselection   {config.catalyst_element}\n")
            f.write("\n")

        f.write("# Output Settings\n")
        f.write(f"outputdir       {config.output_dir}\n")
        f.write(f"loglevel        {config.log_level}\n")
