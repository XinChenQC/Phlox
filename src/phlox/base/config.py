"""
Configuration management for Phlox.
"""

from dataclasses import dataclass, field
from typing import List, Optional, Tuple
from pathlib import Path


@dataclass
class Config:
    """
    Configuration for trajectory analysis.

    Attributes:
        # Computation settings
        n_cores: Number of CPU cores for parallel processing

        # Trajectory settings
        traj_type: Trajectory format ('xyz', 'lmp', 'extxyz')
        traj_files: List of trajectory file paths
        time_per_frame: Time per frame (value, unit)

        # PBC settings
        pbc_molecule: Whether molecules can be split by PBC boundary
        pbc_box: PBC box dimensions ((x0,y0,z0), (x1,y1,z1))

        # Analysis parameters
        bond_multiplier: Multiplier for bond distance criteria
        neighbor_cutoff: Cutoff for neighbor search (Angstroms)
        stability_time: Molecular lifetime threshold (value, unit)
        stability_lag: Frame lag for stability determination

        # LAMMPS-specific
        lammps_atom_types: Mapping of LAMMPS type to element

        # Catalyst settings
        enable_catalyst: Track catalyst atoms separately
        catalyst_selection_type: 1=by element, 2=by indices
        catalyst_indices: List of catalyst atom indices
        catalyst_element: Element symbol for catalyst selection

        # Output settings
        output_dir: Directory for output files
        log_level: Logging verbosity
    """
    # Computation
    n_cores: int = 1

    # Trajectory
    traj_type: str = "lmp"
    traj_files: List[str] = field(default_factory=list)
    time_per_frame: Tuple[float, str] = (1.0, "fs")

    # PBC
    pbc_molecule: bool = False
    pbc_box: Optional[Tuple[Tuple[float, float, float], Tuple[float, float, float]]] = None

    # Analysis
    bond_multiplier: float = 1.5
    neighbor_cutoff: float = 10.0
    neighbor_update_freq: int = 100  # Rebuild neighbor lists every N frames
    stability_time: Tuple[float, str] = (1.0, "ps")
    stability_lag: int = 10

    # LAMMPS
    lammps_atom_types: dict = field(default_factory=dict)

    # Catalyst
    enable_catalyst: bool = False
    catalyst_selection_type: int = 2
    catalyst_indices: List[int] = field(default_factory=list)
    catalyst_element: str = ""

    # Output
    output_dir: Path = field(default_factory=lambda: Path("./output"))
    log_level: str = "INFO"

    def __post_init__(self):
        """Validate configuration after initialization."""
        self.validate()

    def validate(self):
        """Validate configuration parameters."""
        if not self.traj_files:
            raise ValueError("At least one trajectory file must be specified")

        if self.traj_type not in ["xyz", "lmp", "extxyz"]:
            raise ValueError(f"Invalid trajectory type: {self.traj_type}")

        if self.traj_type in ["xyz", "extxyz"] and self.pbc_box is None:
            raise ValueError("PBC box must be specified for XYZ format")

        if self.n_cores < 1:
            raise ValueError("Number of cores must be at least 1")

        if self.bond_multiplier <= 0:
            raise ValueError("Bond multiplier must be positive")

        # Convert string paths to Path objects
        self.traj_files = [str(Path(f)) for f in self.traj_files]
        self.output_dir = Path(self.output_dir)

    def calculate_stability_lag(self) -> int:
        """
        Calculate frame lag for stability threshold.

        Returns:
            Number of frames corresponding to stability_time
        """
        stab_val, stab_unit = self.stability_time
        frame_val, frame_unit = self.time_per_frame

        # Convert to same units
        if stab_unit == frame_unit:
            self.stability_lag = int(stab_val / frame_val)
        elif stab_unit == "ps" and frame_unit == "fs":
            self.stability_lag = int((stab_val * 1000) / frame_val)
        elif stab_unit == "fs" and frame_unit == "ps":
            self.stability_lag = int(stab_val / (frame_val * 1000))
        else:
            raise ValueError(f"Cannot convert {stab_unit} to {frame_unit}")

        return self.stability_lag

    def to_dict(self) -> dict:
        """Convert configuration to dictionary."""
        return {
            'n_cores': self.n_cores,
            'traj_type': self.traj_type,
            'traj_files': self.traj_files,
            'time_per_frame': self.time_per_frame,
            'pbc_molecule': self.pbc_molecule,
            'pbc_box': self.pbc_box,
            'bond_multiplier': self.bond_multiplier,
            'neighbor_cutoff': self.neighbor_cutoff,
            'neighbor_update_freq': self.neighbor_update_freq,
            'stability_time': self.stability_time,
            'stability_lag': self.stability_lag,
            'lammps_atom_types': self.lammps_atom_types,
            'enable_catalyst': self.enable_catalyst,
            'catalyst_selection_type': self.catalyst_selection_type,
            'catalyst_indices': self.catalyst_indices,
            'catalyst_element': self.catalyst_element,
            'output_dir': str(self.output_dir),
            'log_level': self.log_level,
        }
