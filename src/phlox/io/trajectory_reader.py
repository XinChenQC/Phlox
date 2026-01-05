"""
Trajectory file readers for different formats.
"""

from abc import ABC, abstractmethod
from typing import List, Tuple, Optional
import numpy as np
from pathlib import Path

from phlox.base.atom import Atom
from phlox.base.config import Config


class TrajectoryReader(ABC):
    """Base class for trajectory readers."""

    def __init__(self, filepath: str, config: Config):
        """
        Initialize reader.

        Args:
            filepath: Path to trajectory file
            config: Configuration object
        """
        self.filepath = Path(filepath)
        self.config = config
        self.current_frame = 0

        if not self.filepath.exists():
            raise FileNotFoundError(f"Trajectory file not found: {filepath}")

    @abstractmethod
    def read_frame(self) -> Optional[Tuple[List[Atom], Optional[tuple]]]:
        """
        Read next frame from trajectory.

        Returns:
            Tuple of (atom_list, pbc_box) or None if end of file
            pbc_box format: ((x0,y0,z0), (x1,y1,z1))
        """
        pass

    @abstractmethod
    def count_frames(self) -> int:
        """Count total number of frames in trajectory."""
        pass

    def reset(self):
        """Reset reader to beginning of file."""
        self.current_frame = 0

    @abstractmethod
    def skip_frames(self, n_frames: int):
        """
        Skip n frames from current position.

        Args:
            n_frames: Number of frames to skip
        """
        pass


class XYZReader(TrajectoryReader):
    """Reader for XYZ format trajectories."""

    def __init__(self, filepath: str, config: Config):
        super().__init__(filepath, config)
        self.file_handle = None

    def __enter__(self):
        self.file_handle = open(self.filepath, 'r')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.file_handle:
            self.file_handle.close()

    def read_frame(self) -> Optional[Tuple[List[Atom], Optional[tuple]]]:
        """
        Read one frame from XYZ file.

        XYZ format:
            <number of atoms>
            comment line
            <element> <x> <y> <z>
            ...

        Returns:
            Tuple of (atoms, pbc_box) or None if EOF
        """
        if self.file_handle is None:
            raise RuntimeError("File not opened. Use 'with' context manager.")

        # Read number of atoms
        line = self.file_handle.readline()
        if not line:
            return None

        try:
            n_atoms = int(line.strip())
        except ValueError:
            return None

        # Read comment line
        self.file_handle.readline()

        # Read atoms
        atoms = []
        for i in range(n_atoms):
            line = self.file_handle.readline()
            if not line:
                raise ValueError(f"Unexpected EOF while reading frame {self.current_frame}")

            parts = line.split()
            element = parts[0]
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])

            atom = Atom(
                element=element,
                position=np.array([x, y, z], dtype=np.float32),
                index=i
            )
            atoms.append(atom)

        self.current_frame += 1
        return atoms, self.config.pbc_box

    def count_frames(self) -> int:
        """Count frames by reading through file."""
        with open(self.filepath, 'r') as f:
            frame_count = 0
            while True:
                line = f.readline()
                if not line:
                    break
                try:
                    n_atoms = int(line.strip())
                    # Skip comment line and atom lines
                    for _ in range(n_atoms + 1):
                        f.readline()
                    frame_count += 1
                except (ValueError, AttributeError):
                    break
        return frame_count

    def skip_frames(self, n_frames: int):
        """
        Skip n frames from current position.

        Args:
            n_frames: Number of frames to skip
        """
        if self.file_handle is None:
            raise RuntimeError("File not opened. Use 'with' context manager.")

        for _ in range(n_frames):
            # Read number of atoms
            line = self.file_handle.readline()
            if not line:
                return

            try:
                n_atoms = int(line.strip())
            except ValueError:
                return

            # Skip comment line and atom lines
            for _ in range(n_atoms + 1):
                self.file_handle.readline()

            self.current_frame += 1


class LAMMPSReader(TrajectoryReader):
    """Reader for LAMMPS dump format trajectories."""

    def __init__(self, filepath: str, config: Config):
        super().__init__(filepath, config)
        self.file_handle = None
        self.n_atoms = 0

    def __enter__(self):
        self.file_handle = open(self.filepath, 'r')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.file_handle:
            self.file_handle.close()

    def read_frame(self) -> Optional[Tuple[List[Atom], Optional[tuple]]]:
        """
        Read one frame from LAMMPS dump file.

        LAMMPS dump format:
            ITEM: TIMESTEP
            <timestep>
            ITEM: NUMBER OF ATOMS
            <n_atoms>
            ITEM: BOX BOUNDS
            <x0> <x1>
            <y0> <y1>
            <z0> <z1>
            ITEM: ATOMS id type x y z
            <id> <type> <x> <y> <z>
            ...

        Returns:
            Tuple of (atoms, pbc_box) or None if EOF
        """
        if self.file_handle is None:
            raise RuntimeError("File not opened. Use 'with' context manager.")

        atoms = []
        pbc_box = None
        read_state = None

        while True:
            line = self.file_handle.readline()
            if not line:
                return None if not atoms else (atoms, pbc_box)

            line = line.strip()

            if "ITEM: TIMESTEP" in line:
                read_state = "timestep"
                continue

            elif "ITEM: NUMBER OF ATOMS" in line:
                read_state = "natoms"
                continue

            elif "ITEM: BOX BOUNDS" in line:
                read_state = "box"
                box_data = []
                continue

            elif "ITEM: ATOMS" in line:
                read_state = "atoms"
                # Parse column headers
                parts = line.split()
                parts.remove("ITEM:")
                parts.remove("ATOMS")
                idx_id = parts.index('id')
                idx_type = parts.index('type')
                idx_x = parts.index('x')
                idx_y = parts.index('y')
                idx_z = parts.index('z')
                atom_count = 0
                continue

            # Process based on state
            if read_state == "timestep":
                timestep = int(line)
                read_state = None

            elif read_state == "natoms":
                self.n_atoms = int(line)
                read_state = None

            elif read_state == "box":
                box_data.append([float(x) for x in line.split()])
                if len(box_data) == 3:
                    pbc_box = (
                        (box_data[0][0], box_data[1][0], box_data[2][0]),
                        (box_data[0][1], box_data[1][1], box_data[2][1])
                    )
                    read_state = None

            elif read_state == "atoms":
                parts = line.split()
                atom_id = int(parts[idx_id]) - 1  # Convert to 0-indexed
                atom_type = int(parts[idx_type])
                x = float(parts[idx_x])
                y = float(parts[idx_y])
                z = float(parts[idx_z])

                element = self.config.lammps_atom_types.get(atom_type, f"X{atom_type}")

                atom = Atom(
                    element=element,
                    position=np.array([x, y, z], dtype=np.float32),
                    index=atom_id
                )
                atoms.append(atom)
                atom_count += 1

                if atom_count == self.n_atoms:
                    # Sort atoms by index
                    atoms.sort(key=lambda a: a.index)
                    self.current_frame += 1
                    return atoms, pbc_box

        return None

    def count_frames(self) -> int:
        """Count frames in LAMMPS dump file."""
        with open(self.filepath, 'r') as f:
            frame_count = 0
            for line in f:
                if "ITEM: TIMESTEP" in line:
                    frame_count += 1
        return frame_count

    def skip_frames(self, n_frames: int):
        """
        Skip n frames from current position.

        Args:
            n_frames: Number of frames to skip
        """
        if self.file_handle is None:
            raise RuntimeError("File not opened. Use 'with' context manager.")

        for _ in range(n_frames):
            result = self.read_frame()
            if result is None:
                return


def create_reader(filepath: str, config: Config) -> TrajectoryReader:
    """
    Factory function to create appropriate trajectory reader.

    Args:
        filepath: Path to trajectory file
        config: Configuration object

    Returns:
        TrajectoryReader instance

    Raises:
        ValueError: If trajectory type is not supported
    """
    if config.traj_type in ["xyz", "extxyz"]:
        return XYZReader(filepath, config)
    elif config.traj_type == "lmp":
        return LAMMPSReader(filepath, config)
    else:
        raise ValueError(f"Unsupported trajectory type: {config.traj_type}")
