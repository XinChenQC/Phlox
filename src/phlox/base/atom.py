"""
Atom data structure and related functionality.
"""

from dataclasses import dataclass, field
from typing import List
import numpy as np


@dataclass
class Atom:
    """
    Represents a single atom in the simulation.

    Attributes:
        element: Element symbol (e.g., 'C', 'H', 'O')
        position: 3D coordinates [x, y, z] in Angstroms
        index: Atom index in the global atom list
        neighbors: List of neighboring atom indices
        molecule_id: ID of the molecule this atom belongs to (set during analysis)
    """
    element: str
    position: np.ndarray
    index: int
    neighbors: List[int] = field(default_factory=list)
    molecule_id: int = -1

    def __post_init__(self):
        """Ensure position is a numpy array."""
        if not isinstance(self.position, np.ndarray):
            self.position = np.array(self.position, dtype=np.float32)

    def distance_to(self, other: 'Atom', pbc_box: tuple = None) -> float:
        """
        Calculate distance to another atom.

        Args:
            other: Another Atom instance
            pbc_box: Optional PBC box dimensions ((x0,y0,z0), (x1,y1,z1))

        Returns:
            Distance in Angstroms
        """
        delta = other.position - self.position

        if pbc_box is not None:
            box_size = np.array([
                pbc_box[1][0] - pbc_box[0][0],
                pbc_box[1][1] - pbc_box[0][1],
                pbc_box[1][2] - pbc_box[0][2]
            ])
            # Apply minimum image convention
            delta = delta - box_size * np.round(delta / box_size)

        return np.linalg.norm(delta)

    def add_neighbor(self, atom_index: int):
        """Add a neighboring atom index."""
        if atom_index not in self.neighbors:
            self.neighbors.append(atom_index)

    def clear_neighbors(self):
        """Clear the neighbor list."""
        self.neighbors.clear()

    def __repr__(self) -> str:
        return f"Atom({self.element}, idx={self.index}, pos={self.position}, neighbors={self.neighbors})"
