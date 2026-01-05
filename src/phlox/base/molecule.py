"""
Molecule data structure and related functionality.
"""

from dataclasses import dataclass, field
from typing import List, Optional
import numpy as np
import hashlib


@dataclass
class Molecule:
    """
    Represents a molecular fragment.

    Attributes:
        mol_id: Unique molecule ID
        atom_indices: List of atom indices belonging to this molecule
        hash: Molecular structure hash (Weisfeiler-Lehman)
        smiles: SMILES representation
        formula: Molecular formula (e.g., 'C2H4O')
        center: Center of mass coordinates
        radius: Molecular radius (max distance from center)
        fragment_hash: Unique hash for this specific fragment instance
        elements: List of element symbols (stored for SMILES generation)
        coordinates: Atomic coordinates (stored for SMILES generation)
        connectivity: Connectivity matrix (stored for SMILES generation)
    """
    mol_id: int
    atom_indices: List[int] = field(default_factory=list)
    hash: str = ""
    smiles: str = ""
    formula: str = ""
    center: np.ndarray = field(default_factory=lambda: np.zeros(3))
    radius: float = 0.0
    fragment_hash: str = ""
    elements: List[str] = field(default_factory=list)
    coordinates: List[List[float]] = field(default_factory=list)
    connectivity: Optional[np.ndarray] = None

    def calculate_fragment_hash(self) -> str:
        """
        Generate unique hash for this fragment instance.

        Combines atom labels and molecular hash to create unique identifier.
        ReaxANA approach: SHA1(atom_indices_string + speciesID)

        Returns:
            SHA1 hash string (first 20 characters)
        """
        # Sort atom indices first (ReaxANA: rec[1].sort() before hashing)
        sorted_indices = sorted(self.atom_indices)
        atom_label_str = "".join([str(i) for i in sorted_indices])
        # Concatenate with species hash (already has "H" prefix)
        hash_str = atom_label_str + self.hash
        self.fragment_hash = hashlib.sha1(hash_str.encode('utf-8')).hexdigest()[:20]
        return self.fragment_hash

    def calculate_center(self, positions: np.ndarray) -> np.ndarray:
        """
        Calculate center of mass.

        Args:
            positions: Array of atomic positions (N x 3)

        Returns:
            Center of mass coordinates
        """
        if len(self.atom_indices) == 0:
            self.center = np.zeros(3)
            return self.center

        mol_positions = positions[self.atom_indices]
        self.center = np.mean(mol_positions, axis=0)
        return self.center

    def calculate_radius(self, positions: np.ndarray) -> float:
        """
        Calculate molecular radius (max distance from center).

        Args:
            positions: Array of atomic positions (N x 3)

        Returns:
            Radius in Angstroms
        """
        if len(self.atom_indices) == 0:
            self.radius = 0.0
            return self.radius

        mol_positions = positions[self.atom_indices]
        distances = np.linalg.norm(mol_positions - self.center, axis=1)
        self.radius = np.max(distances) if len(distances) > 0 else 0.0
        return self.radius

    def __repr__(self) -> str:
        return f"Molecule(id={self.mol_id}, formula={self.formula}, atoms={len(self.atom_indices)})"
