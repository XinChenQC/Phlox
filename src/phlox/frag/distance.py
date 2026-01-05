"""
Distance calculation utilities with PBC support.
"""

import numpy as np
from scipy.spatial.distance import pdist, squareform, cdist
from typing import Optional, Tuple


class DistanceCalculator:
    """
    Calculate interatomic distances with periodic boundary conditions.
    """

    @staticmethod
    def pairwise_distances(
        positions: np.ndarray,
        pbc_box: Optional[Tuple[Tuple[float, float, float], Tuple[float, float, float]]] = None
    ) -> np.ndarray:
        """
        Calculate pairwise distance matrix.

        Args:
            positions: N x 3 array of atomic positions
            pbc_box: Optional PBC box ((x0,y0,z0), (x1,y1,z1))

        Returns:
            N x N distance matrix
        """
        if pbc_box is None:
            # No PBC - use simple Euclidean distance
            distances = squareform(pdist(positions))
            return distances.astype(np.float32)

        # With PBC - apply minimum image convention
        box_size = np.array([
            pbc_box[1][0] - pbc_box[0][0],
            pbc_box[1][1] - pbc_box[0][1],
            pbc_box[1][2] - pbc_box[0][2]
        ])

        # Calculate displacements
        n = len(positions)
        displacements = np.zeros((3, n * (n - 1) // 2))

        for i, pos in enumerate(positions.T):
            pdist(pos[:, None], 'cityblock', out=displacements[i])

        # Apply minimum image convention
        for i in range(3):
            mask = displacements[i] > box_size[i] / 2
            displacements[i][mask] -= box_size[i]

        # Calculate distances
        distances = np.linalg.norm(displacements, axis=0)
        return squareform(distances).astype(np.float32)

    @staticmethod
    def neighbor_distances(
        atom_position: np.ndarray,
        neighbor_positions: np.ndarray,
        pbc_box: Optional[Tuple[Tuple[float, float, float], Tuple[float, float, float]]] = None
    ) -> np.ndarray:
        """
        Calculate distances from one atom to a list of neighbors.

        Args:
            atom_position: 1D array [x, y, z]
            neighbor_positions: N x 3 array of neighbor positions
            pbc_box: Optional PBC box

        Returns:
            1D array of distances
        """
        if len(neighbor_positions) == 0:
            return np.array([])

        if pbc_box is None:
            # No PBC
            delta = neighbor_positions - atom_position
            return np.linalg.norm(delta, axis=1).astype(np.float32)

        # With PBC
        box_size = np.array([
            pbc_box[1][0] - pbc_box[0][0],
            pbc_box[1][1] - pbc_box[0][1],
            pbc_box[1][2] - pbc_box[0][2]
        ])

        delta = neighbor_positions - atom_position
        # Apply minimum image convention
        delta = delta - box_size * np.round(delta / box_size)

        return np.linalg.norm(delta, axis=1).astype(np.float32)

    @staticmethod
    def build_bond_mask(
        elements: list,
        bond_criteria: dict,
        default_multiplier: float = 1.5
    ) -> np.ndarray:
        """
        Build bond distance threshold matrix.

        Args:
            elements: List of element symbols
            bond_criteria: Dict mapping (elem1, elem2) to threshold
            default_multiplier: Default multiplier if pair not in dict

        Returns:
            N x N matrix of bond distance thresholds
        """
        from phlox.base.constants import ATOMIC_RADII

        n = len(elements)
        mask = np.zeros((n, n), dtype=np.float32)

        for i in range(n):
            for j in range(i + 1, n):
                elem_i = elements[i]
                elem_j = elements[j]

                # Try to get specific threshold
                if elem_i in bond_criteria and elem_j in bond_criteria[elem_i]:
                    threshold = bond_criteria[elem_i][elem_j]
                elif elem_j in bond_criteria and elem_i in bond_criteria[elem_j]:
                    threshold = bond_criteria[elem_j][elem_i]
                else:
                    # Use default: sum of covalent radii * multiplier
                    radius_i = ATOMIC_RADII.get(elem_i, 1.0)
                    radius_j = ATOMIC_RADII.get(elem_j, 1.0)
                    threshold = (radius_i + radius_j) * default_multiplier

                mask[i, j] = threshold
                mask[j, i] = threshold

        return mask
