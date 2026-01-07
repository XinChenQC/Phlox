"""
Molecular fragmentation engine.

Ports functionality from ReaxANA/tool.py
"""

import logging
import time
from typing import List, Tuple, Optional
import numpy as np
import networkx as nx
from collections import Counter
import gc

from phlox.base.atom import Atom
from phlox.base.molecule import Molecule
from phlox.base.config import Config
from phlox.frag.distance import DistanceCalculator
from phlox.frag.union_find import UnionFind

logger = logging.getLogger(__name__)


class Fragmenter:
    """
    Molecular fragmentation engine.

    Identifies molecular fragments from atomic coordinates and builds
    molecular connectivity.

    Ported from ReaxANA/tool.py
    """

    def __init__(self, config: Config):
        """
        Initialize fragmenter.

        Args:
            config: Configuration object
        """
        self.config = config
        self.distance_calc = DistanceCalculator()

        # Global matrices (reused across frames)
        self.global_dist_mat = None
        self.global_mask_mat = None

        # Store atoms across frames to preserve neighbor lists
        self.atoms = None

    def initialize_mask_matrix(self, atoms: List[Atom]):
        """
        Initialize global bond distance mask matrix.

        Based on BuildMaskForXYZ/BuildMaskForLAMMPS from tool.py

        Args:
            atoms: List of atoms
        """
        n_atoms = len(atoms)
        elements = [atom.element for atom in atoms]

        # Build radii matrix (ReaxANA approach: tool.py lines 40-42)
        from phlox.base.constants import ATOMIC_RADII

        radii_array = np.array([ATOMIC_RADII.get(elem, 1.0) for elem in elements], dtype=np.float32)
        radii_array = np.tile(radii_array, (n_atoms, 1))

        # Build mask matrix: (radii[i] + radii[j]) * multiplier
        # This matches ReaxANA's BuildMaskForXYZ approach
        self.global_mask_mat = (radii_array + radii_array.T) * self.config.bond_multiplier
        np.fill_diagonal(self.global_mask_mat, 0)

        # Initialize distance matrix
        self.global_dist_mat = np.zeros((n_atoms, n_atoms), dtype=np.float32)

        logger.debug(f"Initialized mask matrix for {n_atoms} atoms")

    def identify_molecules_initial(
        self,
        atoms: List[Atom],
        pbc_box: Optional[tuple],
        neighbor_cutoff: float = 10.0
    ) -> Tuple[List[Molecule], np.ndarray]:
        """
        Identify molecular fragments (initial frame with full distance calc).

        Ported from buildNeigh_AtomicBased() in tool.py (lines 196-245)

        Args:
            atoms: List of Atom objects
            pbc_box: PBC box dimensions
            neighbor_cutoff: Cutoff for neighbor search (Angstroms)

        Returns:
            Tuple of (list of Molecule objects, distance matrix)
        """
        start_time = time.time()

        # Store atoms for reuse across frames
        self.atoms = atoms
        n_atoms = len(atoms)

        # Extract positions
        positions = np.array([atom.position for atom in atoms], dtype=np.float32)

        # Calculate full distance matrix
        if pbc_box is not None:
            box_size = (
                pbc_box[1][0] - pbc_box[0][0],
                pbc_box[1][1] - pbc_box[0][1],
                pbc_box[1][2] - pbc_box[0][2]
            )
            self.global_dist_mat = self.distance_calc.pairwise_distances(positions, pbc_box)
        else:
            self.global_dist_mat = self.distance_calc.pairwise_distances(positions, None)

        logger.debug(f"Distance matrix calculated in {time.time() - start_time:.3f}s")

        # Build neighbor mask (atoms within cutoff)
        neighbor_mask = self.global_dist_mat < neighbor_cutoff

        # Build connectivity mask (atoms bonded)
        link_mat = self.global_dist_mat < self.global_mask_mat

        # Update atom neighbor lists (only store j < i to avoid duplication)
        for i in range(n_atoms):
            link_mat[i][i] = False
            atoms[i].clear_neighbors()
            for j in range(i):
                if neighbor_mask[i, j]:
                    atoms[i].add_neighbor(j)

        # Use UnionFind to identify fragments
        uf = UnionFind.from_connectivity_matrix(link_mat)
        components = uf.get_components()

        # Clean up
        del neighbor_mask
        del link_mat
        gc.collect()

        # Create Molecule objects
        molecules = []
        for mol_id, atom_indices in enumerate(components):
            mol = Molecule(
                mol_id=mol_id,
                atom_indices=sorted(list(atom_indices))
            )
            molecules.append(mol)

        logger.debug(f"Identified {len(molecules)} molecules in {time.time() - start_time:.3f}s")

        return molecules, self.global_dist_mat

    def identify_molecules_update(
        self,
        atoms: List[Atom],
        pbc_box: Optional[tuple]
    ) -> List[Molecule]:
        """
        Update molecular fragments (subsequent frames with partial distance calc).

        Ported from buildDistMart() in tool.py (lines 247-279)

        Args:
            atoms: List of Atom objects with new positions
            pbc_box: PBC box dimensions

        Returns:
            List of Molecule objects
        """
        start_time = time.time()

        # Update stored atoms with new positions
        if self.atoms is None:
            raise RuntimeError("Atoms not initialized. Call identify_molecules_initial() first.")

        for i, atom in enumerate(atoms):
            self.atoms[i].position = atom.position

        # Update distances only for neighbors
        if pbc_box is not None:
            box_size = (
                pbc_box[1][0] - pbc_box[0][0],
                pbc_box[1][1] - pbc_box[0][1],
                pbc_box[1][2] - pbc_box[0][2]
            )
        else:
            box_size = None
        for idx, atom in enumerate(self.atoms):
            if not atom.neighbors:
                continue

            # Get neighbor positions
            neighbor_positions = np.array([self.atoms[j].position for j in atom.neighbors], dtype=np.float32)

            # Calculate distances
            distances = self.distance_calc.neighbor_distances(atom.position, neighbor_positions, pbc_box)
            # Update global distance matrix
            for j_idx, j in enumerate(atom.neighbors):
                self.global_dist_mat[idx][j] = distances[j_idx]
                self.global_dist_mat[j][idx] = distances[j_idx]

        # Build connectivity matrix
        link_mat = self.global_dist_mat < self.global_mask_mat

        # Use UnionFind to identify fragments
        uf = UnionFind.from_connectivity_matrix(link_mat)
        components = uf.get_components()

        del link_mat
        gc.collect()

        # Create Molecule objects
        molecules = []
        for mol_id, atom_indices in enumerate(components):
            mol = Molecule(
                mol_id=mol_id,
                atom_indices=sorted(list(atom_indices))
            )
            molecules.append(mol)

        logger.debug(f"Updated {len(molecules)} molecules in {time.time() - start_time:.3f}s")

        return molecules

    def update_molecule_info(
        self,
        molecule: Molecule,
        atoms: List[Atom],
        pbc_box: Optional[tuple],
        store_coordinates: bool = True
    ):
        """
        Update molecular properties (hash, formula, center, radius).

        Ported from BlockInfoUpdatePBC() and getBlkInfoPBC() in tool.py

        Args:
            molecule: Molecule object to update
            atoms: Full list of atoms
            pbc_box: PBC box dimensions
            store_coordinates: If True, store unwrapped coordinates in molecule (for later SMILES generation)
        """
        # Get atoms in molecule
        mol_atoms = [atoms[i] for i in molecule.atom_indices]
        mol_positions = np.array([a.position for a in mol_atoms], dtype=np.float32)
        mol_elements = [a.element for a in mol_atoms]

        # Remove PBC wrapping for internal coordinates
        if pbc_box is not None:
            box_size = np.array([
                pbc_box[1][0] - pbc_box[0][0],
                pbc_box[1][1] - pbc_box[0][1],
                pbc_box[1][2] - pbc_box[0][2]
            ])

            # Unwrap relative to first atom
            ref_pos = mol_positions[0]
            for i in range(1, len(mol_positions)):
                delta = mol_positions[i] - ref_pos
                for dim in range(3):
                    if abs(delta[dim]) > box_size[dim] / 2:
                        mol_positions[i, dim] -= box_size[dim] * np.sign(delta[dim])

        # Calculate internal distance matrix
        from scipy.spatial.distance import pdist, squareform
        internal_dist = squareform(pdist(mol_positions))

        # Build connectivity matrix for this molecule
        from phlox.base.constants import ATOMIC_RADII

        n_atoms_mol = len(mol_atoms)
        connectivity = np.zeros((n_atoms_mol, n_atoms_mol), dtype=bool)

        for i in range(n_atoms_mol):
            for j in range(i + 1, n_atoms_mol):
                r_i = ATOMIC_RADII.get(mol_elements[i], 1.0)
                r_j = ATOMIC_RADII.get(mol_elements[j], 1.0)
                threshold = (r_i + r_j) * self.config.bond_multiplier

                if internal_dist[i, j] < threshold:
                    connectivity[i, j] = True
                    connectivity[j, i] = True

        # Calculate molecular hash and formula
        mol_hash, formula = self.calculate_molecular_hash(connectivity, mol_elements)
        molecule.hash = mol_hash
        molecule.formula = formula

        # Store coordinates and elements for later SMILES generation (ReaxANA approach)
        if store_coordinates:
            molecule.elements = mol_elements
            molecule.coordinates = mol_positions.tolist()
            molecule.connectivity = connectivity

        # Don't compute SMILES during trajectory - will be computed later in batch
        # Use hash as placeholder (ReaxANA approach)
        molecule.smiles = mol_hash

        # Calculate center of mass
        molecule.center = np.mean(mol_positions, axis=0)

        # Calculate radius
        distances = np.linalg.norm(mol_positions - molecule.center, axis=1)
        molecule.radius = np.max(distances) if len(distances) > 0 else 0.0

        # Calculate fragment hash
        molecule.calculate_fragment_hash()

    def calculate_molecular_hash(
        self,
        connectivity: np.ndarray,
        elements: List[str]
    ) -> Tuple[str, str]:
        """
        Calculate Weisfeiler-Lehman hash and molecular formula.

        Ported from MolHash() in tool.py (lines 289-330)

        Args:
            connectivity: Boolean connectivity matrix
            elements: List of element symbols

        Returns:
            Tuple of (hash_string, formula_string)
        """
        # Build NetworkX graph
        G = nx.Graph()

        # Count elements
        element_counts = Counter(elements)

        # Add nodes
        for i in range(len(elements)):
            G.add_node(i, label=elements[i])

        # Add edges
        for i in range(len(connectivity)):
            for j in range(i + 1, len(connectivity)):
                if connectivity[i][j]:
                    G.add_edge(i, j)

        # Calculate Weisfeiler-Lehman hash
        hash_str = nx.weisfeiler_lehman_graph_hash(G, node_attr='label')

        # Build formula string (C, H, O, N, F order, then alphabetical)
        formula = ""
        for elem in ['C', 'H', 'O', 'N', 'F']:
            if elem in element_counts:
                count = element_counts[elem]
                formula += elem if count == 1 else f"{elem}{count}"

        # Add other elements
        for elem in sorted(element_counts.keys()):
            if elem not in ['C', 'H', 'O', 'N', 'F']:
                count = element_counts[elem]
                formula += elem if count == 1 else f"{elem}{count}"

        return "H" + hash_str, formula
