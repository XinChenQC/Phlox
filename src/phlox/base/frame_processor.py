"""
Shared frame processing logic for both serial and parallel analysis.
"""

import logging
import time
from typing import List, Optional, Dict, Any, Tuple
from collections import Counter

from phlox.base.config import Config
from phlox.base.atom import Atom
from phlox.base.molecule import Molecule
from phlox.frag.fragmenter import Fragmenter

logger = logging.getLogger(__name__)


class FrameProcessor:
    """
    Handles frame-by-frame molecular fragment identification and tracking.

    Shared by both serial analysis (Analyzer) and parallel workers (Worker).
    """

    def __init__(self, config: Config, worker_id: Optional[int] = None):
        """
        Initialize frame processor.

        Args:
            config: Configuration object
            worker_id: Optional worker identifier (for logging in parallel mode)
        """
        self.config = config
        self.worker_id = worker_id
        self.fragmenter: Optional[Fragmenter] = None
        self.is_first_frame = True

    def process_frame(
        self,
        atoms: List[Atom],
        pbc_box: Optional[tuple],
        frame_number: int,
        struct_dict: Dict[str, Any],
        molecule_info: Dict[str, List[int]],
    ) -> Tuple[List[Molecule], List[str]]:
        """
        Process a single trajectory frame.

        Steps:
        1. Initialize fragmenter on first frame
        2. Identify molecular fragments
        3. Update molecular properties
        4. Update struct_dict and molecule_info
        5. Track species for time series

        Args:
            atoms: List of atoms
            pbc_box: PBC box dimensions
            frame_number: Global frame number
            struct_dict: Dictionary to store species structures (DicStruct)
            molecule_info: Dictionary to store molecule instances (DicMoleInfo)

        Returns:
            Tuple of (molecules, frame_species_ids) where:
            - molecules: List of Molecule objects
            - frame_species_ids: List of species IDs for Counter/time series
        """
        log_prefix = f"Worker {self.worker_id}: " if self.worker_id is not None else ""
        logger.debug(f"{log_prefix}processing frame {frame_number}")

        # Initialize fragmenter on first frame
        if self.fragmenter is None:
            self.fragmenter = Fragmenter(self.config)
            self.fragmenter.initialize_mask_matrix(atoms)
            self.is_first_frame = True
        else:
            self.is_first_frame = False

        # Identify molecules
        if self.is_first_frame:
            # First frame - full distance calculation with neighbor search
            molecules, dist_mat = self.fragmenter.identify_molecules_initial(
                atoms, pbc_box, neighbor_cutoff=10.0
            )
        else:
            # Periodically rebuild neighbor lists to handle large atomic movements
            if frame_number % self.config.neighbor_update_freq == 0:
                molecules, dist_mat = self.fragmenter.identify_molecules_initial(
                    atoms, pbc_box, neighbor_cutoff=10.0
                )
            else:
                # Subsequent frames - incremental update using neighbor lists
                molecules = self.fragmenter.identify_molecules_update(atoms, pbc_box)

        # Update molecular properties for each molecule
        t_update_start = time.time()

        # Track species IDs for this frame (for Counter)
        frame_species_ids = []

        for molecule in molecules:
            self.fragmenter.update_molecule_info(molecule, atoms, pbc_box)

            # speciesID = "H" + molecular structure hash (hashD in ReaxANA)
            species_id = molecule.hash

            # Add to struct_dict (DicStruct) keyed by speciesID
            if species_id not in struct_dict:
                struct_dict[species_id] = {
                    'elements': molecule.elements,
                    'coordinates': molecule.coordinates,
                    'connectivity': molecule.connectivity,
                    'smiles_computed': False,
                    'smiles': '',
                    'formula': molecule.formula
                }

            # moleculeID = "S" + fragment_hash (unique instance identifier)
            # fragment_hash = SHA1(atom_indices + speciesID)
            molecule_id = "S" + molecule.fragment_hash

            # Add to molecule_info (DicMoleInfo) keyed by moleculeID
            if molecule_id not in molecule_info:
                molecule_info[molecule_id] = molecule.atom_indices

            # Collect species ID for counting
            frame_species_ids.append(species_id)

        # Log summary (show formulas instead of SMILES since SMILES not computed yet)
        unique_formulas = set(mol.formula for mol in molecules if mol.formula)
        elapsed = time.time() - t_update_start
        logger.info(f"{log_prefix}[FRAG] Frame {frame_number}: {len(molecules)} molecules, "
                   f"{len(unique_formulas)} unique formulas in {elapsed:.3f}s")
        if unique_formulas and len(unique_formulas) <= 20:  # Only show if not too many
            logger.info(f"{log_prefix}[FRAG] Formulas: {sorted(unique_formulas)}")

        return molecules, frame_species_ids

    def update_neighbor_lists(
        self,
        atoms: List[Atom],
        pbc_box: Optional[tuple],
        cutoff: float = 10.0
    ):
        """
        Update atomic neighbor lists.

        Called periodically to refresh neighbor lists as molecules move.
        Based on BlockNeighborUpdate() from tool.py (lines 424-472)

        Args:
            atoms: List of atoms
            pbc_box: PBC box dimensions
            cutoff: Distance cutoff for neighbors (Angstroms)
        """
        import numpy as np
        from phlox.frag.distance import DistanceCalculator

        n_atoms = len(atoms)
        positions = np.array([atom.position for atom in atoms], dtype=np.float32)

        # Calculate distance matrix
        calc = DistanceCalculator()
        if pbc_box is not None:
            dist_mat = calc.pairwise_distances(positions, pbc_box)
        else:
            dist_mat = calc.pairwise_distances(positions, None)

        # Update neighbor lists
        for i in range(n_atoms):
            atoms[i].clear_neighbors()
            for j in range(i):
                if dist_mat[i, j] < cutoff:
                    atoms[i].add_neighbor(j)
                    atoms[j].add_neighbor(i)
