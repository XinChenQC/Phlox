"""
Batch SMILES generation for all species.

Ported from printUnknowStruc() in ReaxANA/tool.py
"""

import logging
from typing import Dict, Any

logger = logging.getLogger(__name__)


def compute_species_smiles(species_dict: Dict[str, Dict[str, Any]]) -> int:
    """
    Compute SMILES for all species in batch.

    Ported from printUnknowStruc() in tool.py (lines 742-766)

    This is called AFTER trajectory processing to compute SMILES for all
    unique species in batch, which is much more efficient than computing
    SMILES for every molecule in every frame.

    Args:
        species_dict: Dictionary mapping species hash -> species info
                     Each entry must contain:
                     - elements: list of element symbols
                     - coordinates: list of atomic coordinates
                     - formula: molecular formula
                     - smiles_computed: boolean flag
                     - smiles: SMILES string (empty if not computed)

    Returns:
        Number of species for which SMILES was computed

    Side effects:
        Updates species_dict in place:
        - Sets smiles_computed = True
        - Sets smiles = computed SMILES or formula (for large molecules)
    """
    from phlox.utils.smiles import generate_smiles

    count = 0
    total_species = len(species_dict)

    logger.info(f"Computing SMILES for {total_species} unique species...")

    for hash_key, species_data in species_dict.items():
        # Skip if already computed
        if species_data.get('smiles_computed', False):
            continue

        # Mark as computed
        species_data['smiles_computed'] = True
        count += 1

        elements = species_data['elements']
        coordinates = species_data['coordinates']
        formula = species_data['formula']
        connectivity = species_data.get('connectivity', None)

        # Check if too many atoms (>30 heavy atoms)
        heavy_atoms = sum(1 for elem in elements if elem != 'H')
        if heavy_atoms > 30:
            # Use formula for large molecules
            species_data['smiles'] = formula
            continue

        # Generate SMILES
        smiles = generate_smiles(elements, coordinates, connectivity)

        # Post-processing for NO2 groups (ReaxANA approach, tool.py lines 760-763)
        if smiles:
            if "N([O])[O]" in smiles and not smiles.startswith("N([O])[O]") and ".N([O])[O]" not in smiles:
                smiles = smiles.replace("N([O])[O]", "N(=O)=O")
            if "[O]N([O])" in smiles and not smiles.startswith("[O]N([O])") and ".[O]N([O])" not in smiles:
                smiles = smiles.replace("[O]N([O])", "O=N(=O)")

        # Store SMILES (or use formula if generation failed)
        species_data['smiles'] = smiles if smiles else formula

        # Progress logging
        if count % 100 == 0:
            logger.info(f"  Computed SMILES for {count}/{total_species} species...")

    logger.info(f"SMILES computation complete: {count} species processed")

    return count
