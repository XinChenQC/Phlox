"""
SMILES generation and manipulation.

Uses xyz2mol algorithm from ReaxANA for proper bond order determination.
"""

import logging
import time
from typing import List, Optional

logger = logging.getLogger(__name__)


def generate_smiles(elements: List[str], coordinates: List[List[float]], connectivity=None) -> str:
    """
    Generate SMILES string from molecular structure using xyz2mol approach.

    This properly determines bond orders (single, double, triple) from 3D geometry,
    unlike simple connectivity which assumes all bonds are single.

    Args:
        elements: List of element symbols
        coordinates: List of [x, y, z] coordinates
        connectivity: Connectivity matrix (numpy array) - optional, will be computed if not provided

    Returns:
        SMILES string, or empty string if generation fails or molecule is too large
    """
    t_start = time.time()
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit import RDLogger
        import numpy as np

        # Disable RDKit warnings
        RDLogger.DisableLog('rdApp.info')

        # Count heavy atoms (non-hydrogen)
        heavy_atom_count = sum(1 for elem in elements if elem != 'H')

        # Skip SMILES generation for large molecules (>30 heavy atoms)
        if heavy_atom_count > 30:
            return ""

        # Use xyz2mol approach from ReaxANA (trans_smile.py)
        try:
            from phlox.utils.xyz2mol import xyz2mol, int_atom

            # Convert element symbols to atomic numbers
            atoms = [int_atom(elem) for elem in elements]

            # Use xyz2mol to get proper bond orders
            mols = xyz2mol(
                atoms,
                coordinates,
                charge=0,
                use_graph=True,
                allow_charged_fragments=False,
                embed_chiral=False,
                use_huckel=True
            )

            if not mols:
                return ""

            # Get first mol and generate SMILES
            mol = mols[0]

            # Remove explicit hydrogens
            mol = Chem.RemoveHs(mol)

            # Generate SMILES
            smiles = Chem.MolToSmiles(mol, isomericSmiles=True)

            # Canonical hack - parse and regenerate for consistency
            m = Chem.MolFromSmiles(smiles)
            if m:
                smiles = Chem.MolToSmiles(m, isomericSmiles=True)

            return smiles

        except Exception as e:
            logger.debug(f"[SMILES] Failed: {e}")
            return ""

    except ImportError as e:
        logger.debug(f"[SMILES] Import failed: {e}")
        return ""
    except Exception as e:
        logger.debug(f"[SMILES] Generation failed: {e}")
        return ""


def canonicalize_smiles(smiles: str) -> str:
    """
    Canonicalize SMILES string.

    Args:
        smiles: Input SMILES string

    Returns:
        Canonical SMILES string
    """
    try:
        from rdkit import Chem

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return smiles

        return Chem.MolToSmiles(mol)

    except ImportError:
        logger.error("RDKit not installed. Cannot canonicalize SMILES.")
        return smiles
    except Exception as e:
        logger.warning(f"SMILES canonicalization failed: {e}")
        return smiles


def smiles_to_formula(smiles: str) -> str:
    """
    Convert SMILES to molecular formula.

    Args:
        smiles: SMILES string

    Returns:
        Molecular formula (e.g., 'C2H4O')
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return ""

        return Descriptors.MolecularFormula(mol)

    except ImportError:
        logger.error("RDKit not installed.")
        return ""
    except Exception as e:
        logger.warning(f"Formula generation failed: {e}")
        return ""
