"""
Physical and chemical constants used throughout the package.
"""

from typing import Dict

# Element name to atomic number mapping
ELEMENTS: Dict[str, int] = {
    'GHOST': 0, 'H': 1, 'He': 2, 'Li': 3, 'Be': 4,
    'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9,
    'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14,
    'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19,
    'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24,
    'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29,
    'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34,
    'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39,
    'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44,
    'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49,
    'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54,
    'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59,
    'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64,
    'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69,
    'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74,
    'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79,
    'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84,
    'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89,
    'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94,
    'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99,
    'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104,
    'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109,
}

# Atomic radii in Angstroms (covalent radii)
ATOMIC_RADII: Dict[str, float] = {
    'H': 0.31, 'He': 0.28,
    'Li': 1.21, 'Be': 0.96, 'B': 0.84, 'C': 0.69, 'N': 0.71,
    'O': 0.66, 'F': 0.64, 'Ne': 0.58,
    'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 'P': 1.07,
    'S': 1.05, 'Cl': 1.02, 'Ar': 1.06,
    'K': 2.03, 'Ca': 1.76, 'Sc': 1.70, 'Ti': 1.60, 'V': 1.53,
    'Cr': 1.39, 'Mn': 1.39, 'Fe': 1.32, 'Co': 1.26, 'Ni': 1.24,
    'Cu': 1.32, 'Zn': 1.22
}

# Atomic masses in amu
ATOMIC_MASSES: Dict[str, float] = {
    'H': 1.008, 'He': 4.003,
    'Li': 6.941, 'Be': 9.012, 'B': 10.811, 'C': 12.011, 'N': 14.007,
    'O': 15.999, 'F': 18.998, 'Ne': 20.180,
    'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.086, 'P': 30.974,
    'S': 32.065, 'Cl': 35.453, 'Ar': 39.948,
    'K': 39.098, 'Ca': 40.078, 'Sc': 44.956, 'Ti': 47.867, 'V': 50.942,
    'Cr': 51.996, 'Mn': 54.938, 'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693,
    'Cu': 63.546, 'Zn': 65.409
}

# Bond criteria multiplier for radii sum
BOND_CRITERIA: Dict[str, Dict[str, float]] = {
    'C': {'C': 1.34, 'H': 1.07, 'O': 1.27, 'Ni': 1.91},
    'H': {'C': 1.07, 'H': 0.80, 'O': 0.97, 'Ni': 1.64},
    'O': {'C': 1.27, 'H': 0.97, 'O': 1.14, 'Ni': 1.81},
    'Ni': {'C': 1.27, 'H': 0.97, 'O': 1.14, 'Ni': 2.48},
}

# Default bond multiplier if specific pair not found
DEFAULT_BOND_MULTIPLIER: float = 1.5

# Neighbor search cutoff (Angstroms)
NEIGHBOR_CUTOFF: float = 10.0

# Block neighbor update cutoff (Angstroms)
BLOCK_NEIGHBOR_CUTOFF: float = 5.0
