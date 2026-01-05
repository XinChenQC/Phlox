"""
Chemistry utility functions.
"""

from collections import Counter
from typing import List, Dict


def molecular_formula(elements: List[str]) -> str:
    """
    Generate molecular formula from element list.

    Args:
        elements: List of element symbols

    Returns:
        Molecular formula (e.g., 'C2H4O')

    Example:
        >>> molecular_formula(['C', 'C', 'H', 'H', 'H', 'H', 'O'])
        'C2H4O'
    """
    counts = Counter(elements)

    # Order: C, H, then alphabetical
    formula = ""

    # Carbon first
    if 'C' in counts:
        count = counts['C']
        formula += 'C' if count == 1 else f'C{count}'

    # Hydrogen second
    if 'H' in counts:
        count = counts['H']
        formula += 'H' if count == 1 else f'H{count}'

    # Other elements alphabetically
    for element in sorted(counts.keys()):
        if element not in ['C', 'H']:
            count = counts[element]
            formula += element if count == 1 else f'{element}{count}'

    return formula


def element_count(elements: List[str]) -> Dict[str, int]:
    """
    Count occurrences of each element.

    Args:
        elements: List of element symbols

    Returns:
        Dictionary mapping element to count

    Example:
        >>> element_count(['C', 'C', 'H', 'H', 'H', 'H', 'O'])
        {'C': 2, 'H': 4, 'O': 1}
    """
    return dict(Counter(elements))


def mass_from_formula(formula: str) -> float:
    """
    Calculate molecular mass from formula.

    Args:
        formula: Molecular formula

    Returns:
        Molecular mass in amu
    """
    from phlox.base.constants import ATOMIC_MASSES

    # TODO: Parse formula and sum masses
    # This is a simplified placeholder
    return 0.0
