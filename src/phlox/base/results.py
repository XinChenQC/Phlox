"""
Results container for analysis outputs.
"""

from dataclasses import dataclass, field
from typing import Dict, List
from collections import Counter
import networkx as nx

from phlox.base.config import Config


@dataclass
class Results:
    """
    Container for analysis results.

    Attributes:
        config: Configuration used for analysis
        graph: Reaction network graph
        struct_dict: Dictionary of molecular structures (DicStruct in ReaxANA)
                     Keyed by speciesID = "H" + hash(structure)
                     Format: {speciesID: {'elements': [...], 'coordinates': [...],
                                         'smiles_computed': bool, 'smiles': str, 'formula': str}}
        molecule_info: Dictionary of molecular instances (DicMoleInfo in ReaxANA)
                       Keyed by moleculeID = "S" + SHA1(atom_indices + speciesID)
                       Format: {moleculeID: atom_indices}
        reactions: Dictionary of classified reactions
        time_series: Species count per frame (SpeciesCount in ReaxANA)
                     List of Counters, one Counter per frame
                     Format: [Counter({speciesID: count}), ...]
        total_frames: Total number of frames processed
    """
    config: Config
    graph: nx.MultiDiGraph = field(default_factory=nx.MultiDiGraph)
    struct_dict: Dict = field(default_factory=dict)  # DicStruct: speciesID -> structure info
    molecule_info: Dict = field(default_factory=dict)  # DicMoleInfo: moleculeID -> atom indices
    reactions: Dict = field(default_factory=dict)
    time_series: List[Counter] = field(default_factory=list)  # SpeciesCount: [Counter, ...]
    total_frames: int = 0

    def to_dict(self) -> dict:
        """Convert results to dictionary."""
        return {
            'config': self.config.to_dict(),
            'total_frames': self.total_frames,
            'n_species': len(self.struct_dict),
            'n_molecules': len(self.molecule_info),
            'n_reactions': len(self.reactions),
            'n_nodes': self.graph.number_of_nodes(),
            'n_edges': self.graph.number_of_edges(),
        }

    def summary(self) -> str:
        """Generate summary string."""
        return f"""
Phlox Analysis Results
=====================
Total frames: {self.total_frames}
Unique species: {len(self.struct_dict)}
Molecular instances: {len(self.molecule_info)}
Reaction types: {len(self.reactions)}
Graph nodes: {self.graph.number_of_nodes()}
Graph edges: {self.graph.number_of_edges()}
Output directory: {self.config.output_dir}
"""
