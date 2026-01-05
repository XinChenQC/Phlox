"""
Reaction network builder.

Constructs reaction network by comparing consecutive trajectory frames.

Ports functionality from ReaxANA/tool.py::compare2Step()
"""

import logging
import networkx as nx
from typing import List, Set, Tuple
from phlox.base.molecule import Molecule

logger = logging.getLogger(__name__)


class NetworkBuilder:
    """
    Build reaction network from molecular trajectories.

    Compares consecutive frames to detect reactions and builds directed graph.

    Based on compare2Step() from ReaxANA/tool.py
    """

    def __init__(self, config):
        """
        Initialize network builder.

        Args:
            config: Configuration object
        """
        self.config = config
        self.graph = nx.MultiDiGraph()

    def compare_frames(
        self,
        molecules_prev: List[Molecule],
        molecules_curr: List[Molecule],
        timestep: int
    ):
        """
        Compare two consecutive frames to detect reactions.

        Ported from compare2Step() in tool.py (lines 465-601)

        Args:
            molecules_prev: Molecules at timestep - 1
            molecules_curr: Molecules at timestep
            timestep: Current timestep

        Algorithm:
        1. Get fragment hashes from both frames
        2. Find common, disappearing, appearing fragments
        3. Match by atom overlap
        4. Classify reaction type
        5. Add edges to graph
        """
        logger.debug(f"Comparing frames at timestep {timestep}")

        # Extract fragment hashes (index 6 in original blk arrays)
        hashindex_prev = [mol.fragment_hash for mol in molecules_prev]
        hashindex_curr = [mol.fragment_hash for mol in molecules_curr]

        set_prev = set(hashindex_prev)
        set_curr = set(hashindex_curr)

        # Find common and changing fragments
        common_part = set_prev & set_curr
        from_hashes = set_prev - common_part  # Disappearing (reactants)
        to_hashes = set_curr - common_part    # Appearing (products)

        if not from_hashes and not to_hashes:
            return  # No reaction

        logger.debug(f"  Reactants: {len(from_hashes)}, Products: {len(to_hashes)}")

        # Build mappings from hash to molecule and atom sets
        from_atom_sets = []
        from_indices = []
        from_molecules = []

        for fhash in from_hashes:
            idx = hashindex_prev.index(fhash)
            mol = molecules_prev[idx]
            from_atom_sets.append(set(mol.atom_indices))
            from_indices.append(idx)
            from_molecules.append(mol)

        to_atom_sets = []
        to_indices = []
        to_molecules = []

        for fhash in to_hashes:
            idx = hashindex_curr.index(fhash)
            mol = molecules_curr[idx]
            to_atom_sets.append(set(mol.atom_indices))
            to_indices.append(idx)
            to_molecules.append(mol)

        # Build temporary graph to find connected reactions
        G_temp = nx.DiGraph()

        for i in range(len(from_atom_sets)):
            for j in range(len(to_atom_sets)):
                overlap = from_atom_sets[i] & to_atom_sets[j]
                if overlap:
                    # Add edge from reactant to product
                    G_temp.add_edge(f"R{from_indices[i]}", f"P{to_indices[j]}")

        # Find connected components (separate reactions)
        components = list(nx.weakly_connected_components(G_temp))

        # Process each separate reaction
        for component in components:
            r_list = []
            p_list = []
            n_r = 0
            n_p = 0

            for node in component:
                if node[0] == 'R':
                    n_r += 1
                    r_list.append(int(node[1:]))
                elif node[0] == 'P':
                    n_p += 1
                    p_list.append(int(node[1:]))

            # Classify and add reaction to graph
            if n_r > 1 and n_p > 1:
                # Multi-to-multi reaction
                self._add_multi_reaction(
                    r_list, p_list, timestep,
                    molecules_prev, molecules_curr
                )
            elif n_r > 1 and n_p == 1:
                # Combination reaction
                self._add_combination(
                    r_list, p_list[0], timestep,
                    molecules_prev, molecules_curr
                )
            elif n_r == 1 and n_p > 1:
                # Dissociation reaction
                self._add_dissociation(
                    r_list[0], p_list, timestep,
                    molecules_prev, molecules_curr
                )
            elif n_r == 1 and n_p == 1:
                # Isomerization reaction
                self._add_isomerization(
                    r_list[0], p_list[0], timestep,
                    molecules_prev, molecules_curr
                )

    def classify_reaction_type(
        self,
        reactants: List[Molecule],
        products: List[Molecule]
    ) -> str:
        """
        Classify reaction based on stoichiometry.

        Args:
            reactants: List of reactant molecules
            products: List of product molecules

        Returns:
            Reaction type: 'combination', 'dissociation', 'isomerization', 'multi'
        """
        n_reactants = len(reactants)
        n_products = len(products)

        if n_reactants == 1 and n_products == 1:
            return 'isomerization'
        elif n_reactants > 1 and n_products == 1:
            return 'combination'
        elif n_reactants == 1 and n_products > 1:
            return 'dissociation'
        else:
            return 'multi'

    def _add_combination(
        self,
        r_list: List[int],
        p_idx: int,
        timestep: int,
        molecules_prev: List[Molecule],
        molecules_curr: List[Molecule]
    ):
        """
        Add combination reaction edge: A + B → C

        Ported from tool.py lines 562-571

        Args:
            r_list: List of reactant indices
            p_idx: Product index
            timestep: Timestep when reaction occurred
            molecules_prev: Previous frame molecules
            molecules_curr: Current frame molecules
        """
        # Product node
        prod_mol = molecules_curr[p_idx]
        hash_p = "S" + prod_mol.fragment_hash

        self.graph.add_node(
            hash_p,
            hashD=prod_mol.hash,
            label=prod_mol.formula,
            smiles=prod_mol.smiles,
            atom_indices=prod_mol.atom_indices
        )

        # Reactant nodes and edges
        for r_idx in r_list:
            react_mol = molecules_prev[r_idx]
            hash_r = "S" + react_mol.fragment_hash

            self.graph.add_node(
                hash_r,
                hashD=react_mol.hash,
                label=react_mol.formula,
                smiles=react_mol.smiles,
                atom_indices=react_mol.atom_indices
            )

            self.graph.add_edge(
                hash_r,
                hash_p,
                label=str(timestep),
                color='red'
            )

    def _add_dissociation(
        self,
        r_idx: int,
        p_list: List[int],
        timestep: int,
        molecules_prev: List[Molecule],
        molecules_curr: List[Molecule]
    ):
        """
        Add dissociation reaction edge: C → A + B

        Ported from tool.py lines 574-585

        Args:
            r_idx: Reactant index
            p_list: List of product indices
            timestep: Timestep when reaction occurred
            molecules_prev: Previous frame molecules
            molecules_curr: Current frame molecules
        """
        # Reactant node
        react_mol = molecules_prev[r_idx]
        hash_r = "S" + react_mol.fragment_hash

        self.graph.add_node(
            hash_r,
            hashD=react_mol.hash,
            label=react_mol.formula,
            smiles=react_mol.smiles,
            atom_indices=react_mol.atom_indices
        )

        # Product nodes and edges
        for p_idx in p_list:
            prod_mol = molecules_curr[p_idx]
            hash_p = "S" + prod_mol.fragment_hash

            self.graph.add_node(
                hash_p,
                hashD=prod_mol.hash,
                label=prod_mol.formula,
                smiles=prod_mol.smiles,
                atom_indices=prod_mol.atom_indices
            )

            self.graph.add_edge(
                hash_r,
                hash_p,
                label=str(timestep),
                color='blue'
            )

    def _add_isomerization(
        self,
        r_idx: int,
        p_idx: int,
        timestep: int,
        molecules_prev: List[Molecule],
        molecules_curr: List[Molecule]
    ):
        """
        Add isomerization reaction edge: A → A'

        Ported from tool.py lines 588-599

        Args:
            r_idx: Reactant index
            p_idx: Product index
            timestep: Timestep when reaction occurred
            molecules_prev: Previous frame molecules
            molecules_curr: Current frame molecules
        """
        react_mol = molecules_prev[r_idx]
        prod_mol = molecules_curr[p_idx]

        hash_r = "S" + react_mol.fragment_hash
        hash_p = "S" + prod_mol.fragment_hash

        self.graph.add_node(
            hash_r,
            hashD=react_mol.hash,
            label=react_mol.formula,
            smiles=react_mol.smiles,
            atom_indices=react_mol.atom_indices
        )

        self.graph.add_node(
            hash_p,
            hashD=prod_mol.hash,
            label=prod_mol.formula,
            smiles=prod_mol.smiles,
            atom_indices=prod_mol.atom_indices
        )

        self.graph.add_edge(
            hash_r,
            hash_p,
            label=str(timestep),
            color='grey'
        )

    def _add_multi_reaction(
        self,
        r_list: List[int],
        p_list: List[int],
        timestep: int,
        molecules_prev: List[Molecule],
        molecules_curr: List[Molecule]
    ):
        """
        Add multi-molecular reaction: A + B → C + D

        Creates intermediate node based on hash difference.

        Ported from tool.py lines 528-559

        Args:
            r_list: List of reactant indices
            p_list: List of product indices
            timestep: Timestep when reaction occurred
            molecules_prev: Previous frame molecules
            molecules_curr: Current frame molecules
        """
        # Calculate intermediate node hash
        # Sum reactant hashes, subtract product hashes
        num = 0

        for r_idx in r_list:
            mol = molecules_prev[r_idx]
            fhash = mol.fragment_hash
            # Handle special prefixes (CC_, NC_)
            if "CC_" in fhash or "NC_" in fhash:
                num += int(fhash[3:], 16)
            else:
                num += int(fhash, 16)

        for p_idx in p_list:
            mol = molecules_curr[p_idx]
            fhash = mol.fragment_hash
            if "CC_" in fhash or "NC_" in fhash:
                num -= int(fhash[3:], 16)
            else:
                num -= int(fhash, 16)

        # Create intermediate node
        inter_index = "int-" + hex(num)
        inter_label = "int-" + hex(num)[0:5]

        self.graph.add_node(
            inter_index,
            label=inter_label,
            hashD=inter_label
        )

        # Add reactant → intermediate edges
        for r_idx in r_list:
            react_mol = molecules_prev[r_idx]
            hash_r = "S" + react_mol.fragment_hash

            self.graph.add_node(
                hash_r,
                hashD=react_mol.hash,
                label=react_mol.formula,
                smiles=react_mol.smiles,
                atom_indices=react_mol.atom_indices
            )

            self.graph.add_edge(
                hash_r,
                inter_index,
                label=str(timestep),
                color='red'
            )

        # Add intermediate → product edges
        for p_idx in p_list:
            prod_mol = molecules_curr[p_idx]
            hash_p = "S" + prod_mol.fragment_hash

            self.graph.add_node(
                hash_p,
                hashD=prod_mol.hash,
                label=prod_mol.formula,
                smiles=prod_mol.smiles,
                atom_indices=prod_mol.atom_indices
            )

            self.graph.add_edge(
                inter_index,
                hash_p,
                label=str(timestep),
                color='blue'
            )
