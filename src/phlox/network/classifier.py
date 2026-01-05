"""
Reaction classifier.

Classifies and counts different reaction types in the reaction network.

Ported from ReaxANA/tool_reactCount.py
"""

import logging
import hashlib
import networkx as nx
from typing import List, Tuple, Dict
from copy import deepcopy

logger = logging.getLogger(__name__)


def sort_by_index(sub_li: list, index: int) -> list:
    """Sort list by specific index."""
    sub_li.sort(key=lambda x: x[index])
    return sub_li
import networkx as nx
from typing import Dict, List, Tuple, Set
from collections import defaultdict

logger = logging.getLogger(__name__)


class ReactionClassifier:
    """
    Classify reactions from reaction network.

    Extracts discrete reaction events and groups them by type.

    Based on generateReactions() from ReaxANA/tool_reactCount.py
    """

    def __init__(self, stability_threshold: int = 10):
        """
        Initialize classifier.

        Args:
            stability_threshold: Time threshold for short/long reactions
        """
        self.stability_threshold = stability_threshold
        self.reactions = {}

    def classify_all(self, graph: nx.MultiDiGraph) -> Dict:
        """
        Classify all reactions in the network.

        Based on generateReactions() from tool_reactCount.py

        Args:
            graph: Reaction network

        Returns:
            Dictionary of reactions grouped by type
        """
        logger.info("Classifying reactions")

        # Build node information
        node_info = self.build_node_info(graph)

        # Find different reaction types
        short_reactions = self.find_short_reactions(graph, node_info)
        long_reactions = self.find_long_reactions(graph, node_info)
        combinations = self.find_combinations(graph, node_info)
        dissociations = self.find_dissociations(graph, node_info)
        isomerizations = self.find_isomerizations(graph)

        # Combine and classify
        all_reactions = (
            short_reactions + long_reactions +
            combinations + dissociations + isomerizations
        )

        # Group by reaction hash
        reaction_dict = self.group_by_hash(all_reactions)

        logger.info(f"Found {len(reaction_dict)} unique reaction types")
        logger.info(f"Total reaction events: {sum(len(v) for v in reaction_dict.values())}")

        self.reactions = reaction_dict
        return reaction_dict

    def build_node_info(self, graph: nx.MultiDiGraph) -> Dict:
        """
        Build information about each node.

        Based on buildBasicInfo() from tool_reactCount.py

        Returns:
            Dictionary mapping node_id to node information
        """
        node_info = {}

        for node in graph.nodes():
            # Get incoming edges
            in_edges = []
            for pred, _, data in graph.in_edges(node, data=True):
                timestep = int(data.get('label', 0))
                color = data.get('color', 'grey')
                in_edges.append((pred, color, timestep))

            # Get outgoing edges
            out_edges = []
            for _, succ, data in graph.out_edges(node, data=True):
                timestep = int(data.get('label', 0))
                color = data.get('color', 'grey')
                out_edges.append((succ, color, timestep))

            # Sort by timestep
            in_edges.sort(key=lambda x: x[2])
            out_edges.sort(key=lambda x: x[2])

            # Determine node type
            node_type = self.determine_node_type(in_edges, out_edges)

            node_info[node] = {
                'in_edges': in_edges,
                'out_edges': out_edges,
                'type': node_type
            }

        return node_info

    def determine_node_type(self, in_edges: List, out_edges: List) -> str:
        """
        Determine if node is 'center' or 'edge'.

        Center nodes have multiple simultaneous reactions.

        Args:
            in_edges: List of (node, color, timestep) tuples
            out_edges: List of (node, color, timestep) tuples

        Returns:
            'center' or 'edge'
        """
        # Check for multiple red edges at same timestep (combination)
        for i in range(len(in_edges) - 1):
            if (in_edges[i][2] == in_edges[i+1][2] and
                in_edges[i][1] == 'red' and in_edges[i+1][1] == 'red'):
                return 'center'

        # Check for multiple blue edges at same timestep (dissociation)
        for i in range(len(out_edges) - 1):
            if (out_edges[i][2] == out_edges[i+1][2] and
                out_edges[i][1] == 'blue' and out_edges[i+1][1] == 'blue'):
                return 'center'

        return 'edge'

    def find_short_reactions(self, graph: nx.MultiDiGraph, node_info: Dict) -> List:
        """
        Find fast reactions: A + B → C + D (single step).

        Based on SimpleNodeshortReactAbs() from tool_reactCount.py

        Returns:
            List of reaction records
        """
        reactions = []

        # TODO: Implement short reaction detection
        # Look for center nodes with:
        # - Multiple red edges at time t1
        # - Multiple blue edges at time t2
        # - (t2 - t1) < stability_threshold

        return reactions

    def find_long_reactions(self, graph: nx.MultiDiGraph, node_info: Dict) -> List:
        """
        Find slow reactions: A + B → I → C + D (through intermediate).

        Returns:
            List of reaction records
        """
        reactions = []

        # TODO: Implement long reaction detection
        # Look for center nodes with:
        # - Multiple red edges at time t1
        # - Multiple blue edges at time t2
        # - (t2 - t1) > stability_threshold

        return reactions

    def find_combinations(self, graph: nx.MultiDiGraph, node_info: Dict) -> List:
        """
        Find combination reactions: A + B → C

        Based on CombAndSplitReactABS() from tool_reactCount.py

        Returns:
            List of reaction records
        """
        reactions = []

        # Look for nodes with multiple red incoming edges at same timestep
        for node, info in node_info.items():
            in_edges = info['in_edges']

            # Group red edges by timestep
            red_groups = {}
            for pred, color, timestep in in_edges:
                if color == 'red':
                    if timestep not in red_groups:
                        red_groups[timestep] = []
                    red_groups[timestep].append(pred)

            # Find combination reactions (multiple reactants at same time)
            for timestep, reactant_nodes in red_groups.items():
                if len(reactant_nodes) > 1:  # A + B → C
                    # Get node attributes
                    reactant_smiles = []
                    reactant_atoms = []
                    for reactant in reactant_nodes:
                        r_data = graph.nodes.get(reactant, {})
                        reactant_smiles.append(r_data.get('smiles', ''))
                        reactant_atoms.append(r_data.get('atom_indices', []))

                    product_data = graph.nodes.get(node, {})

                    reactions.append({
                        'type': 'combination',
                        'reactants': reactant_nodes,
                        'products': [node],
                        'timestep': timestep,
                        'reactant_smiles': reactant_smiles,
                        'product_smiles': [product_data.get('smiles', '')],
                        'reactant_atoms': reactant_atoms,
                        'product_atoms': [product_data.get('atom_indices', [])]
                    })

        logger.debug(f"Found {len(reactions)} combination reactions")
        return reactions

    def find_dissociations(self, graph: nx.MultiDiGraph, node_info: Dict) -> List:
        """
        Find dissociation reactions: C → A + B

        Returns:
            List of reaction records
        """
        reactions = []

        # Look for nodes with multiple blue outgoing edges at same timestep
        for node, info in node_info.items():
            out_edges = info['out_edges']

            # Group blue edges by timestep
            blue_groups = {}
            for succ, color, timestep in out_edges:
                if color == 'blue':
                    if timestep not in blue_groups:
                        blue_groups[timestep] = []
                    blue_groups[timestep].append(succ)

            # Find dissociation reactions (multiple products at same time)
            for timestep, product_nodes in blue_groups.items():
                if len(product_nodes) > 1:  # C → A + B
                    # Get node attributes
                    product_smiles = []
                    product_atoms = []
                    for product in product_nodes:
                        p_data = graph.nodes.get(product, {})
                        product_smiles.append(p_data.get('smiles', ''))
                        product_atoms.append(p_data.get('atom_indices', []))

                    reactant_data = graph.nodes.get(node, {})

                    reactions.append({
                        'type': 'dissociation',
                        'reactants': [node],
                        'products': product_nodes,
                        'timestep': timestep,
                        'reactant_smiles': [reactant_data.get('smiles', '')],
                        'product_smiles': product_smiles,
                        'reactant_atoms': [reactant_data.get('atom_indices', [])],
                        'product_atoms': product_atoms
                    })

        logger.debug(f"Found {len(reactions)} dissociation reactions")
        return reactions

    def find_isomerizations(self, graph: nx.MultiDiGraph) -> List:
        """
        Find isomerization reactions: A → A'

        Based on GreyReactAbs() from tool_reactCount.py

        Returns:
            List of reaction records
        """
        reactions = []

        # Find all gray edges
        for u, v, data in graph.edges(data=True):
            if data.get('color') == 'grey':
                timestep = int(data.get('label', 0))

                # Get node attributes
                u_data = graph.nodes.get(u, {})
                v_data = graph.nodes.get(v, {})

                reactions.append({
                    'type': 'isomerization',
                    'reactants': [u],
                    'products': [v],
                    'timestep': timestep,
                    'reactant_smiles': [u_data.get('smiles', '')],
                    'product_smiles': [v_data.get('smiles', '')],
                    'reactant_atoms': [u_data.get('atom_indices', [])],
                    'product_atoms': [v_data.get('atom_indices', [])]
                })

        logger.debug(f"Found {len(reactions)} isomerization reactions")

        return reactions

    def generate_reaction_hash(
        self,
        reactants: List[str],
        products: List[str],
        catalysts: List[str] = None
    ) -> str:
        """
        Generate unique hash for reaction type.

        Based on getReactionID() from tool_reactCount.py

        Args:
            reactants: List of reactant SMILES
            products: List of product SMILES
            catalysts: List of catalyst SMILES

        Returns:
            SHA1 hash string
        """
        if catalysts is None:
            catalysts = []

        # Sort for consistency
        reactants = sorted(reactants)
        products = sorted(products)
        catalysts = sorted(catalysts)

        # Build hash string
        hash_str = f"R {len(reactants)} "
        hash_str += " ".join(reactants)
        hash_str += f" P {len(products)} "
        hash_str += " ".join(products)
        if catalysts:
            hash_str += " cata "
            hash_str += " ".join(catalysts)

        return hashlib.sha1(hash_str.encode('utf-8')).hexdigest()

    def group_by_hash(self, reactions: List[dict]) -> Dict:
        """
        Group reactions by unique hash.

        Args:
            reactions: List of reaction records

        Returns:
            Dictionary mapping hash to list of events
        """
        grouped = defaultdict(list)

        for reaction in reactions:
            # Generate hash using SMILES, not node IDs
            rxn_hash = self.generate_reaction_hash(
                reaction.get('reactant_smiles', []),
                reaction.get('product_smiles', []),
                reaction.get('catalyst_smiles', [])
            )

            grouped[rxn_hash].append(reaction)

        return dict(grouped)
