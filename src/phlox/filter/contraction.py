"""
Node contraction for reaction networks.

Ported from ReaxANA/tool_contract.py
"""

import logging
import networkx as nx
from copy import deepcopy
from typing import List, Tuple, Set

logger = logging.getLogger(__name__)


class NodeContractor:
    """
    Contract redundant nodes in reaction network.

    Merges nodes that represent the same species or intermediate states.

    Based on ReaxANA/tool_contract.py
    """

    def __init__(self, stability_threshold: int = 10):
        """Initialize node contractor."""
        self.stability_threshold = stability_threshold

    def contract_nodes(self, graph: nx.MultiDiGraph) -> nx.MultiDiGraph:
        """
        Contract redundant nodes in the graph.

        Merges nodes that:
        - Have identical SMILES
        - Represent the same intermediate state
        - Can be simplified without losing information

        Based on contract_nodes() from tool_contract.py

        Args:
            graph: Reaction network

        Returns:
            Contracted graph
        """
        logger.info("Contracting redundant nodes")

        initial_nodes = graph.number_of_nodes()
        initial_edges = graph.number_of_edges()

        # Ported from tool_contract.py
        graph_copy = deepcopy(graph)
        pass_pair = []
        converged = False

        while not converged:
            prev_nodes = graph_copy.number_of_nodes()
            contr_pairs, curr_pass = self._contract_grey(graph_copy, pass_pair)

            for cont in contr_pairs:
                try:
                    graph_copy = nx.contracted_nodes(graph_copy, cont[1], cont[0], self_loops=False, copy=False)
                except:
                    pass

            for itm in curr_pass:
                pass_pair.append(itm)

            curr_nodes = graph_copy.number_of_nodes()
            if curr_nodes == prev_nodes:
                converged = True

        final_nodes = graph.number_of_nodes()
        final_edges = graph.number_of_edges()

        logger.info(f"Contracted {initial_nodes - graph_copy.number_of_nodes()} nodes")

        return graph_copy

    def _contract_grey(self, graph, pass_pair):
        """Contract grey transformations."""
        grey_list = []

        for u, v, data in graph.edges(data=True):
            if (data.get('color') == 'grey' and
                (u, v) not in grey_list and (v, u) not in grey_list and
                (u, v) not in pass_pair and (v, u) not in pass_pair and u != v):
                grey_list.append((u, v))

        # Build node lifetime
        built_rec = {}
        for pair in grey_list:
            for node in pair:
                if node in built_rec:
                    continue

                min_time = 999999
                max_time = 0

                for _, _, data in graph.in_edges(node, data=True):
                    min_time = min(min_time, int(data.get('label', 0)))
                for _, _, data in graph.out_edges(node, data=True):
                    max_time = max(max_time, int(data.get('label', 0)))

                built_rec[node] = max_time - min_time if max_time > min_time else 0

        # Contract pairs
        contr_result = []
        curr_pass_pair = []

        for pair in grey_list:
            lifetime_a = built_rec.get(pair[0], 999)
            lifetime_b = built_rec.get(pair[1], 999)

            if lifetime_a < lifetime_b and lifetime_a < self.stability_threshold:
                contr_result.append((pair[0], pair[1]))
            elif lifetime_b < lifetime_a and lifetime_b < self.stability_threshold:
                contr_result.append((pair[1], pair[0]))
            else:
                curr_pass_pair.append(pair)

        return contr_result, curr_pass_pair

    def contract_gray_edges(self, graph: nx.MultiDiGraph) -> nx.MultiDiGraph:
        """
        Contract nodes connected only by gray (isomerization) edges.

        Gray edges represent A → A' transformations that may be artifacts.

        Based on contract_nodes_Onlygrey() from tool_contract.py

        Args:
            graph: Reaction network

        Returns:
            Contracted graph
        """
        logger.info("Contracting gray edge chains")

        # TODO: Implement gray edge contraction
        # Algorithm:
        # 1. Find all gray edges (color='grey')
        # 2. Find chains of gray edges: A → B → C (all gray)
        # 3. Contract chain into single node
        # 4. Update edge connections

        # Placeholder
        logger.warning("Gray edge contraction not yet implemented")

        return graph

    def find_contractible_pairs(
        self,
        graph: nx.MultiDiGraph
    ) -> List[Tuple[str, str]]:
        """
        Find pairs of nodes that can be contracted.

        Args:
            graph: Reaction network

        Returns:
            List of (node1, node2) pairs that can be merged
        """
        contractible = []

        # TODO: Implement contractible pair detection
        # Criteria:
        # - Same molecular formula/SMILES
        # - Simple path between them
        # - No branching

        return contractible

    def merge_nodes(
        self,
        graph: nx.MultiDiGraph,
        node1: str,
        node2: str,
        keep: str = None
    ) -> nx.MultiDiGraph:
        """
        Merge two nodes into one.

        Args:
            graph: Reaction network
            node1: First node to merge
            node2: Second node to merge
            keep: Which node to keep (node1, node2, or None for new node)

        Returns:
            Updated graph
        """
        if keep is None:
            keep = node1

        # TODO: Implement node merging
        # 1. Redirect all edges from node2 to keep node
        # 2. Copy attributes
        # 3. Remove node2

        return graph
