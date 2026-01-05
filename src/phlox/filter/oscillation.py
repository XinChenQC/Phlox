"""
Oscillation removal from reaction networks.

Ported from ReaxANA/tool_removeOscill.py
"""

import logging
import networkx as nx
import time
from copy import deepcopy
from typing import Set, Tuple, List

logger = logging.getLogger(__name__)


def sort_list(sub_li: list, index: int) -> list:
    """Sort list by index."""
    sub_li.sort(key=lambda x: x[index])
    return sub_li


def build_basic_info(graph: nx.MultiDiGraph) -> list:
    """Build basic node information."""
    records = []
    for node, data in graph.nodes(data=True):
        node_info = [node, data.get('label', '')]

        # In edges
        in_edges = []
        for edge in graph.in_edges(node, data=True):
            label = edge[2].get('label', 0)
            in_edges.append((edge[0], edge[2].get('color', 'grey'), int(label) if not isinstance(label, int) else label))
        in_edges = sort_list(in_edges, 2)
        node_info.append(in_edges)

        # Out edges
        out_edges = []
        for edge in graph.out_edges(node, data=True):
            label = edge[2].get('label', 0)
            out_edges.append((edge[1], edge[2].get('color', 'grey'), int(label) if not isinstance(label, int) else label))
        out_edges = sort_list(out_edges, 2)
        node_info.append(out_edges)

        # Determine type
        node_type = 'edge'
        for i in range(len(node_info[2]) - 1):
            if (node_info[2][i][2] == node_info[2][i+1][2] and
                node_info[2][i][1] == 'red' and node_info[2][i+1][1] == 'red'):
                node_type = 'center'
                break
        for i in range(len(node_info[3]) - 1):
            if (node_info[3][i][2] == node_info[3][i+1][2] and
                node_info[3][i][1] == 'blue' and node_info[3][i+1][1] == 'blue'):
                node_type = 'center'
                break

        node_info.append(node_type)
        records.append(node_info)

    return records


def rm_single_node(graph: nx.MultiDiGraph):
    """Remove isolated nodes."""
    rm_list = []
    for node in graph.nodes():
        if graph.in_degree(node) == 0 and graph.out_degree(node) == 0:
            rm_list.append(node)
    for node in rm_list:
        graph.remove_node(node)


def rm_edge_list(graph: nx.MultiDiGraph, rm_list: list) -> nx.MultiDiGraph:
    """Remove edges from graph."""
    for itm in rm_list:
        try:
            graph.remove_edge(itm[0], itm[1], key=itm[3])
        except:
            pass
    return graph


def add_key_to_rm_list(graph: nx.MultiDiGraph, rm_list: list) -> list:
    """Add edge keys to remove list."""
    rm_list_upd = []
    for itm in rm_list:
        dicts = graph.get_edge_data(itm[0], itm[1])
        if dicts is None:
            continue
        for key in dicts:
            if dicts[key]['label'] == itm[2]:
                rm_list_upd.append((itm[0], itm[1], itm[2], key))
    return rm_list_upd


class OscillationRemover:
    """
    Remove oscillating species from reaction network.

    Ported from ReaxANA/tool_removeOscill.py
    """

    def __init__(self, stability_threshold: int = 10):
        """Initialize oscillation remover."""
        self.stability_threshold = stability_threshold

    def remove_oscillations(self, graph: nx.MultiDiGraph) -> nx.MultiDiGraph:
        """
        Remove oscillating transformations.

        Ported from remove_useless_trans() in tool_removeOscill.py
        """
        logger.info(f"Removing oscillations (threshold: {self.stability_threshold})")

        graph_copy = deepcopy(graph)
        initial_edges = graph_copy.number_of_edges()

        # Remove all oscillation patterns iteratively
        converged = False
        iteration = 0

        while not converged:
            iteration += 1
            prev_edges = graph_copy.number_of_edges()

            records = build_basic_info(graph_copy)

            # Remove grey oscillations (forward direct)
            rm_list = self._filter_grey_forward_direct(graph_copy, records)
            rm_list_upd = add_key_to_rm_list(graph_copy, rm_list)
            rm_list_upd = list(set(rm_list_upd))
            graph_copy = rm_edge_list(graph_copy, rm_list_upd)

            # Remove grey oscillations (reversed direct)
            records = build_basic_info(graph_copy)
            rm_list = self._filter_grey_reversed_direct(graph_copy, records)
            rm_list_upd = add_key_to_rm_list(graph_copy, rm_list)
            rm_list_upd = list(set(rm_list_upd))
            graph_copy = rm_edge_list(graph_copy, rm_list_upd)

            curr_edges = graph_copy.number_of_edges()
            if curr_edges == prev_edges:
                converged = True

        rm_single_node(graph_copy)

        final_edges = graph_copy.number_of_edges()
        logger.info(f"Removed {initial_edges - final_edges} oscillating edges in {iteration} iterations")

        return graph_copy

    def _filter_grey_forward_direct(self, graph, records):
        """Remove grey edges in forward direct pattern."""
        rm_list = []
        grey_pass = []

        for rec in records:
            if rec[0] in grey_pass:
                continue

            rec_temp = deepcopy(rec)
            # Keep only grey transformations
            for trans in rec[2]:
                if trans[1] in ['red', 'blue']:
                    rec_temp[2].remove(trans)
            for trans in rec[3]:
                if trans[1] in ['red', 'blue']:
                    rec_temp[3].remove(trans)

            if len(rec_temp[2]) == 0 or len(rec_temp[3]) == 0:
                continue

            for i in range(min(len(rec_temp[2]), len(rec_temp[3]))):
                if (rec_temp[2][i][2] - rec_temp[3][i][2] < self.stability_threshold and
                    rec_temp[2][i][2] > rec_temp[3][i][2] and
                    rec_temp[3][i][0] == rec_temp[2][i][0]):
                    rm_list.append((rec_temp[2][i][0], rec_temp[0], rec_temp[2][i][2]))
                    rm_list.append((rec_temp[0], rec_temp[3][i][0], rec_temp[3][i][2]))
                    grey_pass.append(rec_temp[2][i][0])

        return rm_list

    def _filter_grey_reversed_direct(self, graph, records):
        """Remove grey edges in reversed direct pattern."""
        rm_list = []
        grey_pass = []

        for rec in records:
            if rec[0] in grey_pass:
                continue

            rec_temp = deepcopy(rec)
            # Keep only grey transformations
            for trans in rec[2]:
                if trans[1] in ['red', 'blue']:
                    rec_temp[2].remove(trans)
            for trans in rec[3]:
                if trans[1] in ['red', 'blue']:
                    rec_temp[3].remove(trans)

            if len(rec_temp[2]) == 0 or len(rec_temp[3]) == 0:
                continue

            for i in range(min(len(rec_temp[2]), len(rec_temp[3]))):
                if (rec_temp[3][i][2] - rec_temp[2][i][2] < self.stability_threshold and
                    rec_temp[3][i][2] > rec_temp[2][i][2] and
                    rec_temp[3][i][0] == rec_temp[2][i][0]):
                    rm_list.append((rec_temp[2][i][0], rec_temp[0], rec_temp[2][i][2]))
                    rm_list.append((rec_temp[0], rec_temp[3][i][0], rec_temp[3][i][2]))
                    grey_pass.append(rec_temp[2][i][0])

        return rm_list
