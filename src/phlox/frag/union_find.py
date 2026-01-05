"""
Union-Find data structure for molecular fragmentation.

Implements weighted quick-union with path compression for efficient
disjoint set operations.
"""

from typing import List, Set
import numpy as np


class UnionFind:
    """
    Union-Find disjoint set data structure.

    Maintains connected components (molecular fragments) and supports
    efficient union and find operations.

    Time complexity: O(α(N)) amortized for union/find operations
    where α is the inverse Ackermann function (effectively constant).
    """

    def __init__(self, n_elements: int):
        """
        Initialize union-find structure.

        Args:
            n_elements: Number of elements (atoms)
        """
        self.n_elements = n_elements
        self.n_components = n_elements
        self.parent = list(range(n_elements))  # Each element is its own parent initially
        self.size = [1] * n_elements  # Size of each component

    @classmethod
    def from_connectivity_matrix(cls, connectivity: np.ndarray) -> 'UnionFind':
        """
        Create UnionFind from boolean connectivity matrix.

        Args:
            connectivity: N x N boolean array where connectivity[i,j] = True
                         means atoms i and j are bonded

        Returns:
            UnionFind instance with connected components
        """
        n = len(connectivity)
        uf = cls(n)

        # Find all connections and union them
        rows, cols = np.nonzero(connectivity)
        for i, j in zip(rows, cols):
            if i < j:  # Avoid double counting
                uf.union(i, j)

        return uf

    def find(self, x: int) -> int:
        """
        Find root of element's component with path compression.

        Args:
            x: Element index

        Returns:
            Root index of component containing x
        """
        if x < 0 or x >= self.n_elements:
            raise ValueError(f"Element {x} out of range [0, {self.n_elements})")

        # Path compression: make all nodes point directly to root
        root = x
        while root != self.parent[root]:
            root = self.parent[root]

        # Compress path
        while x != root:
            next_x = self.parent[x]
            self.parent[x] = root
            x = next_x

        return root

    def union(self, x: int, y: int):
        """
        Merge components containing x and y using weighted union.

        Args:
            x: First element
            y: Second element
        """
        root_x = self.find(x)
        root_y = self.find(y)

        if root_x == root_y:
            return  # Already in same component

        # Weighted union: attach smaller tree to larger tree
        if self.size[root_x] < self.size[root_y]:
            self.parent[root_x] = root_y
            self.size[root_y] += self.size[root_x]
        else:
            self.parent[root_y] = root_x
            self.size[root_x] += self.size[root_y]

        self.n_components -= 1

    def connected(self, x: int, y: int) -> bool:
        """
        Check if two elements are in the same component.

        Args:
            x: First element
            y: Second element

        Returns:
            True if x and y are connected
        """
        return self.find(x) == self.find(y)

    def get_components(self) -> List[Set[int]]:
        """
        Get all connected components.

        Returns:
            List of sets, where each set contains indices of atoms in a molecule
        """
        components = {}
        for i in range(self.n_elements):
            root = self.find(i)
            if root not in components:
                components[root] = set()
            components[root].add(i)

        return list(components.values())

    def component_size(self, x: int) -> int:
        """
        Get size of component containing x.

        Args:
            x: Element index

        Returns:
            Number of elements in component
        """
        root = self.find(x)
        return self.size[root]
