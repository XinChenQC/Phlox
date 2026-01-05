"""
Phlox - Reaction mechanism discovery from reactive molecular dynamics trajectories.

This package provides tools for analyzing reactive MD simulations to discover
reaction mechanisms, track molecular species, and build reaction networks.
"""

__version__ = "0.1.0"
__author__ = "Chen Xin"
__email__ = "chenxin199261@gmail.com"

from phlox.base.config import Config
from phlox.base.analyzer import Analyzer
from phlox.base.results import Results

__all__ = ["Config", "Analyzer", "Results", "__version__"]
