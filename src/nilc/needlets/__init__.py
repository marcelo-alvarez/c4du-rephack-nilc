"""Needlet decomposition module for NILC component separation.

This module implements the needlet frame as described in Section III-B of
arXiv:2307.01258, using axisymmetric needlet kernels for multi-scale
decomposition of CAR maps.
"""

from .decomposition import (
    compute_needlet_kernel,
    compute_needlet_kernels,
    needlet_decompose,
    needlet_recombine,
    ELL_PEAKS,
)

__all__ = [
    'compute_needlet_kernel',
    'compute_needlet_kernels',
    'needlet_decompose',
    'needlet_recombine',
    'ELL_PEAKS',
]
