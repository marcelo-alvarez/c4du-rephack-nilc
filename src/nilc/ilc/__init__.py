"""Internal Linear Combination (ILC) component separation with bias mitigation.

This module implements the ILC method for separating CMB and kSZ components
from multi-frequency microwave sky maps, following the approach described in
Section III-C of arXiv:2307.01258 (The Atacama Cosmology Telescope:
High-resolution component-separated maps).

Key Features
------------
- In-band covariance computation for each needlet scale
- ILC weight calculation with standard minimum-variance formulation
- Hybrid bias mitigation strategy:
  * Large scales (ell < ~1000): Harmonic-space mode exclusion
  * Small scales (ell >= ~1000): Real-space "donut" smoothing
- Support for multi-component separation (CMB, kSZ, etc.)

Main Functions
--------------
separate_components : Main entry point for component separation
apply_ilc_weights : Apply computed weights to extract components
reconstruct_full_component : Recombine needlet scales to full map

Modules
-------
covariance : Covariance matrix computation
weights : ILC weight calculation
bias_mitigation : Harmonic exclusion and donut smoothing
separation : Main component separation logic

Example
-------
>>> from nilc.ilc import separate_components
>>> # Assuming needlet-decomposed maps are available
>>> needlet_maps = {100: [map_90_ell100, map_150_ell100], ...}
>>> separated = separate_components(needlet_maps, mask, ell_peaks)
>>> cmb_needlets = separated['cmb']
>>> ksz_needlets = separated['ksz']
"""

from .separation import separate_components, apply_ilc_weights, reconstruct_full_component
from .covariance import compute_inband_covariance
from .weights import compute_ilc_weights
from .bias_mitigation import apply_bias_mitigation

__all__ = [
    'separate_components',
    'apply_ilc_weights',
    'reconstruct_full_component',
    'compute_inband_covariance',
    'compute_ilc_weights',
    'apply_bias_mitigation',
]
