"""Map I/O utilities for ACT data."""
import numpy as np
from pixell import enmap

def load_act_map(filepath):
    """
    Load an ACT map from FITS file.

    Parameters
    ----------
    filepath : str
        Path to FITS file

    Returns
    -------
    map : enmap.ndmap
        Loaded map with WCS information
    """
    return enmap.read_map(filepath)

def get_stokes_component(map_data, component='I'):
    """
    Extract a specific Stokes component from map.

    Parameters
    ----------
    map_data : enmap.ndmap
        Full map with potentially multiple components
    component : str
        'I' (temperature), 'Q', or 'U' (polarization)

    Returns
    -------
    component_map : enmap.ndmap
        Single component map
    """
    component_idx = {'I': 0, 'Q': 1, 'U': 2}
    if map_data.ndim == 3:
        return map_data[component_idx[component]]
    return map_data
