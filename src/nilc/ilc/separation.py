"""Main ILC component separation module."""
import numpy as np

try:
    from pixell import enmap
    PIXELL_AVAILABLE = True
except ImportError:
    PIXELL_AVAILABLE = False
    # Mock enmap for testing without pixell
    class enmap:
        @staticmethod
        def zeros_like(template):
            return np.zeros_like(template)

from .covariance import compute_inband_covariance
from .weights import compute_ilc_weights
from .bias_mitigation import apply_bias_mitigation


def separate_components(needlet_maps, mask, ell_peaks, components=['cmb', 'ksz'],
                        transition_ell=1000, exclusion_width=0.3,
                        smooth_scale_factor=2.0):
    """
    Separate components using ILC with bias mitigation.

    This is the main entry point for ILC component separation. For each needlet
    scale, it applies appropriate bias mitigation, computes covariance, calculates
    ILC weights, and extracts the requested components.

    Parameters
    ----------
    needlet_maps : dict
        Dictionary mapping ell_peak to list of frequency maps at that scale.
        Structure: {ell_peak: [map_90GHz, map_150GHz, ...]}
        Each map is an enmap.ndmap at a specific needlet scale.
    mask : enmap.ndmap
        Analysis mask (1 = include, 0 = exclude)
    ell_peaks : list or array
        List of needlet peak multipoles (e.g., [100, 200, 300, ...])
    components : list of str
        Components to separate. Options: 'cmb', 'ksz', 'tsz'
        Default: ['cmb', 'ksz']
    transition_ell : float
        Transition multipole between harmonic and donut bias mitigation
        Default: 1000
    exclusion_width : float
        Width for harmonic exclusion (fraction of ell_peak)
    smooth_scale_factor : float
        Smoothing factor for donut method

    Returns
    -------
    separated_components : dict
        Dictionary with separated components and needlet-scale outputs
        Structure: {
            'cmb': {ell_peak: component_map, ...},
            'ksz': {ell_peak: component_map, ...},
            'weights': {
                'cmb': {ell_peak: weights_array, ...},
                'ksz': {ell_peak: weights_array, ...}
            },
            'metadata': {
                'ell_peaks': ell_peaks,
                'strategies': {ell_peak: 'harmonic' or 'donut', ...}
            }
        }

    Notes
    -----
    Following Section III-C of arXiv:2307.01258:
    1. For each needlet scale:
       a. Apply appropriate bias mitigation (harmonic or donut)
       b. Compute in-band covariance matrix
       c. Compute ILC weights for each component
       d. Apply weights to extract component at that scale
    2. Large scales (ell < transition_ell) use harmonic-space exclusion
    3. Small scales (ell >= transition_ell) use real-space donut smoothing

    For ACT with two frequencies (90 GHz, 150 GHz):
    - CMB SED: [1, 1] (achromatic in thermodynamic units)
    - kSZ SED: [1, 1] (also achromatic)
    - Both can be separated even with same SED due to different noise properties

    Examples
    --------
    >>> # Assuming you have needlet-decomposed maps
    >>> needlet_maps = {
    ...     100: [map_90_scale100, map_150_scale100],
    ...     200: [map_90_scale200, map_150_scale200],
    ... }
    >>> ell_peaks = [100, 200]
    >>> separated = separate_components(needlet_maps, mask, ell_peaks)
    >>> cmb_map_scale100 = separated['cmb'][100]
    >>> ksz_map_scale200 = separated['ksz'][200]
    """
    # Initialize output structure
    separated_components = {comp: {} for comp in components}
    separated_components['weights'] = {comp: {} for comp in components}
    separated_components['metadata'] = {
        'ell_peaks': ell_peaks,
        'strategies': {}
    }

    # Define component SEDs
    # For two frequencies (90 GHz and 150 GHz):
    component_seds = {
        'cmb': np.array([1.0, 1.0]),  # Achromatic in thermodynamic units
        'ksz': np.array([1.0, 1.0]),  # Achromatic
        'tsz': None  # Would need frequency-dependent SED
    }

    # Check that requested components have defined SEDs
    for comp in components:
        if comp not in component_seds:
            raise ValueError(f"Unknown component: {comp}")
        if component_seds[comp] is None:
            raise NotImplementedError(f"SED for {comp} not implemented")

    # Process each needlet scale
    for ell_peak in ell_peaks:
        if ell_peak not in needlet_maps:
            raise ValueError(f"No needlet maps provided for ell_peak={ell_peak}")

        maps_at_scale = needlet_maps[ell_peak]
        n_freq = len(maps_at_scale)

        # Verify SED dimensions match number of frequencies
        for comp in components:
            sed = component_seds[comp]
            if len(sed) != n_freq:
                raise ValueError(
                    f"SED for {comp} has length {len(sed)}, "
                    f"but {n_freq} frequency maps provided"
                )

        # Apply bias mitigation
        processed_maps, original_maps, strategy = apply_bias_mitigation(
            maps_at_scale, ell_peak, mask, transition_ell,
            exclusion_width, smooth_scale_factor
        )
        separated_components['metadata']['strategies'][ell_peak] = strategy

        # Compute covariance from processed maps
        cov = compute_inband_covariance(processed_maps, mask)

        # For each component, compute weights and extract
        for comp in components:
            sed = component_seds[comp]

            # Compute ILC weights
            weights = compute_ilc_weights(cov, sed)
            separated_components['weights'][comp][ell_peak] = weights

            # Apply weights to original maps to extract component
            # For donut method: apply to unsmoothed maps
            # For harmonic method: original_maps = processed_maps
            component_map = apply_ilc_weights(original_maps, weights, mask)

            separated_components[comp][ell_peak] = component_map

    return separated_components


def apply_ilc_weights(frequency_maps, weights, mask=None):
    """
    Apply ILC weights to frequency maps to extract a component.

    Parameters
    ----------
    frequency_maps : list of enmap.ndmap
        Maps at different frequencies
    weights : ndarray
        ILC weights, one per frequency
    mask : enmap.ndmap, optional
        Analysis mask to apply to output

    Returns
    -------
    component_map : enmap.ndmap
        Extracted component map

    Notes
    -----
    The component map is the weighted sum: sum_i(w_i * map_i)
    where w_i are the ILC weights and map_i are the frequency maps.
    """
    n_freq = len(frequency_maps)
    if len(weights) != n_freq:
        raise ValueError(
            f"Number of weights ({len(weights)}) must match "
            f"number of frequency maps ({n_freq})"
        )

    # Initialize output with same shape and WCS as input
    component_map = enmap.zeros_like(frequency_maps[0])

    # Weighted sum
    for i, (map_i, w_i) in enumerate(zip(frequency_maps, weights)):
        component_map += w_i * map_i

    # Apply mask if provided
    if mask is not None:
        component_map *= mask

    return component_map


def reconstruct_full_component(component_needlet_maps, ell_peaks):
    """
    Reconstruct full-resolution component map from needlet scales.

    This combines the separated component maps across all needlet scales
    to produce the final full-resolution component map.

    Parameters
    ----------
    component_needlet_maps : dict
        Dictionary mapping ell_peak to component map at that scale
        Structure: {ell_peak: component_map, ...}
    ell_peaks : list
        List of needlet peak multipoles

    Returns
    -------
    full_map : enmap.ndmap
        Full-resolution reconstructed component map

    Notes
    -----
    The needlet decomposition is nearly orthogonal, so reconstruction is
    simply the sum of all needlet scales:
        component_full = sum_j(component_j)
    where j indexes needlet scales.
    """
    # Get reference map for shape and WCS
    first_ell = ell_peaks[0]
    reference_map = component_needlet_maps[first_ell]

    # Initialize output
    full_map = enmap.zeros_like(reference_map)

    # Sum over all needlet scales
    for ell_peak in ell_peaks:
        if ell_peak in component_needlet_maps:
            full_map += component_needlet_maps[ell_peak]
        else:
            raise ValueError(f"Missing needlet scale ell_peak={ell_peak}")

    return full_map
