"""
Component separation using Internal Linear Combination (ILC).

Main module that orchestrates the NILC pipeline for separating CMB and kSZ
components from multi-frequency needlet coefficients.
"""
import numpy as np
from . import covariance, weights, bias


def apply_ilc_weights(needlet_coeffs, ilc_weights):
    """
    Apply ILC weights to needlet coefficients to extract component.

    Parameters
    ----------
    needlet_coeffs : ndarray
        Needlet coefficients with shape (n_freq, n_y, n_x)
    ilc_weights : ndarray
        ILC weights with shape (n_y, n_x, n_freq)

    Returns
    -------
    component_map : ndarray
        Separated component with shape (n_y, n_x)

    Notes
    -----
    The ILC combination is: s = sum_i w_i * d_i
    where w_i are the weights and d_i are the data at each frequency.
    """
    n_freq, n_y, n_x = needlet_coeffs.shape

    # Initialize output
    component_map = np.zeros((n_y, n_x))

    # Apply weights
    for ifreq in range(n_freq):
        component_map += ilc_weights[:, :, ifreq] * needlet_coeffs[ifreq, :, :]

    return component_map


def separate_components_single_scale(needlet_coeffs, ell_peak, wcs,
                                     component_responses, mask=None,
                                     use_bias_mitigation=True):
    """
    Perform ILC component separation for a single needlet scale.

    Parameters
    ----------
    needlet_coeffs : ndarray
        Needlet coefficients with shape (n_freq, n_y, n_x)
    ell_peak : float
        Peak ell for this needlet scale
    wcs : astropy.wcs.WCS
        WCS information for map geometry
    component_responses : ndarray
        Response matrix with shape (n_freq, n_comp)
        Each column is a component's spectral signature
    mask : enmap.ndmap, optional
        Analysis mask
    use_bias_mitigation : bool
        Whether to apply bias mitigation strategies

    Returns
    -------
    components : ndarray
        Separated components with shape (n_comp, n_y, n_x)
    ilc_weights : ndarray
        ILC weights used, shape (n_y, n_x, n_freq, n_comp)

    Notes
    -----
    This function implements the full NILC procedure for one scale:
    1. Optionally apply bias mitigation (harmonic exclusion or donut)
    2. Compute covariance matrix
    3. Compute ILC weights
    4. Apply weights to separate components
    """
    n_freq, n_y, n_x = needlet_coeffs.shape
    n_comp = component_responses.shape[1]

    # Determine bias mitigation method
    if use_bias_mitigation:
        bias_method = bias.select_bias_method(ell_peak)
        bias_params = bias.get_exclusion_parameters(ell_peak, bias_method)
    else:
        bias_method = None

    # Compute covariance with appropriate bias mitigation
    if bias_method == 'harmonic':
        # Apply harmonic exclusion
        coeffs_for_cov = bias.apply_harmonic_exclusion(
            needlet_coeffs,
            bias_params['ell_min'],
            bias_params['ell_max'],
            wcs
        )
        # Compute covariance from excluded coefficients
        cov_map = covariance.compute_covariance_map(coeffs_for_cov, mask=mask)

    elif bias_method == 'donut':
        # Use donut smoothing for covariance
        cov_map = covariance.compute_covariance_donut(
            needlet_coeffs,
            bias_params['inner_radius'],
            bias_params['outer_radius'],
            wcs,
            mask=mask
        )

    else:
        # No bias mitigation
        cov_map = covariance.compute_covariance_map(needlet_coeffs, mask=mask)

    # Invert covariance
    cov_inv = covariance.invert_covariance(cov_map, regularization=1e-10)

    # Compute ILC weights for all components
    ilc_weights = weights.compute_constrained_ilc_weights(cov_inv, component_responses)

    # Separate components
    components = np.zeros((n_comp, n_y, n_x))
    for icomp in range(n_comp):
        # Extract weights for this component
        w_comp = ilc_weights[:, :, :, icomp]  # (n_y, n_x, n_freq)

        # Apply weights to original (non-excluded) needlet coefficients
        components[icomp] = apply_ilc_weights(needlet_coeffs, w_comp)

    return components, ilc_weights


def separate_components_multiScale(needlet_coeffs_all_scales, ell_peaks, wcs,
                                   component_responses, mask=None,
                                   use_bias_mitigation=True):
    """
    Perform ILC component separation across all needlet scales.

    Parameters
    ----------
    needlet_coeffs_all_scales : list of ndarray
        List of needlet coefficients, one per scale
        Each element has shape (n_freq, n_y, n_x)
    ell_peaks : array_like
        Peak ell values for each needlet scale
    wcs : astropy.wcs.WCS
        WCS information
    component_responses : ndarray
        Response matrix with shape (n_freq, n_comp)
    mask : enmap.ndmap, optional
        Analysis mask
    use_bias_mitigation : bool
        Whether to apply bias mitigation

    Returns
    -------
    components_all_scales : list of ndarray
        Separated components for each scale
        Each element has shape (n_comp, n_y, n_x)
    weights_all_scales : list of ndarray
        ILC weights for each scale

    Notes
    -----
    This applies the ILC procedure independently to each needlet scale.
    The final maps are obtained by summing across scales (done in synthesis).
    """
    n_scales = len(needlet_coeffs_all_scales)
    components_all_scales = []
    weights_all_scales = []

    for iscale in range(n_scales):
        print(f"Processing scale {iscale+1}/{n_scales}, ell_peak={ell_peaks[iscale]}")

        components, ilc_weights = separate_components_single_scale(
            needlet_coeffs_all_scales[iscale],
            ell_peaks[iscale],
            wcs,
            component_responses,
            mask=mask,
            use_bias_mitigation=use_bias_mitigation
        )

        components_all_scales.append(components)
        weights_all_scales.append(ilc_weights)

    return components_all_scales, weights_all_scales


def synthesize_components(components_all_scales):
    """
    Synthesize final component maps from multi-scale needlet coefficients.

    Parameters
    ----------
    components_all_scales : list of ndarray
        Separated components for each scale
        Each element has shape (n_comp, n_y, n_x)

    Returns
    -------
    final_components : ndarray
        Final synthesized component maps with shape (n_comp, n_y, n_x)

    Notes
    -----
    The synthesis is simply a sum across needlet scales, since the
    needlet decomposition is designed such that the original map is
    the sum of all needlet bands.
    """
    # Sum across scales
    final_components = sum(components_all_scales)

    return final_components


def nilc_pipeline(needlet_coeffs_all_scales, ell_peaks, wcs,
                 frequencies_ghz, mask=None, components_to_separate=('CMB', 'kSZ')):
    """
    Full NILC pipeline for component separation.

    Parameters
    ----------
    needlet_coeffs_all_scales : list of ndarray
        Needlet coefficients for all scales
    ell_peaks : array_like
        Peak ell values for needlet scales
    wcs : astropy.wcs.WCS
        WCS information
    frequencies_ghz : array_like
        Observing frequencies in GHz (e.g., [90, 150] for ACT)
    mask : enmap.ndmap, optional
        Analysis mask
    components_to_separate : tuple of str
        Components to separate, e.g., ('CMB', 'kSZ')

    Returns
    -------
    separated_maps : dict
        Dictionary mapping component name to separated map
        Each map has shape (n_y, n_x)

    Examples
    --------
    >>> separated = nilc_pipeline(needlet_coeffs, ell_peaks, wcs,
    ...                          frequencies_ghz=[90, 150],
    ...                          components_to_separate=('CMB', 'kSZ'))
    >>> cmb_map = separated['CMB']
    >>> ksz_map = separated['kSZ']
    """
    # Build response matrix
    n_freq = len(frequencies_ghz)
    n_comp = len(components_to_separate)
    response_matrix = np.zeros((n_freq, n_comp))

    for icomp, comp_name in enumerate(components_to_separate):
        if comp_name.upper() == 'CMB':
            response_matrix[:, icomp] = weights.get_cmb_response(frequencies_ghz)
        elif comp_name.upper() == 'KSZ':
            response_matrix[:, icomp] = weights.get_ksz_response(frequencies_ghz)
        elif comp_name.upper() == 'TSZ':
            response_matrix[:, icomp] = weights.get_tsz_response(frequencies_ghz)
        else:
            raise ValueError(f"Unknown component: {comp_name}")

    # Perform multi-scale ILC
    components_all_scales, _ = separate_components_multiScale(
        needlet_coeffs_all_scales,
        ell_peaks,
        wcs,
        response_matrix,
        mask=mask,
        use_bias_mitigation=True
    )

    # Synthesize final maps
    final_components = synthesize_components(components_all_scales)

    # Package as dictionary
    separated_maps = {}
    for icomp, comp_name in enumerate(components_to_separate):
        separated_maps[comp_name] = final_components[icomp]

    return separated_maps
