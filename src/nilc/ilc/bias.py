"""
Bias mitigation strategies for ILC component separation.

Implements harmonic-space and real-space bias exclusion methods from
Section III-C of arXiv:2307.01258.
"""
import numpy as np
from pixell import enmap, fft


def apply_harmonic_exclusion(needlet_coeffs, ell_min_exclude, ell_max_exclude, wcs):
    """
    Apply harmonic-space bias exclusion for large-scale modes.

    Excludes certain ell ranges from the covariance computation to prevent
    ILC bias. This is used for large angular scales where the harmonic
    exclusion strategy is effective.

    Parameters
    ----------
    needlet_coeffs : ndarray
        Needlet coefficients with shape (n_freq, n_y, n_x)
    ell_min_exclude : float
        Minimum ell to exclude
    ell_max_exclude : float
        Maximum ell to exclude
    wcs : astropy.wcs.WCS
        WCS information for map geometry

    Returns
    -------
    coeffs_excluded : ndarray
        Modified needlet coefficients with excluded modes removed

    Notes
    -----
    For computing covariance on large scales, we exclude the modes that will
    be used in the final ILC application to reduce bias from signal-noise
    correlation.
    """
    n_freq = needlet_coeffs.shape[0]
    shape = needlet_coeffs.shape[1:3]

    # Get modlmap
    modlmap = enmap.modlmap(shape, wcs)

    # Create exclusion mask in Fourier space
    exclusion_mask = np.ones_like(modlmap)
    exclusion_mask[(modlmap >= ell_min_exclude) & (modlmap <= ell_max_exclude)] = 0.0

    # Apply exclusion to each frequency
    coeffs_excluded = np.zeros_like(needlet_coeffs)
    for ifreq in range(n_freq):
        # Transform to Fourier space
        coeff_map = enmap.ndmap(needlet_coeffs[ifreq], wcs)
        fmap = fft.fft(coeff_map, axes=[-2, -1])

        # Apply exclusion mask
        fmap_excluded = fmap * exclusion_mask

        # Transform back
        coeff_excluded = fft.ifft(fmap_excluded, axes=[-2, -1]).real
        coeffs_excluded[ifreq] = coeff_excluded

    return coeffs_excluded


def compute_bias_correction(weights, needlet_coeffs, cov_map):
    """
    Estimate and correct for ILC bias.

    The ILC can be biased when the same data is used to compute weights
    and apply them. This function estimates the bias and can be used
    to apply a correction.

    Parameters
    ----------
    weights : ndarray
        ILC weights with shape (n_y, n_x, n_freq)
    needlet_coeffs : ndarray
        Needlet coefficients with shape (n_freq, n_y, n_x)
    cov_map : ndarray
        Covariance matrix with shape (n_y, n_x, n_freq, n_freq)

    Returns
    -------
    bias_estimate : ndarray
        Estimated bias with shape (n_y, n_x)

    Notes
    -----
    This is a simplified bias estimator. The paper uses more sophisticated
    strategies (harmonic exclusion for large scales, donut smoothing for
    small scales) to avoid bias in the first place.
    """
    n_y, n_x, n_freq = weights.shape

    # Estimate bias as w^T @ C @ w (the variance of the ILC map)
    bias_estimate = np.zeros((n_y, n_x))

    for iy in range(n_y):
        for ix in range(n_x):
            w = weights[iy, ix, :]  # (n_freq,)
            C = cov_map[iy, ix, :, :]  # (n_freq, n_freq)

            # Variance = w^T @ C @ w
            variance = w @ C @ w
            bias_estimate[iy, ix] = variance

    return bias_estimate


def select_bias_method(ell_peak):
    """
    Select appropriate bias mitigation method based on needlet scale.

    Following Section III-C of the paper:
    - Large scales (small ell_peak): use harmonic-space exclusion
    - Small scales (large ell_peak): use real-space donut smoothing

    Parameters
    ----------
    ell_peak : float
        Peak ell value for the needlet scale

    Returns
    -------
    method : str
        'harmonic' for harmonic-space exclusion
        'donut' for real-space donut smoothing
    threshold_ell : float
        Transition ell between methods

    Notes
    -----
    The paper doesn't specify an exact transition point. Based on typical
    CMB analysis, we use ell ~ 2000-3000 as a reasonable transition.
    """
    # Transition threshold (adjust based on testing)
    transition_ell = 2500

    if ell_peak < transition_ell:
        return 'harmonic'
    else:
        return 'donut'


def get_exclusion_parameters(ell_peak, method='harmonic'):
    """
    Get appropriate exclusion parameters for a given needlet scale.

    Parameters
    ----------
    ell_peak : float
        Peak ell for the needlet scale
    method : str
        'harmonic' or 'donut'

    Returns
    -------
    params : dict
        Dictionary of parameters for the exclusion method

    Notes
    -----
    These are reasonable defaults that should be tuned based on the
    specific needlet kernel and desired bias/variance trade-off.
    """
    if method == 'harmonic':
        # For harmonic exclusion: exclude a band around ell_peak
        # Exclude roughly +/- half the needlet width
        delta_ell = ell_peak * 0.5  # Adjust based on needlet width
        params = {
            'ell_min': max(0, ell_peak - delta_ell),
            'ell_max': ell_peak + delta_ell
        }
    else:  # donut
        # For donut smoothing: define inner and outer radii in arcminutes
        # Inner radius ~ wavelength at ell_peak
        # Outer radius ~ 2-3x inner radius
        theta_peak_arcmin = 10800.0 / (np.pi * ell_peak)  # Convert ell to arcmin
        params = {
            'inner_radius': theta_peak_arcmin,
            'outer_radius': theta_peak_arcmin * 2.5
        }

    return params
