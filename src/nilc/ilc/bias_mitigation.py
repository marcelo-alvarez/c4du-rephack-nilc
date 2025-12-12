"""Bias mitigation strategies for ILC component separation."""
import numpy as np

try:
    from pixell import enmap, fft
    PIXELL_AVAILABLE = True
except ImportError:
    PIXELL_AVAILABLE = False
    enmap = None
    fft = None


def apply_harmonic_exclusion(needlet_maps, ell_peak, exclusion_width=0.3):
    """
    Apply harmonic-space exclusion for large-scale bias mitigation.

    This method is used for large scales (low ell_peak) to avoid ILC bias.
    It excludes certain harmonic modes from the covariance computation
    to prevent signal-signal correlation bias.

    Parameters
    ----------
    needlet_maps : list of enmap.ndmap
        Needlet-decomposed maps for each frequency at a specific scale
    ell_peak : float
        Peak multipole of the needlet scale
    exclusion_width : float
        Fractional width of exclusion region around ell_peak
        Default 0.3 means exclude ell in [ell_peak*(1-0.3), ell_peak*(1+0.3)]

    Returns
    -------
    filtered_maps : list of enmap.ndmap
        Maps with excluded harmonic modes set to zero

    Notes
    -----
    Following Section III-C of arXiv:2307.01258, large scales use harmonic
    exclusion to compute weights from modes that don't overlap with the
    signal modes being cleaned.

    The exclusion prevents the ILC from "overfitting" to chance correlations
    between signal and noise at the same modes.
    """
    if not PIXELL_AVAILABLE:
        raise ImportError("pixell is required for harmonic exclusion")

    filtered_maps = []

    for map_i in needlet_maps:
        # Transform to Fourier space
        fmap = fft.fft(map_i, axes=[-2, -1])

        # Create modlmap (2D ell values)
        modlmap = enmap.modlmap(map_i.shape, map_i.wcs)

        # Define exclusion region
        ell_min = ell_peak * (1 - exclusion_width)
        ell_max = ell_peak * (1 + exclusion_width)

        # Create exclusion mask (0 in exclusion region, 1 outside)
        exclusion_mask = np.ones_like(modlmap)
        exclusion_mask[(modlmap >= ell_min) & (modlmap <= ell_max)] = 0.0

        # Apply mask in Fourier space
        fmap_filtered = fmap * exclusion_mask

        # Transform back to real space
        map_filtered = fft.ifft(fmap_filtered, axes=[-2, -1]).real

        filtered_maps.append(enmap.ndmap(map_filtered, map_i.wcs))

    return filtered_maps


def apply_donut_smoothing(needlet_maps, ell_peak, mask, smooth_scale_factor=2.0):
    """
    Apply "donut" smoothing for small-scale bias mitigation.

    This method is used for small scales (high ell_peak). It computes ILC
    weights from smoothed maps but applies them to unsmoothed data, creating
    a "donut" in Fourier space where weights don't see the exact modes they
    clean.

    Parameters
    ----------
    needlet_maps : list of enmap.ndmap
        Needlet-decomposed maps for each frequency
    ell_peak : float
        Peak multipole of the needlet scale
    mask : enmap.ndmap
        Analysis mask
    smooth_scale_factor : float
        Factor to determine smoothing scale relative to needlet scale
        smoothing_ell ~ ell_peak / smooth_scale_factor

    Returns
    -------
    smoothed_maps : list of enmap.ndmap
        Smoothed versions of input maps for covariance computation
    original_maps : list of enmap.ndmap
        Original unsmoothed maps (for applying weights)

    Notes
    -----
    Following Section III-C of arXiv:2307.01258, small scales use real-space
    donut smoothing. The ILC weights are computed from smoothed maps to avoid
    reusing the same modes, but applied to unsmoothed maps.

    This is more efficient than harmonic exclusion for small scales because
    it preserves more information while still mitigating bias.
    """
    if not PIXELL_AVAILABLE:
        raise ImportError("pixell is required for donut smoothing")

    # Determine smoothing scale from ell_peak
    # Smoothing scale in ell-space
    smooth_ell = ell_peak / smooth_scale_factor

    # Convert to real-space smoothing scale (sigma in radians)
    # For a Gaussian: ell ~ 1/sigma, so sigma ~ 1/ell
    # More precisely, FWHM = 2*sqrt(2*ln(2)) * sigma ~ 2.355 * sigma
    # and ell_FWHM ~ 2*pi/FWHM
    if smooth_ell > 0:
        sigma_rad = 1.0 / smooth_ell  # Rough estimate
    else:
        sigma_rad = 0.1  # Default for very large scales

    smoothed_maps = []
    original_maps = list(needlet_maps)  # Keep original references

    for map_i in needlet_maps:
        # Create Gaussian smoothing kernel in Fourier space
        modlmap = enmap.modlmap(map_i.shape, map_i.wcs)

        # Gaussian kernel: exp(-ell^2 * sigma^2 / 2)
        smooth_kernel = np.exp(-0.5 * modlmap**2 * sigma_rad**2)

        # Apply smoothing in Fourier space
        fmap = fft.fft(map_i, axes=[-2, -1])
        fmap_smooth = fmap * smooth_kernel
        map_smooth = fft.ifft(fmap_smooth, axes=[-2, -1]).real

        smoothed_maps.append(enmap.ndmap(map_smooth, map_i.wcs))

    return smoothed_maps, original_maps


def choose_bias_mitigation_strategy(ell_peak, transition_ell=1000):
    """
    Choose between harmonic exclusion and donut smoothing based on scale.

    Parameters
    ----------
    ell_peak : float
        Peak multipole of the needlet scale
    transition_ell : float
        Multipole at which to transition from harmonic to donut
        Default 1000 (suggested in the prompt)

    Returns
    -------
    strategy : str
        Either 'harmonic' or 'donut'

    Notes
    -----
    Large scales (ell_peak < transition_ell) use harmonic-space exclusion.
    Small scales (ell_peak >= transition_ell) use real-space donut smoothing.

    The transition scale is chosen based on where real-space and harmonic-space
    methods are most effective. At large scales (low ell), harmonic exclusion
    is cleaner. At small scales (high ell), donut smoothing is more efficient.
    """
    if ell_peak < transition_ell:
        return 'harmonic'
    else:
        return 'donut'


def apply_bias_mitigation(needlet_maps, ell_peak, mask, transition_ell=1000,
                          exclusion_width=0.3, smooth_scale_factor=2.0):
    """
    Apply appropriate bias mitigation strategy based on scale.

    This is a convenience function that automatically chooses and applies
    the appropriate bias mitigation method.

    Parameters
    ----------
    needlet_maps : list of enmap.ndmap
        Needlet-decomposed maps for each frequency
    ell_peak : float
        Peak multipole of the needlet scale
    mask : enmap.ndmap
        Analysis mask
    transition_ell : float
        Transition multipole between strategies
    exclusion_width : float
        Width parameter for harmonic exclusion
    smooth_scale_factor : float
        Smoothing parameter for donut method

    Returns
    -------
    processed_maps : list of enmap.ndmap
        Processed maps for covariance computation
    original_maps : list of enmap.ndmap
        Original maps for applying weights (relevant for donut method)
    strategy : str
        Which strategy was used

    Notes
    -----
    For harmonic exclusion, processed_maps are the filtered maps and should
    be used for both covariance and weight application.

    For donut smoothing, processed_maps are smoothed (for covariance) and
    original_maps are unsmoothed (for weight application).
    """
    strategy = choose_bias_mitigation_strategy(ell_peak, transition_ell)

    if strategy == 'harmonic':
        processed_maps = apply_harmonic_exclusion(
            needlet_maps, ell_peak, exclusion_width
        )
        original_maps = processed_maps  # Same for harmonic method
    else:  # strategy == 'donut'
        processed_maps, original_maps = apply_donut_smoothing(
            needlet_maps, ell_peak, mask, smooth_scale_factor
        )

    return processed_maps, original_maps, strategy
