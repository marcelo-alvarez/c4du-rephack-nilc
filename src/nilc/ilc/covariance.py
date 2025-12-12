"""Covariance matrix computation for ILC component separation."""
import numpy as np

try:
    from pixell import enmap, fft
    PIXELL_AVAILABLE = True
except ImportError:
    PIXELL_AVAILABLE = False
    enmap = None
    fft = None


def compute_inband_covariance(needlet_maps, mask):
    """
    Compute the in-band covariance matrix for a given needlet scale.

    This function computes the spatial covariance between frequency maps
    at a specific needlet scale, accounting for the analysis mask.

    Parameters
    ----------
    needlet_maps : list of enmap.ndmap
        Needlet-decomposed maps for each frequency at a specific scale.
        List of length n_freq, each element is a 2D map.
    mask : enmap.ndmap
        Analysis mask (1 = include, 0 = exclude). Should have same shape
        as needlet_maps.

    Returns
    -------
    cov : ndarray
        Covariance matrix of shape (n_freq, n_freq)

    Notes
    -----
    Following Section III-C of arXiv:2307.01258, the covariance is computed
    as the spatial correlation between needlet maps, weighted by the mask
    to avoid biasing from masked regions.

    The covariance C_ij = <map_i * map_j> where <> denotes spatial average
    over the unmasked region.
    """
    n_freq = len(needlet_maps)
    cov = np.zeros((n_freq, n_freq))

    # Ensure mask is boolean or binary
    mask_binary = (mask > 0.5).astype(float)
    n_pix = np.sum(mask_binary)

    if n_pix == 0:
        raise ValueError("Mask contains no valid pixels")

    # Compute covariance for each pair of frequencies
    for i in range(n_freq):
        for j in range(i, n_freq):  # Symmetric, only compute upper triangle
            # Extract maps and apply mask
            map_i = needlet_maps[i] * mask_binary
            map_j = needlet_maps[j] * mask_binary

            # Compute spatial covariance (mean of product)
            # This is equivalent to <x_i * x_j> over the unmasked region
            cov_ij = np.sum(map_i * map_j) / n_pix

            cov[i, j] = cov_ij
            if i != j:
                cov[j, i] = cov_ij  # Symmetry

    return cov


def compute_inband_covariance_smoothed(needlet_maps, mask, smoothing_scale=None):
    """
    Compute covariance from smoothed needlet maps for bias mitigation.

    This is used in the "donut" smoothing approach for small-scale bias mitigation,
    where weights are computed from smoothed maps but applied to unsmoothed data.

    Parameters
    ----------
    needlet_maps : list of enmap.ndmap
        Needlet-decomposed maps for each frequency
    mask : enmap.ndmap
        Analysis mask
    smoothing_scale : float, optional
        Smoothing scale in radians. If None, uses a scale appropriate
        for the needlet band.

    Returns
    -------
    cov : ndarray
        Covariance matrix from smoothed maps

    Notes
    -----
    The smoothing prevents weights from being correlated with the same modes
    they will be applied to, reducing ILC bias per Section III-C.
    """
    if not PIXELL_AVAILABLE:
        raise ImportError("pixell is required for covariance smoothing")

    n_freq = len(needlet_maps)

    # If no smoothing scale provided, estimate from map resolution
    if smoothing_scale is None:
        # Use a conservative smoothing of ~2 arcmin
        smoothing_scale = 2.0 / 60.0 * np.pi / 180.0  # Convert arcmin to radians

    # Smooth each map with a Gaussian kernel
    smoothed_maps = []
    for map_i in needlet_maps:
        # Create Gaussian smoothing kernel in Fourier space
        modlmap = enmap.modlmap(map_i.shape, map_i.wcs)
        ell = modlmap.flatten()

        # Gaussian beam: B(ell) = exp(-ell^2 * sigma^2 / 2)
        sigma = smoothing_scale
        smooth_kernel = np.exp(-0.5 * modlmap**2 * sigma**2)

        # Apply smoothing in Fourier space
        fmap = fft.fft(map_i, axes=[-2, -1])
        fmap_smooth = fmap * smooth_kernel
        map_smooth = fft.ifft(fmap_smooth, axes=[-2, -1]).real

        smoothed_maps.append(enmap.ndmap(map_smooth, map_i.wcs))

    # Compute covariance from smoothed maps
    cov = compute_inband_covariance(smoothed_maps, mask)

    return cov
