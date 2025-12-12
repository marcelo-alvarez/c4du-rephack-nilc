"""
Covariance matrix computation for ILC component separation.

Implements in-band covariance computation per Section III-C of arXiv:2307.01258.
"""
import numpy as np
from pixell import enmap, fft


def compute_covariance_map(needlet_coeffs, mask=None):
    """
    Compute spatial covariance matrix from needlet coefficients.

    For each pixel (or region), computes the covariance between frequency channels
    from the needlet coefficients at that scale.

    Parameters
    ----------
    needlet_coeffs : ndarray
        Needlet coefficients with shape (n_freq, n_y, n_x)
        First dimension indexes frequency channels (2 for ACT 90/150 GHz)
    mask : enmap.ndmap, optional
        Analysis mask (1=use, 0=exclude)

    Returns
    -------
    cov_map : ndarray
        Covariance matrix at each pixel, shape (n_y, n_x, n_freq, n_freq)
        For 2 frequencies, this is (n_y, n_x, 2, 2)

    Notes
    -----
    This computes the simple pixel-wise covariance:
        C_ij = <d_i * d_j>
    where d_i is the needlet coefficient at frequency i.

    For bias mitigation, use compute_covariance_with_exclusion or
    compute_covariance_donut instead.
    """
    n_freq = needlet_coeffs.shape[0]
    n_y, n_x = needlet_coeffs.shape[1:3]

    # Initialize covariance map
    cov_map = np.zeros((n_y, n_x, n_freq, n_freq))

    # Compute covariance at each pixel
    for i in range(n_freq):
        for j in range(n_freq):
            cov_map[:, :, i, j] = needlet_coeffs[i] * needlet_coeffs[j]

    # Apply mask if provided
    if mask is not None:
        mask_2d = mask if mask.ndim == 2 else mask[0]  # Handle (ncomp, ny, nx) or (ny, nx)
        for i in range(n_freq):
            for j in range(n_freq):
                cov_map[:, :, i, j] *= mask_2d

    return cov_map


def compute_covariance_smoothed(needlet_coeffs, smooth_scale_arcmin, wcs, mask=None):
    """
    Compute spatially-smoothed covariance matrix.

    Smooths the covariance estimate by averaging over a region, reducing
    noise in the covariance estimate at the cost of spatial resolution.

    Parameters
    ----------
    needlet_coeffs : ndarray
        Needlet coefficients with shape (n_freq, n_y, n_x)
    smooth_scale_arcmin : float
        Smoothing scale in arcminutes (FWHM of Gaussian kernel)
    wcs : astropy.wcs.WCS
        WCS information for the map geometry
    mask : enmap.ndmap, optional
        Analysis mask

    Returns
    -------
    cov_smoothed : ndarray
        Smoothed covariance matrix, shape (n_y, n_x, n_freq, n_freq)
    """
    n_freq = needlet_coeffs.shape[0]
    shape = needlet_coeffs.shape[1:3]

    # Compute pixel-wise covariance
    cov_map = compute_covariance_map(needlet_coeffs, mask=mask)

    # Create Gaussian smoothing kernel in Fourier space
    modlmap = enmap.modlmap(shape, wcs)
    sigma_rad = (smooth_scale_arcmin / 2.355) / 60.0 * np.pi / 180.0
    smooth_kernel = np.exp(-0.5 * modlmap**2 * sigma_rad**2)

    # Smooth each component of covariance matrix
    cov_smoothed = np.zeros_like(cov_map)
    for i in range(n_freq):
        for j in range(n_freq):
            # Create enmap for this component
            cov_component = enmap.ndmap(cov_map[:, :, i, j], wcs)

            # Fourier transform
            fmap = fft.fft(cov_component, axes=[-2, -1])

            # Apply smoothing
            fmap_smooth = fmap * smooth_kernel

            # Inverse transform
            cov_smooth = fft.ifft(fmap_smooth, axes=[-2, -1]).real
            cov_smoothed[:, :, i, j] = cov_smooth

    return cov_smoothed


def compute_covariance_donut(needlet_coeffs, inner_radius_arcmin, outer_radius_arcmin,
                             wcs, mask=None):
    """
    Compute covariance using "donut" exclusion for bias mitigation.

    For each pixel, computes covariance using data from an annular region
    (donut) around that pixel, excluding the central region. This prevents
    the ILC from using the same modes it's applied to, reducing bias.

    Parameters
    ----------
    needlet_coeffs : ndarray
        Needlet coefficients with shape (n_freq, n_y, n_x)
    inner_radius_arcmin : float
        Inner radius of donut in arcminutes (excludes central region)
    outer_radius_arcmin : float
        Outer radius of donut in arcminutes
    wcs : astropy.wcs.WCS
        WCS information for map geometry
    mask : enmap.ndmap, optional
        Analysis mask

    Returns
    -------
    cov_donut : ndarray
        Covariance matrix computed with donut exclusion, shape (n_y, n_x, n_freq, n_freq)

    Notes
    -----
    This implements the real-space bias mitigation strategy for small angular scales
    described in Section III-C of arXiv:2307.01258.

    For computational efficiency, this uses a simplified approach where we smooth
    the covariance then subtract the central contribution. A more rigorous
    implementation would explicitly compute the annular average.
    """
    # For now, implement simplified version:
    # smooth at outer scale, subtract inner contribution

    # Compute smoothed covariance at outer scale
    cov_outer = compute_covariance_smoothed(needlet_coeffs, outer_radius_arcmin * 2.355,
                                           wcs, mask=mask)

    # Compute smoothed covariance at inner scale
    if inner_radius_arcmin > 0:
        cov_inner = compute_covariance_smoothed(needlet_coeffs, inner_radius_arcmin * 2.355,
                                               wcs, mask=mask)

        # Donut = outer - inner (approximate)
        # This is a simplified version; proper implementation needs area weighting
        cov_donut = cov_outer - cov_inner
    else:
        cov_donut = cov_outer

    return cov_donut


def invert_covariance(cov_map, regularization=1e-10):
    """
    Invert covariance matrices with regularization.

    Parameters
    ----------
    cov_map : ndarray
        Covariance matrices with shape (n_y, n_x, n_freq, n_freq)
    regularization : float
        Regularization parameter added to diagonal for numerical stability

    Returns
    -------
    cov_inv : ndarray
        Inverted covariance matrices, same shape as input

    Notes
    -----
    Adds small regularization to diagonal before inversion to handle
    singular or near-singular matrices.
    """
    n_y, n_x, n_freq, _ = cov_map.shape
    cov_inv = np.zeros_like(cov_map)

    # Add regularization to diagonal
    cov_reg = cov_map.copy()
    for i in range(n_freq):
        cov_reg[:, :, i, i] += regularization

    # Invert at each pixel
    for iy in range(n_y):
        for ix in range(n_x):
            try:
                cov_inv[iy, ix] = np.linalg.inv(cov_reg[iy, ix])
            except np.linalg.LinAlgError:
                # If inversion fails, use pseudoinverse
                cov_inv[iy, ix] = np.linalg.pinv(cov_reg[iy, ix])

    return cov_inv
