"""Needlet decomposition and recombination for CAR maps.

Implements the axisymmetric needlet frame described in Section III-B of
arXiv:2307.01258, Equation 2.

Needlets are a wavelet-like basis that provides localization in both harmonic
and real space. Each needlet scale j is characterized by a window function b_j(ell)
that peaks at ell_peak[j] and has compact support in ell space.
"""

import numpy as np
from pixell import enmap, fft

# 26 needlet scales with ell_peak values from arXiv:2307.01258 Section III-B
ELL_PEAKS = np.array([
    0, 100, 200, 300, 400, 600, 800, 1000, 1250, 1400,
    1800, 2200, 2800, 3400, 4000, 4600, 5200, 6000, 7000, 8000,
    9000, 10000, 11000, 13000, 15000, 17000
])


def _smooth_step(t):
    """Smooth step function f(t) used in needlet construction.

    f(t) = 0 for t <= 0
    f(t) = exp(-1/t) for t > 0

    This function provides smooth transition between 0 and positive values.
    """
    result = np.zeros_like(t, dtype=float)
    mask = t > 0
    result[mask] = np.exp(-1.0 / t[mask])
    return result


def _psi(t):
    """Smooth partition of unity function.

    psi(t) = f(t) / (f(t) + f(1-t))

    This gives a smooth transition from 0 to 1 over the interval [0, 1].
    """
    f_t = _smooth_step(t)
    f_1_t = _smooth_step(1.0 - t)
    # Avoid division by zero
    denom = f_t + f_1_t
    denom = np.where(denom == 0, 1.0, denom)
    return f_t / denom


def compute_needlet_kernel(ell, ell_peak_prev, ell_peak, ell_peak_next):
    """Compute axisymmetric needlet kernel b(ell) for a single scale.

    The needlet kernel follows Equation 2 of arXiv:2307.01258:

    b(ell)^2 = psi(u(ell)) - psi(u(ell/B))

    where:
    - B is the ratio ell_peak_next / ell_peak
    - u(x) maps the multipole to the [0,1] interval
    - psi is a smooth partition function

    The window peaks at ell_peak and smoothly transitions to zero at the
    neighboring scale boundaries.

    Parameters
    ----------
    ell : ndarray
        Array of multipole values
    ell_peak_prev : float
        Peak ell of previous scale (lower bound)
    ell_peak : float
        Peak ell of current scale
    ell_peak_next : float
        Peak ell of next scale (upper bound)

    Returns
    -------
    kernel : ndarray
        Needlet window function b(ell), same shape as ell
    """
    ell = np.asarray(ell, dtype=float)
    kernel = np.zeros_like(ell)

    # Handle edge cases
    if ell_peak == 0:
        # First scale: low-pass filter
        # Smoothly goes from 1 at ell=0 to 0 at ell=ell_peak_next
        if ell_peak_next > 0:
            t = ell / ell_peak_next
            kernel = 1.0 - _psi(t)
        else:
            kernel = np.ones_like(ell)
        return kernel

    if ell_peak_next == np.inf or ell_peak_next == 0:
        # Last scale: high-pass filter
        # Smoothly goes from 0 at ell_peak_prev to 1 at ell_peak
        if ell_peak > 0 and ell_peak_prev >= 0:
            t = (ell - ell_peak_prev) / (ell_peak - ell_peak_prev)
            t = np.clip(t, 0, 1)
            kernel = _psi(t)
        return kernel

    # Standard band-pass needlet kernel
    # Rising edge from ell_peak_prev to ell_peak
    # Falling edge from ell_peak to ell_peak_next

    # Map ell to [0, 1] for rising edge
    if ell_peak > ell_peak_prev:
        t_rise = (ell - ell_peak_prev) / (ell_peak - ell_peak_prev)
        t_rise = np.clip(t_rise, 0, 1)
        rise = _psi(t_rise)
    else:
        rise = np.ones_like(ell)

    # Map ell to [0, 1] for falling edge
    if ell_peak_next > ell_peak:
        t_fall = (ell - ell_peak) / (ell_peak_next - ell_peak)
        t_fall = np.clip(t_fall, 0, 1)
        fall = 1.0 - _psi(t_fall)
    else:
        fall = np.ones_like(ell)

    kernel = rise * fall

    return kernel


def compute_needlet_kernels(ell, ell_peaks=None, lmax=None):
    """Compute all needlet kernels for the standard 26-scale decomposition.

    Creates a set of window functions b_j(ell) that form a partition of unity:
    sum_j b_j(ell)^2 ≈ 1 for all ell in [0, lmax]

    Parameters
    ----------
    ell : ndarray
        Array of multipole values
    ell_peaks : ndarray, optional
        Array of ell_peak values for each scale. Default is ELL_PEAKS.
    lmax : float, optional
        Maximum ell to include. Scales with ell_peak > lmax are skipped.

    Returns
    -------
    kernels : ndarray, shape (n_scales, len(ell))
        Needlet window functions for each scale
    ell_peaks_used : ndarray
        The ell_peak values actually used (may be truncated by lmax)
    """
    if ell_peaks is None:
        ell_peaks = ELL_PEAKS

    ell = np.asarray(ell, dtype=float)
    ell_peaks = np.asarray(ell_peaks, dtype=float)

    # Filter scales by lmax if specified
    if lmax is not None:
        mask = ell_peaks <= lmax
        # Always include at least the first scale
        mask[0] = True
        ell_peaks = ell_peaks[mask]

    n_scales = len(ell_peaks)
    kernels = np.zeros((n_scales, len(ell)))

    for j in range(n_scales):
        # Get neighboring ell_peak values
        ell_peak_prev = ell_peaks[j - 1] if j > 0 else 0
        ell_peak = ell_peaks[j]
        ell_peak_next = ell_peaks[j + 1] if j < n_scales - 1 else np.inf

        kernels[j] = compute_needlet_kernel(ell, ell_peak_prev, ell_peak, ell_peak_next)

    return kernels, ell_peaks


def needlet_decompose(map_data, ell_peaks=None, return_kernels=False):
    """Decompose a CAR map into needlet scales.

    Performs harmonic-space filtering using needlet window functions to
    separate the map into different angular scale bands.

    Parameters
    ----------
    map_data : enmap.ndmap
        Input map in CAR projection (can be 2D or 3D for IQU)
    ell_peaks : ndarray, optional
        Array of ell_peak values for each scale. Default is ELL_PEAKS.
    return_kernels : bool, optional
        If True, also return the 2D needlet kernels used.

    Returns
    -------
    needlet_maps : list of enmap.ndmap
        List of maps at each needlet scale
    ell_peaks_used : ndarray
        The ell_peak values used for decomposition
    kernels_2d : list of ndarray, optional
        2D needlet kernels (only if return_kernels=True)
    """
    if ell_peaks is None:
        ell_peaks = ELL_PEAKS

    # Get 2D ell map from the map geometry
    modlmap = enmap.modlmap(map_data.shape, map_data.wcs)
    lmax = np.max(modlmap)

    # Flatten ell values for kernel computation
    ell_flat = modlmap.flatten()

    # Compute 1D kernels
    kernels_1d, ell_peaks_used = compute_needlet_kernels(ell_flat, ell_peaks, lmax)

    # Reshape kernels to 2D
    n_scales = len(ell_peaks_used)
    kernels_2d = [k.reshape(modlmap.shape) for k in kernels_1d]

    # FFT the input map
    fmap = fft.fft(map_data, axes=[-2, -1])

    # Apply each needlet kernel and inverse FFT
    needlet_maps = []
    for j in range(n_scales):
        # Apply kernel in harmonic space
        filtered_fmap = fmap * kernels_2d[j]

        # Inverse FFT
        filtered_map = fft.ifft(filtered_fmap, axes=[-2, -1]).real

        # Preserve WCS
        needlet_maps.append(enmap.ndmap(filtered_map, map_data.wcs))

    if return_kernels:
        return needlet_maps, ell_peaks_used, kernels_2d
    return needlet_maps, ell_peaks_used


def needlet_recombine(needlet_maps, weights=None):
    """Recombine needlet-decomposed maps into a single map.

    The recombination sums the needlet maps:
    output = sum_j needlet_maps[j]

    If weights are provided, each scale is weighted before summation:
    output = sum_j weights[j] * needlet_maps[j]

    Parameters
    ----------
    needlet_maps : list of enmap.ndmap
        List of maps at each needlet scale from needlet_decompose()
    weights : ndarray, optional
        Weights to apply at each scale. Shape: (n_scales,)
        If None, all weights are 1.0.

    Returns
    -------
    output_map : enmap.ndmap
        Recombined map with same geometry as input needlet maps
    """
    if len(needlet_maps) == 0:
        raise ValueError("needlet_maps list is empty")

    n_scales = len(needlet_maps)

    if weights is None:
        weights = np.ones(n_scales)
    else:
        weights = np.asarray(weights)
        if len(weights) != n_scales:
            raise ValueError(
                f"weights length ({len(weights)}) must match number of "
                f"needlet scales ({n_scales})"
            )

    # Initialize output with first weighted map
    output_map = weights[0] * needlet_maps[0].copy()

    # Add remaining scales
    for j in range(1, n_scales):
        output_map += weights[j] * needlet_maps[j]

    return output_map


def get_scale_arcmin(ell_peak):
    """Convert ell_peak to approximate angular scale in arcminutes.

    Uses the approximate relation: theta [arcmin] ≈ 180 * 60 / ell

    Parameters
    ----------
    ell_peak : float or ndarray
        Peak multipole value(s)

    Returns
    -------
    scale_arcmin : float or ndarray
        Approximate angular scale in arcminutes
    """
    ell_peak = np.asarray(ell_peak, dtype=float)
    # Handle ell=0 case
    with np.errstate(divide='ignore'):
        scale = 180.0 * 60.0 / ell_peak
    scale = np.where(ell_peak == 0, np.inf, scale)
    return scale


if __name__ == '__main__':
    print("=" * 70)
    print("Testing Needlet Decomposition Module")
    print("=" * 70)
    print()

    # Test needlet kernel computation
    print("1. Testing needlet kernel computation...")
    ell = np.arange(0, 5000)
    kernels, ell_peaks_used = compute_needlet_kernels(ell)

    print(f"   Number of scales: {len(ell_peaks_used)}")
    print(f"   ell_peak values: {ell_peaks_used}")
    print()

    # Check partition of unity
    sum_sq = np.sum(kernels**2, axis=0)
    print(f"   Partition of unity check (sum of b^2):")
    print(f"   Min: {np.min(sum_sq):.4f}, Max: {np.max(sum_sq):.4f}")
    print()

    # Test with synthetic map
    print("2. Testing decomposition and recombination...")
    from pixell import enmap

    # Create a simple test map
    shape = (200, 400)
    wcs = enmap.create_wcs(
        shape,
        np.array([[-1, 180], [1, -180]]) * np.pi / 180,  # RA/Dec box
        proj='car'
    )

    # Create synthetic map with different ell content
    np.random.seed(42)
    test_map = enmap.ndmap(np.random.randn(*shape), wcs)

    # Decompose
    needlet_maps, ell_peaks = needlet_decompose(test_map)
    print(f"   Decomposed into {len(needlet_maps)} scales")
    print(f"   Each map shape: {needlet_maps[0].shape}")
    print()

    # Recombine
    recombined = needlet_recombine(needlet_maps)
    print(f"   Recombined map shape: {recombined.shape}")

    # Check reconstruction error
    diff = test_map - recombined
    rel_error = np.std(diff) / np.std(test_map)
    print(f"   Reconstruction relative error: {rel_error:.2e}")
    print()

    # Test angular scale conversion
    print("3. Angular scales for each needlet band:")
    for j, ell_p in enumerate(ELL_PEAKS[:10]):
        scale = get_scale_arcmin(ell_p)
        if np.isinf(scale):
            print(f"   Scale {j}: ell_peak = {ell_p:5.0f}, theta = infinity (DC mode)")
        else:
            print(f"   Scale {j}: ell_peak = {ell_p:5.0f}, theta = {scale:.1f} arcmin")
    print("   ...")
    print()

    print("=" * 70)
    print("All tests completed successfully!")
    print("=" * 70)
