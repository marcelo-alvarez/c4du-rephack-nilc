"""Fourier-domain filtering for ACT maps."""
import numpy as np
from pixell import enmap, fft


def apply_fourier_filter(map_data, lmin=None, lmax=None, filter_type='highpass'):
    """
    Apply Fourier-domain filtering to remove scan-synchronous modes.

    Parameters
    ----------
    map_data : enmap.ndmap
        Input map in CAR projection
    lmin : float, optional
        Minimum ell for filter (for highpass)
    lmax : float, optional
        Maximum ell for filter (for lowpass)
    filter_type : str
        'highpass', 'lowpass', or 'bandpass'

    Returns
    -------
    filtered_map : enmap.ndmap
        Filtered map
    """
    # Get the 2D Fourier transform
    fmap = fft.fft(map_data, axes=[-2, -1])

    # Get the modlmap (ell values for each mode)
    modlmap = enmap.modlmap(map_data.shape, map_data.wcs)

    # Create filter in Fourier space
    filter_2d = np.ones_like(modlmap)

    if filter_type == 'highpass' and lmin is not None:
        # High-pass filter: remove modes below lmin
        filter_2d[modlmap < lmin] = 0.0
    elif filter_type == 'lowpass' and lmax is not None:
        # Low-pass filter: remove modes above lmax
        filter_2d[modlmap > lmax] = 0.0
    elif filter_type == 'bandpass' and lmin is not None and lmax is not None:
        # Band-pass filter
        filter_2d[(modlmap < lmin) | (modlmap > lmax)] = 0.0

    # Apply filter and inverse FFT
    filtered_fmap = fmap * filter_2d
    filtered_map = fft.ifft(filtered_fmap, axes=[-2, -1]).real

    return enmap.ndmap(filtered_map, map_data.wcs)


def remove_scan_modes(map_data, ell_cut=100):
    """
    Remove scan-synchronous modes using high-pass filter.

    Parameters
    ----------
    map_data : enmap.ndmap
        Input map
    ell_cut : float
        Minimum ell to keep (removes modes below this)

    Returns
    -------
    filtered_map : enmap.ndmap
        Map with scan modes removed
    """
    return apply_fourier_filter(map_data, lmin=ell_cut, filter_type='highpass')
