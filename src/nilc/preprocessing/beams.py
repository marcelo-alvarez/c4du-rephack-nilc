"""Beam correction utilities for ACT maps."""
import numpy as np
from pixell import enmap, fft
import os


def load_beam_profile(beam_file):
    """
    Load beam transfer function from file.

    Parameters
    ----------
    beam_file : str
        Path to beam FITS file

    Returns
    -------
    ell : ndarray
        Multipole values
    beam : ndarray
        Beam transfer function B(ell)
    """
    # Placeholder - actual implementation depends on beam file format
    # For now, create a simple Gaussian beam model
    ell = np.arange(0, 10000)
    fwhm_arcmin = 1.4  # Example FWHM for ACT
    sigma = fwhm_arcmin / 2.355 / 60.0 * np.pi / 180.0  # Convert to radians
    beam = np.exp(-0.5 * ell * (ell + 1) * sigma**2)
    return ell, beam


def apply_beam_correction(map_data, beam_ell, beam_profile, correction_type='deconvolve'):
    """
    Apply beam correction in harmonic space.

    Parameters
    ----------
    map_data : enmap.ndmap
        Input map
    beam_ell : ndarray
        Multipole values for beam
    beam_profile : ndarray
        Beam transfer function B(ell)
    correction_type : str
        'deconvolve' to sharpen, 'convolve' to smooth

    Returns
    -------
    corrected_map : enmap.ndmap
        Beam-corrected map
    """
    # Get 2D Fourier transform
    fmap = fft.fft(map_data, axes=[-2, -1])

    # Get modlmap
    modlmap = enmap.modlmap(map_data.shape, map_data.wcs)

    # Interpolate beam to 2D ell space
    beam_2d = np.interp(modlmap.flatten(), beam_ell, beam_profile).reshape(modlmap.shape)

    # Apply correction
    if correction_type == 'deconvolve':
        # Deconvolve beam (with regularization to avoid division by small numbers)
        epsilon = 1e-6
        correction = 1.0 / (beam_2d + epsilon)
        # Limit correction for high ell where beam -> 0
        correction[modlmap > beam_ell[-1]] = 0.0
    else:  # convolve
        correction = beam_2d

    # Apply and inverse FFT
    corrected_fmap = fmap * correction
    corrected_map = fft.ifft(corrected_fmap, axes=[-2, -1]).real

    return enmap.ndmap(corrected_map, map_data.wcs)


def match_beams(map1, map2, beam1_ell, beam1_prof, beam2_ell, beam2_prof,
                target_fwhm_arcmin=None):
    """
    Match beams between two frequency maps.

    Parameters
    ----------
    map1, map2 : enmap.ndmap
        Input maps at two frequencies
    beam1_ell, beam1_prof : ndarray
        Beam for map1
    beam2_ell, beam2_prof : ndarray
        Beam for map2
    target_fwhm_arcmin : float, optional
        Target FWHM in arcminutes. If None, uses larger of the two beams.

    Returns
    -------
    map1_matched, map2_matched : enmap.ndmap
        Beam-matched maps
    """
    # For now, deconvolve both to infinite resolution, then convolve to common beam
    # This is a simplified approach; real implementation should be more careful

    # Deconvolve original beams
    map1_deconv = apply_beam_correction(map1, beam1_ell, beam1_prof, 'deconvolve')
    map2_deconv = apply_beam_correction(map2, beam2_ell, beam2_prof, 'deconvolve')

    # Create target beam
    if target_fwhm_arcmin is None:
        # Use the larger FWHM as target
        target_fwhm_arcmin = 1.4  # Example value

    ell = np.arange(0, 10000)
    sigma = target_fwhm_arcmin / 2.355 / 60.0 * np.pi / 180.0
    target_beam = np.exp(-0.5 * ell * (ell + 1) * sigma**2)

    # Convolve to target beam
    map1_matched = apply_beam_correction(map1_deconv, ell, target_beam, 'convolve')
    map2_matched = apply_beam_correction(map2_deconv, ell, target_beam, 'convolve')

    return map1_matched, map2_matched
