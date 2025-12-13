"""End-to-end needlet-ILC pipeline for CMB extraction.

Implements the NILC pipeline described in arXiv:2307.01258:
1. Load multi-frequency ACT maps and analysis mask
2. Apply preprocessing (filtering, beam corrections)
3. Perform needlet decomposition
4. Compute ILC weights per needlet scale for CMB extraction
5. Apply ILC weights and recombine to produce CMB map
6. Save outputs in CAR and HEALPix formats
"""

import os
import numpy as np
from pixell import enmap, fft, reproject
import healpy as hp

from .decomposition import (
    needlet_decompose,
    needlet_recombine,
    ELL_PEAKS,
    get_scale_arcmin,
)
from ..ilc.weights import compute_ilc_weights, apply_ilc_weights
from ..ilc.covariance import compute_ilc_covariance
from ..ilc.frequency_response import compute_response_matrix


# Data paths
MAP_DIR = "/global/cfs/cdirs/act/data/act_dr6/dr6.02/maps/published"
MASK_PATH = "/global/cfs/cdirs/act/data/act_dr6/dr6.02/nilc/published/ilc_footprint_mask.fits"
BEAM_DIR = "/global/cfs/cdirs/act/data/act_dr6/dr6.02/beams/main_beams/nominal"

# Map file patterns
MAP_90_FILE = "act_dr4dr6_coadd_AA_daynight_f090_map_srcfree.fits"
MAP_150_FILE = "act_dr4dr6_coadd_AA_daynight_f150_map_srcfree.fits"

# Frequencies in GHz
FREQ_90 = 90.0
FREQ_150 = 150.0


def load_act_maps(map_dir=MAP_DIR, component='I'):
    """Load ACT 90 and 150 GHz source-free maps.

    Parameters
    ----------
    map_dir : str
        Directory containing ACT maps
    component : str
        Stokes component to extract: 'I', 'Q', or 'U'

    Returns
    -------
    map_90 : enmap.ndmap
        90 GHz map
    map_150 : enmap.ndmap
        150 GHz map
    """
    print(f"Loading ACT maps from {map_dir}")

    map_90_path = os.path.join(map_dir, MAP_90_FILE)
    map_150_path = os.path.join(map_dir, MAP_150_FILE)

    print(f"  Loading 90 GHz: {MAP_90_FILE}")
    map_90_full = enmap.read_map(map_90_path)

    print(f"  Loading 150 GHz: {MAP_150_FILE}")
    map_150_full = enmap.read_map(map_150_path)

    # Extract requested Stokes component
    component_idx = {'I': 0, 'Q': 1, 'U': 2}
    idx = component_idx.get(component, 0)

    if map_90_full.ndim == 3:
        map_90 = map_90_full[idx]
        map_150 = map_150_full[idx]
        print(f"  Extracted Stokes {component}")
    else:
        map_90 = map_90_full
        map_150 = map_150_full

    print(f"  Map shape: {map_90.shape}")
    return map_90, map_150


def load_mask(mask_path=MASK_PATH, geometry=None):
    """Load analysis mask and optionally reproject to map geometry.

    Parameters
    ----------
    mask_path : str
        Path to mask FITS file
    geometry : tuple, optional
        (shape, wcs) to reproject mask to

    Returns
    -------
    mask : enmap.ndmap or ndarray
        Analysis mask (1=good, 0=bad)
    """
    print(f"Loading mask from {mask_path}")

    # Check if mask is HEALPix or CAR
    try:
        mask = enmap.read_map(mask_path)
        print(f"  Loaded CAR mask, shape: {mask.shape}")
    except Exception:
        # Try loading as HEALPix
        mask_hp = hp.read_map(mask_path)
        print(f"  Loaded HEALPix mask, nside: {hp.npix2nside(len(mask_hp))}")

        if geometry is not None:
            shape, wcs = geometry
            # Reproject HEALPix to CAR
            mask = reproject.healpix2map(
                mask_hp, shape, wcs,
                rot=None, order=0, ncomp=1
            )
            if mask.ndim > 2:
                mask = mask[0]
            print(f"  Reprojected to shape: {mask.shape}")
        else:
            mask = mask_hp

    # Ensure binary mask
    mask = np.where(mask > 0.5, 1.0, 0.0)
    frac = np.mean(mask)
    print(f"  Sky fraction: {frac*100:.1f}%")

    return mask


def load_beam(freq_ghz, beam_dir=BEAM_DIR):
    """Load beam transfer function for a given frequency.

    Parameters
    ----------
    freq_ghz : float
        Frequency in GHz (90 or 150)
    beam_dir : str
        Directory containing beam files

    Returns
    -------
    ell : ndarray
        Multipole values
    beam : ndarray
        Beam transfer function B(ell)
    """
    freq_str = f"f{int(freq_ghz):03d}"
    beam_file = os.path.join(
        beam_dir,
        f"coadd_pa5_{freq_str}_night_beam_tform_jitter_cmb.txt"
    )

    if os.path.exists(beam_file):
        print(f"  Loading beam from {os.path.basename(beam_file)}")
        data = np.loadtxt(beam_file)
        ell = data[:, 0]
        beam = data[:, 1]
    else:
        print(f"  Beam file not found, using Gaussian approximation")
        # Use approximate Gaussian beam
        ell = np.arange(0, 20000)
        fwhm_arcmin = 1.4 if freq_ghz > 100 else 2.1
        sigma = fwhm_arcmin / 2.355 / 60.0 * np.pi / 180.0
        beam = np.exp(-0.5 * ell * (ell + 1) * sigma**2)

    return ell, beam


def apply_preprocessing(map_90, map_150, mask, lmin_cut=100):
    """Apply preprocessing steps to frequency maps.

    Currently implements:
    - High-pass Fourier filter to remove scan-synchronous modes

    Parameters
    ----------
    map_90 : enmap.ndmap
        90 GHz map
    map_150 : enmap.ndmap
        150 GHz map
    mask : ndarray
        Analysis mask
    lmin_cut : int
        Minimum ell to keep (removes ell < lmin_cut)

    Returns
    -------
    map_90_proc : enmap.ndmap
        Preprocessed 90 GHz map
    map_150_proc : enmap.ndmap
        Preprocessed 150 GHz map
    """
    print(f"Applying preprocessing (lmin_cut={lmin_cut})")

    # Get modlmap for filtering
    modlmap = enmap.modlmap(map_90.shape, map_90.wcs)

    # Create high-pass filter
    filter_2d = np.ones_like(modlmap)
    filter_2d[modlmap < lmin_cut] = 0.0

    # Apply to 90 GHz
    fmap_90 = fft.fft(map_90, axes=[-2, -1])
    fmap_90_filt = fmap_90 * filter_2d
    map_90_proc = fft.ifft(fmap_90_filt, axes=[-2, -1]).real
    map_90_proc = enmap.ndmap(map_90_proc, map_90.wcs)

    # Apply to 150 GHz
    fmap_150 = fft.fft(map_150, axes=[-2, -1])
    fmap_150_filt = fmap_150 * filter_2d
    map_150_proc = fft.ifft(fmap_150_filt, axes=[-2, -1]).real
    map_150_proc = enmap.ndmap(map_150_proc, map_150.wcs)

    print("  High-pass filter applied")
    return map_90_proc, map_150_proc


def run_nilc_pipeline(
    output_dir="./output",
    component='I',
    lmin_cut=100,
    target_component='cmb',
    regularization=1e-10,
    save_intermediate=False,
):
    """Run the full needlet-ILC pipeline for CMB extraction.

    Parameters
    ----------
    output_dir : str
        Directory to save output maps
    component : str
        Stokes component to process: 'I', 'Q', or 'U'
    lmin_cut : int
        Minimum ell to keep in preprocessing
    target_component : str
        Component to extract: 'cmb' or 'ksz'
    regularization : float
        Regularization for ILC weight computation
    save_intermediate : bool
        If True, save intermediate needlet maps

    Returns
    -------
    cmb_map : enmap.ndmap
        Extracted CMB map in CAR projection
    """
    print("=" * 70)
    print("NEEDLET-ILC PIPELINE")
    print(f"Target component: {target_component.upper()}")
    print(f"Stokes component: {component}")
    print("=" * 70)
    print()

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Step 1: Load maps
    print("Step 1: Loading frequency maps")
    print("-" * 40)
    map_90, map_150 = load_act_maps(component=component)
    geometry = (map_90.shape, map_90.wcs)
    print()

    # Step 2: Load mask
    print("Step 2: Loading analysis mask")
    print("-" * 40)
    mask = load_mask(geometry=geometry)
    if isinstance(mask, np.ndarray) and not hasattr(mask, 'wcs'):
        mask = enmap.ndmap(mask, map_90.wcs)
    print()

    # Step 3: Apply preprocessing
    print("Step 3: Preprocessing")
    print("-" * 40)
    map_90_proc, map_150_proc = apply_preprocessing(
        map_90, map_150, mask, lmin_cut=lmin_cut
    )
    print()

    # Step 4: Needlet decomposition
    print("Step 4: Needlet decomposition")
    print("-" * 40)
    print("  Decomposing 90 GHz map...")
    needlet_90, ell_peaks = needlet_decompose(map_90_proc)
    print("  Decomposing 150 GHz map...")
    needlet_150, _ = needlet_decompose(map_150_proc)
    n_scales = len(ell_peaks)
    print(f"  Decomposed into {n_scales} needlet scales")
    print()

    # Step 5: ILC weights and component separation per scale
    print("Step 5: Computing ILC weights per needlet scale")
    print("-" * 40)

    freq_list = [FREQ_90, FREQ_150]

    # Build frequency response matrix for CMB only (no deprojection)
    # With only 2 frequencies, we can only constrain 1 component
    F, comp_names = compute_response_matrix(freq_list, components=[target_component])
    e = np.array([1.0])  # Preserve target component

    cmb_needlet_maps = []

    for j in range(n_scales):
        ell_peak = ell_peaks[j]
        scale_arcmin = get_scale_arcmin(ell_peak)

        # Stack the two frequency needlet maps
        maps_j = np.stack([
            np.asarray(needlet_90[j]),
            np.asarray(needlet_150[j])
        ], axis=0)

        # Get mask as numpy array
        mask_arr = np.asarray(mask)

        # Compute covariance with bias mitigation
        if np.isinf(scale_arcmin):
            # DC mode, use simple covariance
            C, method = compute_ilc_covariance(
                maps_j, mask_arr, scale_arcmin=60.0,
                method='none'
            )
        else:
            C, method = compute_ilc_covariance(
                maps_j, mask_arr, scale_arcmin=min(scale_arcmin, 120.0),
                method='auto'
            )

        # Compute ILC weights
        # For simple CMB extraction with 2 frequencies and 1 constraint,
        # the weights are determined by the covariance and unit response
        try:
            w = compute_ilc_weights(C, F, e, regularization=regularization)
        except Exception as err:
            print(f"  Scale {j} (ell={ell_peak}): weight computation failed ({err})")
            # Fall back to equal weights
            w = np.array([0.5, 0.5])

        # Apply weights to extract CMB at this scale
        cmb_j = apply_ilc_weights(maps_j, w, mask=mask_arr)
        cmb_needlet_maps.append(enmap.ndmap(cmb_j, map_90.wcs))

        if j < 5 or j >= n_scales - 3:
            print(f"  Scale {j:2d}: ell_peak={ell_peak:5.0f}, "
                  f"w=[{w[0]:.3f}, {w[1]:.3f}], method={method}")
        elif j == 5:
            print("  ...")

    print()

    # Step 6: Recombine needlet scales
    print("Step 6: Recombining needlet scales")
    print("-" * 40)
    cmb_map = needlet_recombine(cmb_needlet_maps)
    print(f"  Output map shape: {cmb_map.shape}")
    print()

    # Step 7: Save outputs
    print("Step 7: Saving outputs")
    print("-" * 40)

    # Save CAR map
    car_filename = f"nilc_{target_component}_{component}.fits"
    car_path = os.path.join(output_dir, car_filename)
    enmap.write_map(car_path, cmb_map)
    print(f"  CAR map: {car_path}")

    # Convert to HEALPix and save
    nside = 4096  # Standard HEALPix resolution
    print(f"  Converting to HEALPix (nside={nside})...")
    cmb_hp = reproject.map2healpix(cmb_map, nside=nside, rot=None)
    hp_filename = f"nilc_{target_component}_{component}_healpix.fits"
    hp_path = os.path.join(output_dir, hp_filename)
    hp.write_map(hp_path, cmb_hp, overwrite=True)
    print(f"  HEALPix map: {hp_path}")

    # Save intermediate needlet maps if requested
    if save_intermediate:
        inter_dir = os.path.join(output_dir, "needlet_scales")
        os.makedirs(inter_dir, exist_ok=True)
        for j, nm in enumerate(cmb_needlet_maps):
            scale_path = os.path.join(inter_dir, f"needlet_scale_{j:02d}.fits")
            enmap.write_map(scale_path, nm)
        print(f"  Intermediate maps: {inter_dir}/")

    print()
    print("=" * 70)
    print("PIPELINE COMPLETE")
    print("=" * 70)

    return cmb_map


if __name__ == '__main__':
    # Run pipeline with default settings
    cmb = run_nilc_pipeline(
        output_dir="./output",
        component='I',
        target_component='cmb',
    )
