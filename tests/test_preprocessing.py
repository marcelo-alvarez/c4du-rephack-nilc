"""Test preprocessing functions."""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from nilc.io.map_io import load_act_map, get_stokes_component
from nilc.preprocessing.filters import remove_scan_modes
from nilc.preprocessing.beams import load_beam_profile, apply_beam_correction


def test_preprocessing():
    """Test Fourier filtering and beam corrections."""
    # Load test map
    map_dir = "/scratch/jiaqu/hack_data/maps/"
    map_file = "act_dr6.02_std_AA_night_pa4_f220_4way_coadd_map_srcfree.fits"
    map_path = os.path.join(map_dir, map_file)

    # Load mask
    mask_path = "/scratch/jiaqu/hack_data/masks/ilc_footprint_mask.fits"

    print(f"Loading map: {map_path}")
    map_data = load_act_map(map_path)
    temp_map = get_stokes_component(map_data, 'I')

    print(f"Loading mask: {mask_path}")
    mask = load_act_map(mask_path)

    print(f"Map shape: {temp_map.shape}")
    print(f"Mask shape: {mask.shape}")

    # Apply mask
    masked_map = temp_map * mask
    print(f"Applied mask, non-zero pixels: {np.sum(mask > 0)}")

    # Test Fourier filtering
    print("\nApplying Fourier filtering to remove scan modes...")
    filtered_map = remove_scan_modes(masked_map, ell_cut=100)
    print(f"Filtered map statistics:")
    print(f"  Mean: {np.mean(filtered_map[mask > 0]):.3e}")
    print(f"  Std: {np.std(filtered_map[mask > 0]):.3e}")

    # Test beam correction (using example beam)
    print("\nTesting beam correction...")
    beam_ell, beam_prof = load_beam_profile(None)  # Using default Gaussian
    # Just test the function works
    beam_test = apply_beam_correction(filtered_map, beam_ell, beam_prof, 'deconvolve')
    print(f"Beam correction applied successfully")

    # Save visualization
    output_dir = "data/output/test_plots"
    os.makedirs(output_dir, exist_ok=True)

    # Downsample for plotting
    ds = 20

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    axes[0].imshow(temp_map[::ds, ::ds], cmap='RdBu_r', aspect='auto')
    axes[0].set_title('Original Map')
    axes[0].set_xlabel('X pixel')
    axes[0].set_ylabel('Y pixel')

    axes[1].imshow(masked_map[::ds, ::ds], cmap='RdBu_r', aspect='auto')
    axes[1].set_title('Masked Map')
    axes[1].set_xlabel('X pixel')

    axes[2].imshow(filtered_map[::ds, ::ds], cmap='RdBu_r', aspect='auto')
    axes[2].set_title('Filtered Map (ell > 100)')
    axes[2].set_xlabel('X pixel')

    plt.tight_layout()
    output_file = os.path.join(output_dir, 'preprocessing_test.png')
    plt.savefig(output_file, dpi=100, bbox_inches='tight')
    plt.close()
    print(f"\nSaved visualization to: {output_file}")

    return True


if __name__ == "__main__":
    success = test_preprocessing()
    if success:
        print("\nTest PASSED: Preprocessing functions working")
    else:
        print("\nTest FAILED")
