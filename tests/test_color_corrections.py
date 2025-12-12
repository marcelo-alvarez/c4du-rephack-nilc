"""
Test color corrections for ACT maps.

Loads 90 GHz and 150 GHz maps, applies color corrections, and visualizes results.
"""
import numpy as np
import matplotlib.pyplot as plt
from pixell import enmap
import sys
import os

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from nilc.preprocessing.color_corrections import (
    load_passband,
    compute_color_correction,
    apply_color_correction
)


def test_color_corrections():
    """Test color correction functionality."""
    print("=" * 70)
    print("Testing Color Corrections for ACT Maps")
    print("=" * 70)

    # Paths
    map_dir = '/scratch/jiaqu/hack_data/maps/'
    mask_file = '/scratch/jiaqu/hack_data/masks/ilc_footprint_mask.fits'
    passband_dir = '/home/jiaqu/NILC/data/'

    # Map files
    map_90_file = os.path.join(
        map_dir,
        'act_dr6.02_std_AA_night_pa5_f090_4way_coadd_map_srcfree.fits'
    )
    map_150_file = os.path.join(
        map_dir,
        'act_dr6.02_std_AA_night_pa4_f150_4way_coadd_map_srcfree.fits'
    )

    # Passband files (truncated versions with normalized responses)
    passband_90_file = os.path.join(passband_dir, 'PA5_avg_passband_090_wErr_trunc.txt')
    passband_150_file = os.path.join(passband_dir, 'PA5_avg_passband_150_wErr_trunc.txt')

    print("\n1. Loading passbands...")
    print("-" * 70)

    # Load and verify passbands
    freq_90, response_90 = load_passband(passband_90_file)
    freq_150, response_150 = load_passband(passband_150_file)

    print(f"90 GHz passband:")
    print(f"  Frequency range: {freq_90.min():.2f} - {freq_90.max():.2f} GHz")
    print(f"  Number of points: {len(freq_90)}")
    print(f"  Peak response frequency: {freq_90[np.argmax(response_90)]:.2f} GHz")

    print(f"\n150 GHz passband:")
    print(f"  Frequency range: {freq_150.min():.2f} - {freq_150.max():.2f} GHz")
    print(f"  Number of points: {len(freq_150)}")
    print(f"  Peak response frequency: {freq_150[np.argmax(response_150)]:.2f} GHz")

    print("\n2. Computing color corrections for CMB...")
    print("-" * 70)

    # Compute color corrections for CMB
    corr_90, eff_freq_90 = compute_color_correction(
        freq_90, response_90, 'CMB'
    )
    corr_150, eff_freq_150 = compute_color_correction(
        freq_150, response_150, 'CMB'
    )

    print(f"90 GHz band:")
    print(f"  Effective frequency: {eff_freq_90:.2f} GHz")
    print(f"  Color correction factor: {corr_90:.6f}")

    print(f"\n150 GHz band:")
    print(f"  Effective frequency: {eff_freq_150:.2f} GHz")
    print(f"  Color correction factor: {corr_150:.6f}")

    print("\n3. Loading ACT maps...")
    print("-" * 70)

    # Load maps (temperature component only, [0] index)
    print(f"Loading 90 GHz map: {os.path.basename(map_90_file)}")
    map_90 = enmap.read_map(map_90_file, sel=np.s_[0, :, :])
    print(f"  Shape: {map_90.shape}")
    print(f"  WCS: {map_90.wcs}")

    print(f"\nLoading 150 GHz map: {os.path.basename(map_150_file)}")
    map_150 = enmap.read_map(map_150_file, sel=np.s_[0, :, :])
    print(f"  Shape: {map_150.shape}")
    print(f"  WCS: {map_150.wcs}")

    print("\n4. Loading ILC footprint mask...")
    print("-" * 70)
    mask = enmap.read_map(mask_file)
    print(f"  Mask shape: {mask.shape}")
    print(f"  Mask coverage: {np.sum(mask > 0.5) / mask.size * 100:.1f}%")

    # Reproject mask to map geometry if needed
    if mask.shape != map_90.shape or not np.allclose(mask.wcs.wcs.cdelt, map_90.wcs.wcs.cdelt):
        print("  Reprojecting mask to map geometry...")
        mask = enmap.project(mask, map_90.shape, map_90.wcs, order=1)

    # Apply mask
    map_90_masked = map_90 * mask
    map_150_masked = map_150 * mask

    print("\n5. Map statistics before color correction...")
    print("-" * 70)
    valid_90 = map_90_masked[mask > 0.5]
    valid_150 = map_150_masked[mask > 0.5]

    print(f"90 GHz map (masked region):")
    print(f"  Mean: {np.mean(valid_90):.6e} K")
    print(f"  Std:  {np.std(valid_90):.6e} K")
    print(f"  Min:  {np.min(valid_90):.6e} K")
    print(f"  Max:  {np.max(valid_90):.6e} K")

    print(f"\n150 GHz map (masked region):")
    print(f"  Mean: {np.mean(valid_150):.6e} K")
    print(f"  Std:  {np.std(valid_150):.6e} K")
    print(f"  Min:  {np.min(valid_150):.6e} K")
    print(f"  Max:  {np.max(valid_150):.6e} K")

    print("\n6. Applying color corrections...")
    print("-" * 70)

    # Apply color corrections
    (map_90_corr, map_150_corr,
     correction_90, correction_150,
     eff_f90, eff_f150) = apply_color_correction(
        map_90_masked, map_150_masked, component='CMB'
    )

    print(f"Applied corrections:")
    print(f"  90 GHz:  factor = {correction_90:.6f}, eff_freq = {eff_f90:.2f} GHz")
    print(f"  150 GHz: factor = {correction_150:.6f}, eff_freq = {eff_f150:.2f} GHz")

    print("\n7. Map statistics after color correction...")
    print("-" * 70)
    valid_90_corr = map_90_corr[mask > 0.5]
    valid_150_corr = map_150_corr[mask > 0.5]

    print(f"90 GHz map (color-corrected, masked):")
    print(f"  Mean: {np.mean(valid_90_corr):.6e} K")
    print(f"  Std:  {np.std(valid_90_corr):.6e} K")
    print(f"  Min:  {np.min(valid_90_corr):.6e} K")
    print(f"  Max:  {np.max(valid_90_corr):.6e} K")

    print(f"\n150 GHz map (color-corrected, masked):")
    print(f"  Mean: {np.mean(valid_150_corr):.6e} K")
    print(f"  Std:  {np.std(valid_150_corr):.6e} K")
    print(f"  Min:  {np.min(valid_150_corr):.6e} K")
    print(f"  Max:  {np.max(valid_150_corr):.6e} K")

    print("\n8. Creating visualization...")
    print("-" * 70)

    # Downsample maps for visualization
    downsample_factor = 20
    print(f"Downsampling maps by factor {downsample_factor} for visualization")

    map_90_down = map_90_masked[::downsample_factor, ::downsample_factor]
    map_150_down = map_150_masked[::downsample_factor, ::downsample_factor]
    map_90_corr_down = map_90_corr[::downsample_factor, ::downsample_factor]
    map_150_corr_down = map_150_corr[::downsample_factor, ::downsample_factor]
    mask_down = mask[::downsample_factor, ::downsample_factor]

    # Create figure with subplots
    fig = plt.figure(figsize=(16, 10))

    # Subplot 1: Passbands
    ax1 = plt.subplot(2, 3, 1)
    ax1.plot(freq_90, response_90, 'b-', label='90 GHz', linewidth=2)
    ax1.plot(freq_150, response_150, 'r-', label='150 GHz', linewidth=2)
    ax1.axvline(eff_freq_90, color='b', linestyle='--', alpha=0.5,
                label=f'90 GHz eff: {eff_freq_90:.1f} GHz')
    ax1.axvline(eff_freq_150, color='r', linestyle='--', alpha=0.5,
                label=f'150 GHz eff: {eff_freq_150:.1f} GHz')
    ax1.set_xlabel('Frequency [GHz]')
    ax1.set_ylabel('Normalized Response')
    ax1.set_title('ACT PA5 Passbands')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Common colorbar limits
    vmin = min(np.percentile(map_90_down[mask_down > 0.5], 1),
               np.percentile(map_150_down[mask_down > 0.5], 1))
    vmax = max(np.percentile(map_90_down[mask_down > 0.5], 99),
               np.percentile(map_150_down[mask_down > 0.5], 99))

    # Subplot 2: 90 GHz original
    ax2 = plt.subplot(2, 3, 2)
    im2 = ax2.imshow(map_90_down, cmap='RdBu_r', vmin=vmin, vmax=vmax,
                     origin='lower', aspect='auto')
    ax2.set_title('90 GHz Map (Original)')
    ax2.set_xlabel('X pixel')
    ax2.set_ylabel('Y pixel')
    plt.colorbar(im2, ax=ax2, label='Temperature [K]', fraction=0.046)

    # Subplot 3: 150 GHz original
    ax3 = plt.subplot(2, 3, 3)
    im3 = ax3.imshow(map_150_down, cmap='RdBu_r', vmin=vmin, vmax=vmax,
                     origin='lower', aspect='auto')
    ax3.set_title('150 GHz Map (Original)')
    ax3.set_xlabel('X pixel')
    ax3.set_ylabel('Y pixel')
    plt.colorbar(im3, ax=ax3, label='Temperature [K]', fraction=0.046)

    # Subplot 4: Color correction factors
    ax4 = plt.subplot(2, 3, 4)
    factors = [correction_90, correction_150]
    labels = ['90 GHz', '150 GHz']
    colors = ['blue', 'red']
    bars = ax4.bar(labels, factors, color=colors, alpha=0.7)
    ax4.set_ylabel('Correction Factor')
    ax4.set_title('Color Correction Factors (CMB)')
    ax4.grid(True, alpha=0.3, axis='y')
    # Add value labels on bars
    for i, (bar, val) in enumerate(zip(bars, factors)):
        ax4.text(bar.get_x() + bar.get_width()/2, val + 0.0002,
                f'{val:.5f}', ha='center', va='bottom')

    # Subplot 5: 90 GHz corrected
    ax5 = plt.subplot(2, 3, 5)
    im5 = ax5.imshow(map_90_corr_down, cmap='RdBu_r', vmin=vmin, vmax=vmax,
                     origin='lower', aspect='auto')
    ax5.set_title('90 GHz Map (Color-Corrected)')
    ax5.set_xlabel('X pixel')
    ax5.set_ylabel('Y pixel')
    plt.colorbar(im5, ax=ax5, label='Temperature [K]', fraction=0.046)

    # Subplot 6: 150 GHz corrected
    ax6 = plt.subplot(2, 3, 6)
    im6 = ax6.imshow(map_150_corr_down, cmap='RdBu_r', vmin=vmin, vmax=vmax,
                     origin='lower', aspect='auto')
    ax6.set_title('150 GHz Map (Color-Corrected)')
    ax6.set_xlabel('X pixel')
    ax6.set_ylabel('Y pixel')
    plt.colorbar(im6, ax=ax6, label='Temperature [K]', fraction=0.046)

    plt.tight_layout()

    # Save figure
    output_file = 'color_correction_test_results.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Saved visualization to: {output_file}")

    print("\n" + "=" * 70)
    print("Color correction test completed successfully!")
    print("=" * 70)

    return True


if __name__ == '__main__':
    success = test_color_corrections()
    sys.exit(0 if success else 1)
