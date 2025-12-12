"""Test script for reading and visualizing ACT maps."""
import os
import numpy as np
import matplotlib.pyplot as plt
from pixell import enmap, enplot
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from nilc.io.map_io import load_act_map, get_stokes_component

def test_read_and_plot_map():
    """Test reading and plotting a single ACT source-free map."""
    # Define map path
    map_dir = "/scratch/jiaqu/hack_data/maps/"
    map_file = "act_dr6.02_std_AA_night_pa4_f220_4way_coadd_map_srcfree.fits"
    map_path = os.path.join(map_dir, map_file)

    print(f"Loading map from: {map_path}")

    # Check if file exists
    if not os.path.exists(map_path):
        print(f"ERROR: Map file not found at {map_path}")
        return False

    # Load the map
    map_data = load_act_map(map_path)
    print(f"Map shape: {map_data.shape}")
    print(f"Map dtype: {map_data.dtype}")
    print(f"Map WCS: {map_data.wcs}")

    # Extract temperature (Stokes I) if multi-component
    if map_data.ndim == 3:
        temp_map = get_stokes_component(map_data, 'I')
        print(f"Extracted Stokes I, shape: {temp_map.shape}")
    else:
        temp_map = map_data

    # Print statistics
    print(f"Temperature map statistics:")
    print(f"  Min: {np.min(temp_map):.3e}")
    print(f"  Max: {np.max(temp_map):.3e}")
    print(f"  Mean: {np.mean(temp_map):.3e}")
    print(f"  Std: {np.std(temp_map):.3e}")

    # Create output directory
    output_dir = "data/output/test_plots"
    os.makedirs(output_dir, exist_ok=True)

    # Plot and save using matplotlib
    plt.figure(figsize=(12, 8))
    plt.imshow(temp_map, cmap='RdBu_r', aspect='auto')
    plt.colorbar(label='Temperature [K or uK]')
    plt.title(f'ACT 220 GHz Map (Source-Free)\n{map_file}')
    plt.xlabel('X pixel')
    plt.ylabel('Y pixel')
    output_file = os.path.join(output_dir, 'act_f220_srcfree_test.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved plot to: {output_file}")

    # Also save using enplot for comparison
    output_file_enplot = os.path.join(output_dir, 'act_f220_srcfree_test_enplot.png')
    enplot.write(output_file_enplot, enplot.plot(temp_map, colorbar=True))
    print(f"Saved enplot version to: {output_file_enplot}")

    return True

if __name__ == "__main__":
    success = test_read_and_plot_map()
    if success:
        print("\nTest PASSED: Successfully read and plotted ACT map")
    else:
        print("\nTest FAILED")
