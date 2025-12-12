"""Utilities for generating dummy benchmark data."""
import numpy as np
import os
from pixell import enmap
from astropy import wcs
from nilc.io.map_io import load_act_map


def generate_dummy_benchmark(output_path=None, shape=(10320, 43200), seed=42, benchmark_path=None):
    """
    Generate a dummy CMB temperature map in CAR projection matching the benchmark format.
    
    Parameters
    ----------
    output_path : str, optional
        Path to save the dummy benchmark. If None, saves to data/benchmark/dummy_benchmark_T.fits
    shape : tuple, default=(10320, 43200)
        Shape of the output map (ny, nx)
    seed : int, default=42
        Random seed for reproducibility
    benchmark_path : str, optional
        Path to the real benchmark map. If provided, will use the same mask (zero/NaN pixels) as the benchmark.
        If None, uses: /secret/path/to/target/output/target_T.fits
        
    Returns
    -------
    str
        Path to the saved dummy benchmark file
    """
    print(f"Starting dummy benchmark generation...")
    print(f"  Shape: {shape}")
    print(f"  Seed: {seed}")
    
    np.random.seed(seed)
    
    # Create a map with the specified shape
    # Generate simple random noise
    ny, nx = shape
    print(f"Generating random noise for {ny}x{nx} map...")
    
    # Load benchmark first to get statistics and mask
    if benchmark_path is None:
        benchmark_path = "/secret/path/to/target/output/target_T.fits"
    
    print(f"Loading benchmark for mask and statistics from: {benchmark_path}")
    try:
        benchmark_map = load_act_map(benchmark_path)
        
        # Load the footprint mask
        mask_path = "/secret/path/to/target/output/ilc_footprint_mask.fits"
        print(f"Loading footprint mask from: {mask_path}")
        footprint_mask = load_act_map(mask_path)
        
        # Mask is 1 for valid pixels, 0 for masked pixels
        valid_mask = (footprint_mask == 1)
        masked_mask = (footprint_mask == 0)
        print(f"  Valid pixels (mask==1): {np.sum(valid_mask)} ({100*np.sum(valid_mask)/footprint_mask.size:.1f}%)")
        print(f"  Masked pixels (mask==0): {np.sum(masked_mask)} ({100*np.sum(masked_mask)/footprint_mask.size:.1f}%)")
        
        # Get benchmark statistics on valid (unmasked) pixels only
        benchmark_valid = benchmark_map[valid_mask]
        benchmark_mean = np.mean(benchmark_valid)
        benchmark_std = np.std(benchmark_valid)
        benchmark_min = np.min(benchmark_valid)
        benchmark_max = np.max(benchmark_valid)
        print(f"  Benchmark stats on valid pixels: mean={benchmark_mean:.3e}, std={benchmark_std:.3e}, "
              f"min={benchmark_min:.3e}, max={benchmark_max:.3e}")
        
        # Get the value to use for masked pixels (median of masked pixels in benchmark)
        benchmark_masked_values = benchmark_map[masked_mask]
        masked_value = np.median(benchmark_masked_values)
        print(f"  Masked pixel value (median): {masked_value:.6e}")
        
    except Exception as e:
        print(f"  Warning: Could not load benchmark or mask: {e}")
        print("  Using default statistics.")
        valid_mask = np.ones((ny, nx), dtype=bool)
        masked_mask = np.zeros((ny, nx), dtype=bool)
        benchmark_mean = 0.0
        benchmark_std = 77.66  # Default from real benchmark
        benchmark_min = -1355.0
        benchmark_max = 2857.0
        masked_value = -1.284839  # Default median of masked pixels
    
    # Generate random noise (use float32 to match original benchmark format)
    temp_map = np.random.normal(0, 1, size=(ny, nx)).astype(np.float32)
    print("  Random noise generated.")
    
    # Scale to match benchmark statistics on valid pixels only
    print("Scaling to match benchmark statistics...")
    temp_map = temp_map - np.mean(temp_map)  # Center at zero
    temp_map = temp_map * (benchmark_std / np.std(temp_map))  # Scale to match std
    temp_map = temp_map + benchmark_mean  # Shift to match mean
    print("  Scaling complete.")
    
    # Apply mask to dummy data: set masked pixels to the masked value
    temp_map = temp_map * valid_mask.astype(np.float32) + masked_value * masked_mask.astype(np.float32)
    print("  Mask applied to dummy data (masked pixels set to constant value).")
    
    # Create WCS matching the benchmark format
    # WCS from benchmark: car:{cdelt:[-0.008333,0.008333],crval:[0,0],crpix:[21601.00,7560.50]}
    print("Creating WCS...")
    w = wcs.WCS(naxis=2)
    w.wcs.cdelt = [-0.008333, 0.008333]  # degrees per pixel
    w.wcs.crval = [0.0, 0.0]  # reference pixel value (degrees)
    w.wcs.crpix = [nx/2.0 + 0.5, ny/2.0 + 0.5]  # reference pixel (1-indexed)
    w.wcs.ctype = ['RA---CAR', 'DEC--CAR']
    w.wcs.cunit = ['deg', 'deg']
    print("  WCS created.")
    
    # Create enmap with proper WCS
    print("Creating enmap...")
    dummy_map = enmap.ndmap(temp_map, wcs=w)
    print("  Enmap created.")
    
    # Set output path
    if output_path is None:
        # Get project root (assuming we're in src/nilc/utils/)
        project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
        output_path = os.path.join(project_root, "data", "benchmark", "dummy_benchmark_T.fits")
    
    print(f"Output path: {output_path}")
    
    # Ensure output directory exists
    print("Creating output directory...")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    print("  Directory ready.")
    
    # Write to FITS file
    print("Writing FITS file (this may take a while for large maps)...")
    enmap.write_map(output_path, dummy_map)
    print("  FITS file written.")
    
    print(f"Generated dummy benchmark map:")
    print(f"  Shape: {dummy_map.shape}")
    print(f"  WCS: {dummy_map.wcs}")
    print(f"  Statistics: min={np.nanmin(dummy_map):.3e}, max={np.nanmax(dummy_map):.3e}, "
          f"mean={np.nanmean(dummy_map):.3e}, std={np.nanstd(dummy_map):.3e}")
    print(f"  Saved to: {output_path}")
    
    return output_path


if __name__ == "__main__":
    output_path = generate_dummy_benchmark()
    print(f"\nDummy benchmark created at: {output_path}")

