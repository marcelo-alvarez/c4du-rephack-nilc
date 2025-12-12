# Validation Instructions for NILC Pipeline Outputs

This document specifies how to validate the final map and data products from the NILC pipeline. Validation should be performed after completing Phase 5 (Output and validation) and before using the maps for scientific analysis.

## Output Products to Validate

The pipeline produces:
- **CMB component maps**: Temperature (T) and polarization (Q, U) in both CAR and HEALPix formats
- **kSZ component maps**: Temperature (T) and polarization (Q, U) in both CAR and HEALPix formats

## Validation Checklist

### 1. Basic Data Quality

For each output map (CMB and kSZ, both formats):
- No NaN or Inf values present
- Physically reasonable pixel value ranges (compare to input map ranges)
- CAR and HEALPix versions of the same map agree within numerical precision
- Temperature and polarization maps have consistent coverage masks

### 2. Power Spectrum Validation

Compute angular power spectra and verify:
- **CMB TT**: Shows characteristic acoustic peak structure (peaks around ℓ ~200, 500, 800)
- **kSZ TT**: Featureless and subdominant to CMB on large scales (ℓ < 3000)
- **CMB EE**: Shows reionization bump at low ℓ and acoustic peaks
- **kSZ EE**: Subdominant (kSZ is primarily temperature-only)
- High-ℓ noise levels are reasonable (compare to input map noise)

### 3. Visual Inspection

Generate full-sky visualizations and check:
- CMB temperature map shows large-scale structure (no obvious artifacts or striping)
- kSZ temperature map shows localized structure (not pure noise)
- No discontinuities or systematic biases
- CAR and HEALPix versions appear consistent (accounting for projection differences)

### 4. Benchmark Comparison

**Note**: Benchmark (correct output) maps will be provided for comparison. When available, compare pipeline outputs to benchmarks:

- **Map-level comparison**:
  - Compute pixel-by-pixel difference maps (output - benchmark) for CMB and kSZ components
  - Residuals should be consistent with noise levels (not systematic biases)
  - Compute correlation coefficient between output and benchmark maps (should be > 0.95 for CMB on large scales)

- **Power spectrum comparison**:
  - Compute power spectra for both output and benchmark maps
  - Compare CMB TT, EE, TE spectra: output should match benchmark within expected noise variance
  - Compare kSZ TT spectra: output should match benchmark, especially on scales where kSZ is significant

- **Statistics comparison**:
  - Compare pixel value distributions (histograms) between output and benchmark
  - Compare map statistics (mean, std, min, max) per component
  - Verify that differences are consistent with expected numerical precision and noise

- **Component separation quality**:
  - Compare CMB/kSZ ratio in power spectra between output and benchmark
  - Verify that component separation quality (e.g., CMB dominance on large scales) matches benchmark

Example benchmark comparison commands:
```bash
# Compare output to benchmark maps
python scripts/compare_to_benchmark.py \
    --output-cmb path/to/output_cmb.fits \
    --benchmark-cmb path/to/benchmark_cmb.fits \
    --output-ksz path/to/output_ksz.fits \
    --benchmark-ksz path/to/benchmark_ksz.fits \
    --format [car|healpix]
```

### 5. Input Reconstruction Check

- Verify that combining CMB and kSZ components (with frequency-dependent scaling) reconstructs input frequency maps within noise levels
- Check residuals between input and reconstructed maps (should be consistent with noise, not systematic biases)

## Validation Tools

- **Power spectra**: `healpy.anafast()` for HEALPix maps, equivalent for CAR maps
- **Visualization**: `matplotlib`, `healpy.mollview()`, or `pixell` tools
- **FITS I/O**: `healpy.read_map()` and `pixell.enmap.read_map()`

## Success Criteria

Validation passes if:
1. All basic data quality checks pass
2. Power spectra show physically reasonable structure (CMB peaks visible, kSZ subdominant on large scales)
3. Visual inspection reveals no obvious artifacts
4. **Benchmark comparison** (when available): Output maps match benchmarks within expected noise/precision limits
5. Input maps can be reconstructed from outputs within noise levels

If validation fails, document whether failures indicate implementation bugs, expected method limitations, or need for parameter adjustments.