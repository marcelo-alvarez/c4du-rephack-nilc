# Project Context

## Scientific background
- The guiding reference is "The Atacama Cosmology Telescope: High-resolution component-separated maps across one-third of the sky" (arXiv:2307.01258). Section III of the paper describes the NILC pipeline we must emulate: preprocessing, needlet decomposition, and ILC-based component separation with explicit ILC bias mitigation.
- Section III-A details preprocessing steps for ACT-like maps: inpaint around high-significance detections, apply Fourier-domain filtering to remove scan-synchronous modes, account for beam differences with harmonic corrections, and apply color corrections tied to each component's spectral response to match beams by scale.
- Section III-B defines the needlet frame: axisymmetric kernel described by Eq. (2) with ell_peak values [0, 100, 200, 300, 400, 600, 800, 1000, 1250, 1400, 1800, 2200, 2800, 3400, 4000, 4600, 5200, 6000, 7000, 8000, 9000, 10000, 11000, 13000, 15000, 17000]; each scale localizes modes in ell space, letting us adaptively weight spatial regions.
- Section III-C shows that ILC weights come from in-band covariance matrices and need specialized bias control: large scales use a harmonic-space exclusion strategy and small scales rely on a real-space "donut" smoothing so the weights do not reuse the same modes they are applied to; this reduces the classic ILC bias caused by chance correlation of the signal and noise.

## Immediate implementation goal
We have two ACT frequency maps (90 GHz and 150 GHz) in CAR projection with both temperature and polarization. The work is to:
1. Preprocess the input maps following Section III-A (inpainting, Fourier filtering, beam corrections, color corrections).
2. Decompose the maps into needlet scales per Section III-B and compute ILC weights for each scale.
3. Perform component separation (NILC) on each needlet band to reconstruct two components: (a) CMB (temperature and polarization), and (b) kinetic Sunyaev-Zel'dovich (kSZ).
4. Mitigate ILC bias using the hybrid harmonic/real-space exclusion strategy described in Section III-C.
5. Recombine the needlet outputs and save both CAR (via `pixell`) and HEALPix (via `healpy`) format maps.

## Constraints & assumptions
- Input data are CAR-projected temperature and polarization maps at 90 GHz and 150 GHz.
- Two frequencies only; all covariance matrices, ILC weights, and component separation use this pair.
- Dependencies available: `numpy`, `scipy`, `matplotlib`, `pixell`, `healpy`.
- Processing must respect ACT beam profiles with scale-dependent color corrections for frequency-dependent beams.
- Frequency response functions implemented for CMB (f=1), kSZ (f=1), tSZ (f(x) = x(e^x+1)/(e^x-1) - 4), and CIB (modified blackbody: ν^β B_ν(T_dust) with T_dust=20K, β=1.5, normalized at 353 GHz).

## Data locations
- **Maps**: `/global/cfs/cdirs/act/data/act_dr6/dr6.02/maps/published/`
  - Use non-HEALPix format (CAR projection)
  - Use source-free (`srcfree`) versions: `act_dr4dr6_coadd_AA_daynight_f090_map_srcfree.fits`, `act_dr4dr6_coadd_AA_daynight_f150_map_srcfree.fits`
  - Use AA (Atacama Array) maps
  - Frequencies: 90 GHz and 150 GHz
- **Beams**: `/global/cfs/cdirs/act/data/act_dr6/dr6.02/beams/main_beams/nominal/`
  - PA5 90 GHz: `coadd_pa5_f090_night_beam_tform_jitter_cmb.txt`
  - PA5 150 GHz: `coadd_pa5_f150_night_beam_tform_jitter_cmb.txt`
- **Mask**: `/global/cfs/cdirs/act/data/act_dr6/dr6.02/nilc/published/ilc_footprint_mask.fits` (analysis footprint)
- **Data documentation**: https://lambda.gsfc.nasa.gov/product/act/act_dr6.02/act_dr6.02_maps_info.html

## Module structure
- `src/nilc/needlets/` - Needlet decomposition module
  - `decomposition.py` - Axisymmetric needlet kernels (Eq. 2 of arXiv:2307.01258), decomposition and recombination
  - `pipeline.py` - End-to-end NILC pipeline for CMB extraction
- `src/nilc/ilc/` - ILC weight computation
  - `weights.py` - Constrained ILC weight computation with Lagrange multipliers
  - `covariance.py` - Bias-mitigated covariance matrix computation
  - `frequency_response.py` - CMB, kSZ, tSZ, CIB frequency response functions
- `src/nilc/preprocessing/` - Map preprocessing
  - `filters.py` - Fourier-domain filtering
  - `beams.py` - Beam correction utilities
  - `color_corrections.py` - Color correction for component spectral response
- `scripts/run_nilc.py` - Runner script for the NILC pipeline

## Running the pipeline
```bash
cd /global/cfs/cdirs/cosmosim/slac/malvarez/hackathon/c4du-rephack-nilc
python scripts/run_nilc.py --output-dir ./output --component I --target cmb
```

Options:
- `--component I|Q|U` - Stokes component to process (default: I)
- `--target cmb|ksz` - Target component to extract (default: cmb)
- `--lmin INT` - Minimum ell for high-pass filter (default: 100)
- `--save-intermediate` - Save intermediate needlet scale maps

## Operational notes for workers
- Anchor any worker prompt with: `Read agent/worker-instructions.md and agent/project-context.md`.
- Provide clear command lists, prioritize raster operations via `pixell` transforms, and include steps to write both CAR and HEALPix maps.
- Output artifacts should include CAR FITS files (e.g., via `pixell.enmap.write_map`) and HEALPix FITS (e.g., via `healpy.write_map`) for the CMB and kSZ components.
- Include any validation (e.g., visual spot checks or simple power spectrum comparisons) only if the prompt explicitly requests it.
