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

## Data locations
- **Maps**: `/scratch/jiaqu/hack_data/maps/`
  - Use non-HEALPix format (CAR projection)
  - Use source-free (`srcfree`) versions
  - Use AA (Atacama Array) maps
  - Frequencies: 90 GHz and 150 GHz
- **Beams**: `/scratch/jiaqu/hack_data/maps/beams/main_beams/`
- **Mask**: `/scratch/jiaqu/hack_data/masks/ilc_footprint_mask.fits` (analysis footprint)
- **Passbands**: `/home/jiaqu/NILC/data/ACT_ancillary/`
  - PA5 90 GHz: `PA5_avg_passband_90_wErr.txt` (3 columns: frequency [GHz], response, error)
  - PA5 150 GHz: `PA5_avg_passband_150_wErr.txt` (3 columns: frequency [GHz], response, error)
  - Passbands are needed for color corrections to account for component spectral response
- **Data documentation**: https://lambda.gsfc.nasa.gov/product/act/act_dr6.02/act_dr6.02_maps_info.html

## Operational notes for workers
- Anchor any worker prompt with: `Read .agent/worker-instructions.md and .agent/project-context.md`.
- Provide clear command lists, prioritize raster operations via `pixell` transforms, and include steps to write both CAR and HEALPix maps.
- Output artifacts should include CAR FITS files (e.g., via `pixell.enmap.write_map`) and HEALPix FITS (e.g., via `healpy.write_map`) for the CMB and kSZ components.
- Include any validation (e.g., visual spot checks or simple power spectrum comparisons) only if the prompt explicitly requests it.
