# Manager Status

## Project overview
Implement a NILC (Needlet Internal Linear Combination) pipeline for ACT cosmology data processing, following arXiv:2307.01258 Section III. The pipeline will process two ACT frequency maps (90 GHz and 150 GHz) in CAR projection to separate CMB and kSZ components.

## Current state
- **Status**: Project initialization
- **Last updated**: 2025-12-12
- **Repository**: Not yet initialized

## Milestones

### Phase 1: Project setup
- [ ] Initialize git repository
- [ ] Set up Python project structure (directories, requirements)
- [ ] Create placeholder test framework

### Phase 2: Preprocessing implementation (Section III-A)
- [ ] Implement inpainting for high-significance detections
- [ ] Implement Fourier-domain filtering for scan-synchronous modes
- [ ] Implement beam difference corrections (harmonic space)
- [ ] Implement color corrections for component spectral response

### Phase 3: Needlet decomposition (Section III-B)
- [ ] Implement axisymmetric needlet kernel (Eq. 2 from paper)
- [ ] Set up needlet scale configuration (26 ell_peak values)
- [ ] Implement needlet decomposition for CAR maps
- [ ] Validate needlet localization in ell-space

### Phase 4: ILC weights and component separation (Section III-C)
- [ ] Implement in-band covariance matrix computation
- [ ] Implement harmonic-space bias exclusion (large scales)
- [ ] Implement real-space "donut" smoothing (small scales)
- [ ] Compute ILC weights per needlet scale
- [ ] Apply weights to separate CMB and kSZ components

### Phase 5: Output and validation
- [ ] Recombine needlet bands into full-sky maps
- [ ] Export CAR FITS files (via pixell)
- [ ] Export HEALPix FITS files (via healpy)
- [ ] Implement basic validation (power spectra, visual checks)

## Outstanding tasks
None yet; project is in initialization phase.

## Recent completions
- 2025-12-12: Created `.agent/` directory structure
- 2025-12-12: Moved instruction files to `.agent/`
- 2025-12-12: Initialized manager-status.md

## Notes and decisions
- Input: Two CAR-projected maps (90 GHz and 150 GHz) with temperature and polarization
- Two frequencies only; all covariance matrices and ILC weights computed from this pair
- Component separation targets: CMB (T+P) and kSZ
- Dependencies: numpy, scipy, matplotlib, pixell, healpy
- Output both CAR (pixell.enmap.write_map) and HEALPix (healpy.write_map) formats
