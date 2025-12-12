# Manager Status

## Project overview
Implement a NILC (Needlet Internal Linear Combination) pipeline for ACT cosmology data processing, following arXiv:2307.01258 Section III. The pipeline will process two ACT frequency maps (90 GHz and 150 GHz) in CAR projection to separate CMB and kSZ components.

## Current state
- **Status**: Phase 2 complete - ready for Phase 3 (needlet decomposition)
- **Last updated**: 2025-12-12
- **Repository**: https://github.com/marcelo-alvarez/c4du-rephack-nilc
- **Next worker**: Worker #5 - implement needlet decomposition (Section III-B)

## Milestones

### Phase 1: Project setup
- [x] Initialize git repository
- [x] Set up Python project structure (directories, requirements)
- [x] Create placeholder test framework
- [x] Push project structure to GitHub

### Phase 2: Preprocessing implementation (Section III-A)
- [x] Implement basic map I/O and verify data access
- [ ] Implement inpainting for high-significance detections (deferred)
- [x] Implement Fourier-domain filtering for scan-synchronous modes
- [x] Implement beam difference corrections (harmonic space)
- [x] Implement color corrections for component spectral response

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
- Implement needlet decomposition (Phase 3)
- Implement ILC component separation (Phase 4)
- Implement output and validation (Phase 5)

## Recent completions
- 2025-12-12: Created `.agent/` directory structure and initialized manager-status.md
- 2025-12-12: Worker #1 completed - initialized git repo, created Python package structure (src/nilc with 5 modules), test framework, requirements.txt, setup.py, README.md
- 2025-12-12: Worker #1c completed - successfully pushed all project files to GitHub (17 files total)
- 2025-12-12: Worker #1d completed - updated and pushed manager-status.md to reflect Phase 1 completion
- 2025-12-12: Worker #2 completed - implemented map I/O utilities and successfully tested reading ACT 220 GHz source-free map (3 Stokes components, 10320x43200 CAR projection)
- 2025-12-12: Worker #2b completed - fixed visualization by downsampling maps (20x factor), reduced enplot output from 799 MB to 584 KB
- 2025-12-12: Worker #3 completed - implemented preprocessing functions (Fourier filtering, beam corrections), tested with ACT mask (231M pixels)
- 2025-12-12: Worker #4 completed - implemented color corrections using actual PA5 passbands; effective frequencies: 94.96 GHz (90 band), 147.02 GHz (150 band); correction factors ~1.002-1.005

## Notes and decisions
- Input: Two CAR-projected maps (90 GHz and 150 GHz) with temperature and polarization
- Two frequencies only; all covariance matrices and ILC weights computed from this pair
- Component separation targets: CMB (T+P) and kSZ
- Dependencies: numpy, scipy, matplotlib, pixell, healpy
- Output both CAR (pixell.enmap.write_map) and HEALPix (healpy.write_map) formats
