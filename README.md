# ACT NILC Pipeline

NILC (Needlet Internal Linear Combination) pipeline implementation for Atacama Cosmology Telescope (ACT) data. This package performs component separation to extract CMB and kinetic Sunyaev-Zel'dovich (kSZ) signals from multi-frequency microwave maps.

## Overview

This implementation follows the methodology described in Section III of "The Atacama Cosmology Telescope: High-resolution component-separated maps across one-third of the sky" ([arXiv:2307.01258](https://arxiv.org/abs/2307.01258)).

The pipeline consists of three main stages:

1. **Preprocessing** (Section III-A): Inpainting, Fourier-domain filtering, beam corrections, and color corrections
2. **Needlet Decomposition** (Section III-B): Multi-scale decomposition using axisymmetric needlet kernels
3. **ILC Component Separation** (Section III-C): Internal Linear Combination with bias mitigation strategies

### Key Features

- Processes CAR-projected temperature and polarization maps
- Multi-frequency analysis (90 GHz and 150 GHz)
- Needlet decomposition with 25+ scales
- ILC bias mitigation using hybrid harmonic/real-space exclusion
- Outputs in both CAR (pixell) and HEALPix (healpy) formats

## Installation

### Prerequisites

- Python 3.8 or higher
- Git

### Setup

1. Clone the repository:
```bash
git clone git@github.com:marcelo-alvarez/c4du-rephack-nilc.git
cd c4du-rephack-nilc
```

2. Create a virtual environment (recommended):
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install the package and dependencies:
```bash
pip install -r requirements.txt
pip install -e .
```

## Project Structure

```
.
├── src/nilc/              # Main package
│   ├── preprocessing/     # Map preprocessing utilities
│   ├── needlets/          # Needlet decomposition and synthesis
│   ├── ilc/               # ILC weights and component separation
│   ├── io/                # Input/output for CAR and HEALPix maps
│   └── utils/             # General utilities
├── tests/                 # Unit and integration tests
│   ├── unit/              # Unit tests
│   └── integration/       # Integration tests
├── data/                  # Data directory (excluded from git)
│   ├── input/             # Input frequency maps
│   └── output/            # Output component maps
├── notebooks/             # Jupyter notebooks for analysis
├── requirements.txt       # Python dependencies
└── setup.py               # Package configuration
```

## Usage

### Input Data

Place your ACT frequency maps in the `data/input/` directory. Expected format:
- CAR-projected FITS files
- Two frequencies: 90 GHz and 150 GHz
- Temperature and polarization (T, Q, U)

### Running the Pipeline

```python
from nilc.preprocessing import preprocess_maps
from nilc.needlets import needlet_decompose, needlet_synthesize
from nilc.ilc import compute_ilc_weights, apply_ilc

# Load input maps
# ... (implementation coming soon)

# Preprocess
preprocessed = preprocess_maps(input_maps)

# Needlet decomposition
needlet_bands = needlet_decompose(preprocessed)

# Compute ILC weights and separate components
weights = compute_ilc_weights(needlet_bands)
cmb, ksz = apply_ilc(needlet_bands, weights)

# Synthesize final maps
cmb_map = needlet_synthesize(cmb)
ksz_map = needlet_synthesize(ksz)
```

## Testing

Run the test suite:
```bash
pytest tests/
```

Run with coverage:
```bash
pytest --cov=nilc tests/
```

## Dependencies

- **numpy**: Numerical computations
- **scipy**: Scientific computing and signal processing
- **matplotlib**: Visualization
- **pixell**: CAR map manipulation and spherical harmonic transforms
- **healpy**: HEALPix map format support
- **pytest**: Testing framework

## Scientific Background

The NILC method combines:
- **Needlets**: Localized basis functions in both pixel and harmonic space
- **ILC**: Optimal linear combination to minimize variance while preserving signal
- **Bias mitigation**: Prevents signal-noise correlation artifacts

For details, see Section III of [arXiv:2307.01258](https://arxiv.org/abs/2307.01258).

## Development Status

This is an active development project. Core components are being implemented following the ACT DR6 methodology.

## License

(License information to be added)

## Citation

If you use this code, please cite:
```
Madhavacheril et al. 2023, "The Atacama Cosmology Telescope:
High-resolution component-separated maps across one-third of the sky"
arXiv:2307.01258
```

## Contact

(Contact information to be added)
