"""Internal Linear Combination (ILC) component separation with bias mitigation."""

from .frequency_response import (
    freq_response_cmb,
    freq_response_ksz,
    freq_response_tsz,
    freq_response_cib,
    compute_response_matrix,
    h_PLANCK,
    k_BOLTZMANN,
    c_LIGHT,
    T_CMB,
)

from .covariance import (
    compute_covariance_matrix,
    apply_harmonic_exclusion,
    apply_donut_smoothing,
    compute_ilc_covariance,
)

__all__ = [
    # Frequency response functions
    'freq_response_cmb',
    'freq_response_ksz',
    'freq_response_tsz',
    'freq_response_cib',
    'compute_response_matrix',
    'h_PLANCK',
    'k_BOLTZMANN',
    'c_LIGHT',
    'T_CMB',
    # Covariance computation and bias mitigation
    'compute_covariance_matrix',
    'apply_harmonic_exclusion',
    'apply_donut_smoothing',
    'compute_ilc_covariance',
]
