"""Internal Linear Combination (ILC) component separation with bias mitigation."""

from .covariance import (
    compute_covariance_map,
    compute_covariance_smoothed,
    compute_covariance_donut,
    invert_covariance
)

from .weights import (
    compute_ilc_weights,
    compute_constrained_ilc_weights,
    get_cmb_response,
    get_ksz_response,
    get_tsz_response
)

from .bias import (
    apply_harmonic_exclusion,
    compute_bias_correction,
    select_bias_method,
    get_exclusion_parameters
)

from .separation import (
    apply_ilc_weights,
    separate_components_single_scale,
    separate_components_multiScale,
    synthesize_components,
    nilc_pipeline
)

__all__ = [
    # Covariance
    'compute_covariance_map',
    'compute_covariance_smoothed',
    'compute_covariance_donut',
    'invert_covariance',
    # Weights
    'compute_ilc_weights',
    'compute_constrained_ilc_weights',
    'get_cmb_response',
    'get_ksz_response',
    'get_tsz_response',
    # Bias mitigation
    'apply_harmonic_exclusion',
    'compute_bias_correction',
    'select_bias_method',
    'get_exclusion_parameters',
    # Separation
    'apply_ilc_weights',
    'separate_components_single_scale',
    'separate_components_multiScale',
    'synthesize_components',
    'nilc_pipeline',
]
