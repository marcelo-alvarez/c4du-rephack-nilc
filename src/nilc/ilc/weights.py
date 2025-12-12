"""ILC weight computation for component separation."""
import numpy as np


def compute_ilc_weights(covariance, component_sed):
    """
    Compute ILC weights from covariance and component SED.

    Implements the standard ILC formula:
        w = (C^{-1} · s) / (s^T · C^{-1} · s)

    where C is the covariance matrix and s is the spectral energy distribution
    (SED) vector for the target component.

    Parameters
    ----------
    covariance : ndarray
        Covariance matrix of shape (n_freq, n_freq)
    component_sed : array_like
        SED vector for the component of interest, shape (n_freq,)
        For CMB: [1, 1, ...] (achromatic in thermodynamic units)
        For kSZ: [1, 1, ...] (also achromatic)
        For tSZ: frequency-dependent SED

    Returns
    -------
    weights : ndarray
        ILC weights of shape (n_freq,)
        These weights sum to preserve the component with unit response

    Notes
    -----
    Following Section III-C of arXiv:2307.01258, the ILC weights minimize
    variance while preserving the target component with unit gain.

    The formula ensures:
    1. Minimum variance: weights minimize total output variance
    2. Unit gain constraint: sum_i(w_i * s_i) = 1 for target component
    3. Zero response to uncorrelated noise

    For CMB and kSZ (achromatic), the SED is [1, 1] for two frequencies.
    """
    component_sed = np.asarray(component_sed)

    # Check dimensions
    n_freq = covariance.shape[0]
    if covariance.shape != (n_freq, n_freq):
        raise ValueError(f"Covariance must be square, got shape {covariance.shape}")
    if component_sed.shape[0] != n_freq:
        raise ValueError(
            f"SED length {component_sed.shape[0]} must match "
            f"covariance dimension {n_freq}"
        )

    # Add regularization to avoid numerical issues with near-singular matrices
    # This is standard practice in ILC to handle noise-dominated modes
    epsilon = 1e-10 * np.trace(covariance) / n_freq
    cov_reg = covariance + epsilon * np.eye(n_freq)

    try:
        # Compute C^{-1}
        cov_inv = np.linalg.inv(cov_reg)
    except np.linalg.LinAlgError:
        # If still singular, use pseudo-inverse
        cov_inv = np.linalg.pinv(cov_reg)

    # Compute numerator: C^{-1} · s
    numerator = cov_inv @ component_sed

    # Compute denominator: s^T · C^{-1} · s
    denominator = component_sed @ numerator

    if np.abs(denominator) < 1e-15:
        raise ValueError(
            "Denominator in ILC formula is zero. Check covariance and SED."
        )

    # Final weights
    weights = numerator / denominator

    return weights


def validate_weights(weights, component_sed, tolerance=1e-6):
    """
    Validate that ILC weights satisfy the unit gain constraint.

    Parameters
    ----------
    weights : ndarray
        ILC weights
    component_sed : array_like
        Component SED vector
    tolerance : float
        Acceptable deviation from unity

    Returns
    -------
    is_valid : bool
        True if weights satisfy constraint
    gain : float
        Actual gain (should be ~1.0)

    Notes
    -----
    The constraint sum_i(w_i * s_i) = 1 ensures the component is preserved
    with unit response while other components are suppressed.
    """
    component_sed = np.asarray(component_sed)
    gain = np.sum(weights * component_sed)
    is_valid = np.abs(gain - 1.0) < tolerance

    return is_valid, gain


def compute_ilc_weights_multicomponent(covariance, component_seds):
    """
    Compute ILC weights for multiple components simultaneously.

    For separating N components from N frequencies, we need to solve for
    N sets of weights, one per component.

    Parameters
    ----------
    covariance : ndarray
        Covariance matrix of shape (n_freq, n_freq)
    component_seds : list of array_like
        List of SED vectors, one per component to separate
        Each SED has shape (n_freq,)

    Returns
    -------
    weights_dict : dict
        Dictionary mapping component index to weights array
        {0: weights_for_component_0, 1: weights_for_component_1, ...}

    Notes
    -----
    For the ACT NILC pipeline with 2 frequencies and 2 components (CMB, kSZ),
    both SEDs are [1, 1], but the weights will differ due to different
    noise properties and correlations captured in the covariance.

    In practice, the covariance varies by needlet scale, leading to
    scale-dependent weights even when SEDs are constant.
    """
    weights_dict = {}

    for i, sed in enumerate(component_seds):
        weights_dict[i] = compute_ilc_weights(covariance, sed)

    return weights_dict
