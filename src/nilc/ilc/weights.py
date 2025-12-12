"""
ILC weight computation for component separation.

Implements the Internal Linear Combination weight calculation per
Section III-C of arXiv:2307.01258.
"""
import numpy as np


def compute_ilc_weights(cov_inv, component_vector):
    """
    Compute ILC weights from inverse covariance matrix.

    The ILC weights minimize variance while preserving the desired component:
        w = C^(-1) @ e / (e^T @ C^(-1) @ e)

    where:
        C is the covariance matrix
        e is the component response vector (spectral signature)
        w is the weight vector

    Parameters
    ----------
    cov_inv : ndarray
        Inverse covariance matrices with shape (n_y, n_x, n_freq, n_freq)
    component_vector : ndarray
        Component response vector with shape (n_freq,)
        For CMB: [1, 1] (same signal in both frequencies)
        For kSZ: [1, 1] (also achromatic)
        For tSZ: spectral dependence on frequency

    Returns
    -------
    weights : ndarray
        ILC weights with shape (n_y, n_x, n_freq)

    Notes
    -----
    For CMB and kSZ (both achromatic), we need additional constraints or
    priors to separate them, since they have the same spectral signature.
    This can be done with:
    - Multi-component ILC (constrained ILC)
    - External templates or priors
    - Different needlet scales (CMB dominates large scales, kSZ small scales)
    """
    n_y, n_x, n_freq, _ = cov_inv.shape
    weights = np.zeros((n_y, n_x, n_freq))

    # Ensure component vector is column vector for matrix operations
    e = component_vector.reshape(-1, 1)  # (n_freq, 1)

    # Compute weights at each pixel
    for iy in range(n_y):
        for ix in range(n_x):
            C_inv = cov_inv[iy, ix]  # (n_freq, n_freq)

            # Numerator: C^(-1) @ e
            numerator = C_inv @ e  # (n_freq, 1)

            # Denominator: e^T @ C^(-1) @ e
            denominator = (e.T @ C_inv @ e)[0, 0]  # scalar

            # Weights: w = numerator / denominator
            if abs(denominator) > 1e-15:
                weights[iy, ix, :] = (numerator / denominator).flatten()
            else:
                # If denominator is too small, set equal weights
                weights[iy, ix, :] = 1.0 / n_freq

    return weights


def compute_constrained_ilc_weights(cov_inv, response_matrix):
    """
    Compute constrained ILC weights for multi-component separation.

    For separating multiple components (e.g., CMB + kSZ), we need constrained ILC
    that preserves all desired components while minimizing variance.

    The constrained ILC weights are:
        W = C^(-1) @ A @ (A^T @ C^(-1) @ A)^(-1)

    where:
        C is covariance matrix
        A is response matrix (each column is a component's spectral signature)
        W is weight matrix (each column gives weights for one component)

    Parameters
    ----------
    cov_inv : ndarray
        Inverse covariance matrices with shape (n_y, n_x, n_freq, n_freq)
    response_matrix : ndarray
        Component response matrix with shape (n_freq, n_comp)
        Each column represents one component's spectral signature
        For CMB+kSZ: [[1, 1], [1, 1]]  # Both achromatic

    Returns
    -------
    weights : ndarray
        Weight matrix with shape (n_y, n_x, n_freq, n_comp)
        weights[:,:,:,i] gives the ILC weights for component i

    Notes
    -----
    For CMB and kSZ which are both achromatic, additional constraints or
    regularization is needed. Options include:
    - Use different scales (CMB power at large ell, kSZ at small ell)
    - Add deprojection of other components (e.g., deproject tSZ)
    - Use external templates or priors
    """
    n_y, n_x, n_freq, _ = cov_inv.shape
    n_comp = response_matrix.shape[1]
    weights = np.zeros((n_y, n_x, n_freq, n_comp))

    A = response_matrix  # (n_freq, n_comp)

    # Compute weights at each pixel
    for iy in range(n_y):
        for ix in range(n_x):
            C_inv = cov_inv[iy, ix]  # (n_freq, n_freq)

            # Compute A^T @ C^(-1) @ A
            AtCinvA = A.T @ C_inv @ A  # (n_comp, n_comp)

            try:
                # Invert (A^T @ C^(-1) @ A)
                AtCinvA_inv = np.linalg.inv(AtCinvA)

                # Compute weights: W = C^(-1) @ A @ (A^T @ C^(-1) @ A)^(-1)
                weights[iy, ix, :, :] = C_inv @ A @ AtCinvA_inv

            except np.linalg.LinAlgError:
                # If inversion fails, use pseudoinverse
                AtCinvA_inv = np.linalg.pinv(AtCinvA)
                weights[iy, ix, :, :] = C_inv @ A @ AtCinvA_inv

    return weights


def get_cmb_response(frequencies_ghz):
    """
    Get CMB spectral response (achromatic in thermodynamic temperature).

    Parameters
    ----------
    frequencies_ghz : array_like
        Observing frequencies in GHz

    Returns
    -------
    response : ndarray
        CMB response vector (all ones for thermodynamic temperature)
    """
    return np.ones(len(frequencies_ghz))


def get_ksz_response(frequencies_ghz):
    """
    Get kSZ spectral response (achromatic, same as CMB).

    Parameters
    ----------
    frequencies_ghz : array_like
        Observing frequencies in GHz

    Returns
    -------
    response : ndarray
        kSZ response vector (all ones, achromatic)
    """
    return np.ones(len(frequencies_ghz))


def get_tsz_response(frequencies_ghz, T_cmb=2.7255):
    """
    Get tSZ spectral response.

    The tSZ has a characteristic frequency-dependent signature:
        f(nu) = x * coth(x/2) - 4
    where x = h*nu / (k*T_CMB)

    Parameters
    ----------
    frequencies_ghz : array_like
        Observing frequencies in GHz
    T_cmb : float
        CMB temperature in Kelvin

    Returns
    -------
    response : ndarray
        tSZ spectral response

    Notes
    -----
    This is the non-relativistic tSZ spectrum. For high-accuracy work,
    relativistic corrections should be included.
    """
    # Physical constants
    h = 6.62607015e-34  # Planck constant (J·s)
    k = 1.380649e-23    # Boltzmann constant (J/K)

    # Convert frequencies to Hz
    nu_hz = np.array(frequencies_ghz) * 1e9

    # Compute x = h*nu / (k*T_CMB)
    x = h * nu_hz / (k * T_cmb)

    # tSZ spectral function: f(x) = x * coth(x/2) - 4
    response = x / np.tanh(x / 2.0) - 4.0

    return response
