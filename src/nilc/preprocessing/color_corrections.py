"""
Color correction utilities for ACT maps.

Implements color corrections to account for component spectral responses
following Section III-A of arXiv:2307.01258.

Different astrophysical components (CMB, kSZ, tSZ, CIB) have different spectral
energy distributions (SEDs), which means their effective observing frequencies
differ from the nominal bandpass center. Color corrections normalize maps to
a common reference (e.g., thermodynamic temperature for CMB).
"""
import numpy as np


# Physical constants
h_PLANCK = 6.62607015e-34  # J·s
k_BOLTZMANN = 1.380649e-23  # J/K
T_CMB = 2.7255  # K, CMB temperature


def load_passband(passband_file):
    """
    Load ACT passband transmission data from file.

    Parameters
    ----------
    passband_file : str
        Path to passband file with columns: frequency [GHz], response, error

    Returns
    -------
    freq : ndarray
        Frequency array in GHz
    response : ndarray
        Normalized passband response (dimensionless)

    Notes
    -----
    ACT passband files contain three columns:
    - Column 0: Frequency in GHz
    - Column 1: Transmission response
    - Column 2: Uncertainty/error on response

    The response is normalized such that the integral over frequency equals 1.
    """
    data = np.loadtxt(passband_file)
    freq = data[:, 0]
    response = data[:, 1]

    # Filter out zero or negative responses
    valid_mask = response > 0
    freq = freq[valid_mask]
    response = response[valid_mask]

    # Normalize the passband
    integral = np.trapz(response, freq)
    if integral > 0:
        response = response / integral
    else:
        raise ValueError(f"Passband has zero or negative integral: {passband_file}")

    return freq, response


def compute_color_correction(passband_freq, passband_response, component_name):
    """
    Compute color correction factor for a given component.

    The color correction accounts for the fact that different components have
    different spectral energy distributions (SEDs), so the effective frequency
    and calibration differ from the nominal bandpass center.

    Parameters
    ----------
    passband_freq : ndarray
        Frequency array in GHz
    passband_response : ndarray
        Normalized passband transmission
    component_name : str
        Component type: 'CMB', 'kSZ', 'tSZ', 'CIB'

    Returns
    -------
    correction_factor : float
        Multiplicative color correction factor to convert measured brightness
        to thermodynamic temperature units
    effective_freq : float
        Effective frequency in GHz for this component

    Notes
    -----
    The effective frequency is defined as:

    .. math::
        \\nu_{eff} = \\frac{\\int \\nu \\, T(\\nu) \\, S_{comp}(\\nu) \\, d\\nu}
                          {\\int T(\\nu) \\, S_{comp}(\\nu) \\, d\\nu}

    where T(ν) is the passband transmission and S_comp(ν) is the component SED.

    For CMB and kSZ, we use the thermodynamic temperature derivative:

    .. math::
        \\frac{dB_\\nu}{dT} = \\frac{2h\\nu^3}{c^2} \\frac{x^2 e^x}{(e^x - 1)^2}

    where x = hν/(k_B T_CMB) and B_ν is the Planck function.

    References
    ----------
    arXiv:2307.01258, Section III-A (ACT DR6 component separation)
    """
    # Convert frequency to Hz for physical calculations
    freq_hz = passband_freq * 1e9

    # Compute x = h*nu / (k_B * T_CMB)
    x = h_PLANCK * freq_hz / (k_BOLTZMANN * T_CMB)

    if component_name.upper() == 'CMB':
        # CMB: dB_nu/dT in thermodynamic temperature
        # dB/dT = (2h nu^3/c^2) * x^2 * exp(x) / (exp(x) - 1)^2
        sed = x**2 * np.exp(x) / (np.exp(x) - 1)**2

    elif component_name.upper() == 'KSZ':
        # kSZ has the same spectral signature as CMB
        # (it's a Doppler shift of CMB photons)
        sed = x**2 * np.exp(x) / (np.exp(x) - 1)**2

    elif component_name.upper() == 'TSZ':
        # tSZ: thermal Sunyaev-Zeldovich effect
        # Placeholder - full implementation requires tSZ spectral function
        # f(x) = x * (exp(x)+1)/(exp(x)-1) - 4
        raise NotImplementedError(
            "tSZ color correction not yet implemented. "
            "Will be added in future iteration."
        )

    elif component_name.upper() == 'CIB':
        # CIB: cosmic infrared background
        # Placeholder - requires modified blackbody spectrum
        raise NotImplementedError(
            "CIB color correction not yet implemented. "
            "Will be added in future iteration."
        )

    else:
        raise ValueError(
            f"Unknown component: {component_name}. "
            "Supported: 'CMB', 'kSZ', 'tSZ', 'CIB'"
        )

    # Compute effective frequency
    numerator = np.trapz(passband_freq * passband_response * sed, passband_freq)
    denominator = np.trapz(passband_response * sed, passband_freq)
    effective_freq = numerator / denominator

    # Compute color correction as the ratio of SED at effective frequency
    # to the passband-averaged SED
    # For simplicity, we compute the SED at the effective frequency
    x_eff = h_PLANCK * (effective_freq * 1e9) / (k_BOLTZMANN * T_CMB)
    sed_eff = x_eff**2 * np.exp(x_eff) / (np.exp(x_eff) - 1)**2

    # The correction factor normalizes the passband-averaged response
    # to the value at the effective frequency
    sed_avg = np.trapz(passband_response * sed, passband_freq) / np.trapz(passband_response, passband_freq)
    correction_factor = sed_eff / sed_avg

    return correction_factor, effective_freq


def apply_color_correction(map_90ghz, map_150ghz, component='CMB',
                          passband_dir='/home/jiaqu/NILC/data/'):
    """
    Apply color corrections to multi-frequency maps.

    Loads passbands and applies color corrections to normalize maps to the
    same thermodynamic temperature units for the specified component.

    Parameters
    ----------
    map_90ghz : enmap.ndmap or ndarray
        90 GHz map
    map_150ghz : enmap.ndmap or ndarray
        150 GHz map
    component : str, optional
        Component to correct for: 'CMB', 'kSZ', 'tSZ', 'CIB'
        Default is 'CMB'
    passband_dir : str, optional
        Directory containing passband files
        Default is '/home/jiaqu/NILC/data/'

    Returns
    -------
    map_90ghz_corrected : same type as input
        Color-corrected 90 GHz map
    map_150ghz_corrected : same type as input
        Color-corrected 150 GHz map
    correction_90 : float
        Color correction factor applied to 90 GHz map
    correction_150 : float
        Color correction factor applied to 150 GHz map
    eff_freq_90 : float
        Effective frequency (GHz) for 90 GHz map
    eff_freq_150 : float
        Effective frequency (GHz) for 150 GHz map

    Notes
    -----
    The color correction ensures that both maps are calibrated in the same
    units for the specified component. For CMB and kSZ, this is thermodynamic
    temperature.

    Examples
    --------
    >>> map_90_corr, map_150_corr, c90, c150, f90, f150 = apply_color_correction(
    ...     map_90, map_150, component='CMB'
    ... )
    >>> print(f"90 GHz effective frequency: {f90:.2f} GHz")
    >>> print(f"150 GHz effective frequency: {f150:.2f} GHz")
    """
    import os

    # Load PA5 passbands (truncated versions with normalized responses)
    passband_90_file = os.path.join(passband_dir, 'PA5_avg_passband_090_wErr_trunc.txt')
    passband_150_file = os.path.join(passband_dir, 'PA5_avg_passband_150_wErr_trunc.txt')

    freq_90, response_90 = load_passband(passband_90_file)
    freq_150, response_150 = load_passband(passband_150_file)

    # Compute color corrections
    correction_90, eff_freq_90 = compute_color_correction(
        freq_90, response_90, component
    )
    correction_150, eff_freq_150 = compute_color_correction(
        freq_150, response_150, component
    )

    # Apply corrections
    map_90ghz_corrected = map_90ghz * correction_90
    map_150ghz_corrected = map_150ghz * correction_150

    return (map_90ghz_corrected, map_150ghz_corrected,
            correction_90, correction_150,
            eff_freq_90, eff_freq_150)
