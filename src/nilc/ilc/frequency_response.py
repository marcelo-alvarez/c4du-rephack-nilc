"""Frequency response functions for astrophysical components.

These functions describe how each component's brightness varies with frequency,
needed for constrained ILC to preserve or null specific components.
"""

import numpy as np

# Physical constants (SI units)
h_PLANCK = 6.62607015e-34  # Planck constant [J·s]
k_BOLTZMANN = 1.380649e-23  # Boltzmann constant [J/K]
c_LIGHT = 299792458.0  # Speed of light [m/s]
T_CMB = 2.7255  # CMB temperature [K]


def freq_response_cmb(freq_ghz):
    """CMB frequency response in thermodynamic temperature units.

    Parameters
    ----------
    freq_ghz : float or array-like
        Frequency in GHz

    Returns
    -------
    float or array-like
        CMB response (always 1.0 in thermodynamic temperature units)
    """
    freq_ghz = np.asarray(freq_ghz)
    return np.ones_like(freq_ghz, dtype=float)


def freq_response_ksz(freq_ghz):
    """Kinetic SZ frequency response in thermodynamic temperature units.

    The kSZ effect has the same spectral dependence as the CMB.

    Parameters
    ----------
    freq_ghz : float or array-like
        Frequency in GHz

    Returns
    -------
    float or array-like
        kSZ response (always 1.0 in thermodynamic temperature units)
    """
    freq_ghz = np.asarray(freq_ghz)
    return np.ones_like(freq_ghz, dtype=float)


def freq_response_tsz(freq_ghz):
    """Thermal SZ frequency response.

    Implements the spectral function:
    f(x) = x(e^x + 1)/(e^x - 1) - 4
    where x = hν/(k_B T_CMB)

    Parameters
    ----------
    freq_ghz : float or array-like
        Frequency in GHz

    Returns
    -------
    float or array-like
        tSZ response
    """
    freq_ghz = np.asarray(freq_ghz)
    freq_hz = freq_ghz * 1e9  # Convert GHz to Hz

    # Compute dimensionless frequency x = hν/(k_B T_CMB)
    x = h_PLANCK * freq_hz / (k_BOLTZMANN * T_CMB)

    # Compute tSZ spectral function
    # f(x) = x(e^x + 1)/(e^x - 1) - 4
    exp_x = np.exp(x)
    f_tsz = x * (exp_x + 1) / (exp_x - 1) - 4

    return f_tsz


def freq_response_cib(freq_ghz, T_dust=20.0, beta=1.5, freq_ref_ghz=353.0):
    """CIB (Cosmic Infrared Background) frequency response.

    Implements a modified blackbody spectrum:
    f(ν) ∝ ν^beta * B_ν(T_dust)
    normalized at a reference frequency (353 GHz by default).

    Parameters
    ----------
    freq_ghz : float or array-like
        Frequency in GHz
    T_dust : float, optional
        Dust temperature in K (default: 20.0 K)
    beta : float, optional
        Spectral index (default: 1.5)
    freq_ref_ghz : float, optional
        Reference frequency for normalization in GHz (default: 353.0)

    Returns
    -------
    float or array-like
        CIB response normalized to 1.0 at freq_ref_ghz
    """
    freq_ghz = np.asarray(freq_ghz)
    freq_hz = freq_ghz * 1e9  # Convert GHz to Hz
    freq_ref_hz = freq_ref_ghz * 1e9

    # Compute dimensionless frequency x = hν/(k_B T_dust)
    x = h_PLANCK * freq_hz / (k_BOLTZMANN * T_dust)
    x_ref = h_PLANCK * freq_ref_hz / (k_BOLTZMANN * T_dust)

    # Modified blackbody: ν^beta * B_ν(T_dust)
    # B_ν(T) ∝ ν^3 / (e^x - 1) where x = hν/(k_B T)
    # So the full expression is: ν^(beta+3) / (e^x - 1)

    # Compute at target frequency
    numerator = freq_hz**(beta + 3) / (np.exp(x) - 1)

    # Compute at reference frequency for normalization
    denominator = freq_ref_hz**(beta + 3) / (np.exp(x_ref) - 1)

    f_cib = numerator / denominator

    return f_cib


def compute_response_matrix(freq_ghz_list, components=['cmb', 'ksz', 'tsz', 'cib'],
                           cib_params=None):
    """Compute the frequency response matrix F.

    Builds a matrix where F[i, j] is the response of component j at frequency i.

    Parameters
    ----------
    freq_ghz_list : array-like
        List of frequencies in GHz
    components : list of str, optional
        List of components to include (default: ['cmb', 'ksz', 'tsz', 'cib'])
    cib_params : dict, optional
        Parameters for CIB model (T_dust, beta, freq_ref_ghz)

    Returns
    -------
    F : ndarray, shape (n_freq, n_comp)
        Frequency response matrix
    component_names : list of str
        Names of components in the same order as columns of F
    """
    freq_ghz_array = np.asarray(freq_ghz_list)
    n_freq = len(freq_ghz_array)
    n_comp = len(components)

    F = np.zeros((n_freq, n_comp))

    # Set default CIB parameters if not provided
    if cib_params is None:
        cib_params = {'T_dust': 20.0, 'beta': 1.5, 'freq_ref_ghz': 353.0}

    for j, comp in enumerate(components):
        if comp.lower() == 'cmb':
            F[:, j] = freq_response_cmb(freq_ghz_array)
        elif comp.lower() == 'ksz':
            F[:, j] = freq_response_ksz(freq_ghz_array)
        elif comp.lower() == 'tsz':
            F[:, j] = freq_response_tsz(freq_ghz_array)
        elif comp.lower() == 'cib':
            F[:, j] = freq_response_cib(freq_ghz_array, **cib_params)
        else:
            raise ValueError(f"Unknown component: {comp}")

    return F, components


if __name__ == '__main__':
    # Test with ACT frequencies
    freqs_act = np.array([90.0, 150.0])  # GHz

    print("=" * 60)
    print("Frequency Response Functions for ACT Frequencies")
    print("=" * 60)
    print()

    print(f"Frequencies: {freqs_act} GHz")
    print()

    # Test individual response functions
    print("CMB response:")
    cmb_resp = freq_response_cmb(freqs_act)
    for freq, resp in zip(freqs_act, cmb_resp):
        print(f"  {freq:6.1f} GHz: {resp:8.5f}")
    print()

    print("kSZ response:")
    ksz_resp = freq_response_ksz(freqs_act)
    for freq, resp in zip(freqs_act, ksz_resp):
        print(f"  {freq:6.1f} GHz: {resp:8.5f}")
    print()

    print("tSZ response:")
    tsz_resp = freq_response_tsz(freqs_act)
    for freq, resp in zip(freqs_act, tsz_resp):
        print(f"  {freq:6.1f} GHz: {resp:8.5f}")
    print()

    print("CIB response (T_dust=20K, beta=1.5, normalized at 353 GHz):")
    cib_resp = freq_response_cib(freqs_act)
    for freq, resp in zip(freqs_act, cib_resp):
        print(f"  {freq:6.1f} GHz: {resp:8.5f}")
    print()

    # Test response matrix builder
    print("-" * 60)
    print("Full Response Matrix F (frequencies × components)")
    print("-" * 60)
    F, comp_names = compute_response_matrix(freqs_act)

    print(f"\nComponents: {comp_names}")
    print(f"\nMatrix shape: {F.shape} (n_freq={F.shape[0]}, n_comp={F.shape[1]})")
    print("\nF =")

    # Print header
    header = "Freq(GHz) |"
    for comp in comp_names:
        header += f" {comp:>10s} |"
    print(header)
    print("-" * len(header))

    # Print matrix rows
    for i, freq in enumerate(freqs_act):
        row = f"  {freq:6.1f}  |"
        for j in range(len(comp_names)):
            row += f" {F[i, j]:10.5f} |"
        print(row)

    print()
    print("=" * 60)
    print("Physical Constants Used:")
    print("=" * 60)
    print(f"h (Planck constant):     {h_PLANCK:.6e} J·s")
    print(f"k_B (Boltzmann):         {k_BOLTZMANN:.6e} J/K")
    print(f"c (speed of light):      {c_LIGHT:.6e} m/s")
    print(f"T_CMB:                   {T_CMB} K")
    print()
