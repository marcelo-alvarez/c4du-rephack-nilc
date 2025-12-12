"""Plot frequency response functions for astrophysical components.

Creates a visualization similar to the reference plot showing how different
components vary with frequency, normalized at 150 GHz.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Add src to path to import our module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.nilc.ilc.frequency_response import (
    freq_response_cmb,
    freq_response_ksz,
    freq_response_tsz,
    freq_response_cib,
)


def plot_frequency_responses(output_file='frequency_responses_150GHz.png'):
    """Create and save frequency response plot normalized at 150 GHz.

    Parameters
    ----------
    output_file : str
        Output filename for the plot
    """
    # Frequency range from ~50 GHz to 700 GHz
    freqs = np.linspace(50, 700, 500)

    # Normalization frequency
    freq_norm = 150.0  # GHz

    # Compute responses
    cmb = freq_response_cmb(freqs)
    ksz = freq_response_ksz(freqs)
    tsz = freq_response_tsz(freqs)
    cib = freq_response_cib(freqs, T_dust=20.0, beta=1.5, freq_ref_ghz=353.0)

    # Normalize all responses by their value at 150 GHz
    cmb_norm = freq_response_cmb(freq_norm)
    ksz_norm = freq_response_ksz(freq_norm)
    tsz_norm = freq_response_tsz(freq_norm)
    cib_norm = freq_response_cib(freq_norm, T_dust=20.0, beta=1.5, freq_ref_ghz=353.0)

    # Normalize (use absolute value for tSZ to match reference plot style)
    cmb_plot = cmb / cmb_norm
    ksz_plot = ksz / ksz_norm
    tsz_plot = tsz / np.abs(tsz_norm)
    cib_plot = cib / cib_norm

    # ACT frequency bands to highlight
    act_freqs = [90, 150]

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 7))

    # Plot frequency responses
    ax.plot(freqs, cmb_plot, '-', color='purple', linewidth=2.5,
            label='CMB', alpha=0.8)
    ax.plot(freqs, tsz_plot, '--', color='darkorange', linewidth=2.5,
            label='Nonrelativistic tSZ', alpha=0.8)
    ax.plot(freqs, cib_plot, '-.', color='darkgreen', linewidth=2.5,
            label='CIB', alpha=0.8)

    # Add vertical bands for ACT frequencies
    for freq in act_freqs:
        ax.axvspan(freq - 10, freq + 10, alpha=0.15, color='green', zorder=0)

    # Add horizontal line at y=1 (normalization reference)
    ax.axhline(1, color='purple', linestyle='-', linewidth=1.5, alpha=0.5)

    # Formatting
    ax.set_xlabel(r'$\nu$ (GHz)', fontsize=16)
    ax.set_ylabel(r'$f_x(\nu) / |f_x(150)|$', fontsize=16)
    ax.set_xlim(50, 700)
    ax.set_ylim(-2.5, 4.5)

    # Grid
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)

    # Legend
    ax.legend(fontsize=13, loc='upper left', framealpha=0.9)

    # Title
    ax.set_title('Frequency Response Functions (normalized at 150 GHz)',
                 fontsize=14, pad=15)

    # Tick label sizes
    ax.tick_params(labelsize=12)

    # Tight layout
    plt.tight_layout()

    # Save figure
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {output_file}")

    # Also display some numerical values
    print(f"\nNormalization values at {freq_norm} GHz:")
    print(f"  CMB:  {cmb_norm:8.5f}")
    print(f"  kSZ:  {ksz_norm:8.5f}")
    print(f"  tSZ:  {tsz_norm:8.5f} (abs: {np.abs(tsz_norm):8.5f})")
    print(f"  CIB:  {cib_norm:8.5f}")

    print(f"\nValues at ACT frequencies:")
    for freq in act_freqs:
        print(f"\n{freq} GHz:")
        print(f"  CMB:  {freq_response_cmb(freq):8.5f}")
        print(f"  tSZ:  {freq_response_tsz(freq):8.5f}")
        print(f"  CIB:  {freq_response_cib(freq):8.5f}")

    return fig, ax


if __name__ == '__main__':
    # Allow custom output filename from command line
    output_file = sys.argv[1] if len(sys.argv) > 1 else 'frequency_responses_150GHz.png'
    plot_frequency_responses(output_file)
    print("\nDone!")
