# -*- coding: utf-8 -*-
"""Plot frequency response functions for astrophysical components.

Creates a visualization similar to the reference plot showing how different
components vary with frequency, normalized at 150 GHz.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from scipy.interpolate import interp1d

# Add src to path to import our module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.nilc.ilc.frequency_response import (
    freq_response_cmb,
    freq_response_ksz,
    freq_response_tsz,
    freq_response_cib,
)


def estimate_target_curves(freqs):
    """Estimate curves from the target plot (Fig 6) based on key points.
    
    Returns estimated values for:
    - Non-relativistic tSZ (gold/orange)
    - Relativistic tSZ (T_e = 8 keV) (blue)
    - CIB (green, dash-dot)
    - CIB-δβ (orange-brown, dash-dot)
    - CMB (pink/magenta, constant at 1)
    """
    # Key points from the image description
    # Non-relativistic tSZ (gold/orange, solid)
    tsz_nr_freqs = np.array([0, 100, 150, 217, 280, 330, 370, 400, 500, 600, 700])
    tsz_nr_vals = np.array([-2.1, -1.5, -1.0, 0.0, 1.0, 2.0, 3.0, 3.9, 4.5, 5.0, 5.5])
    tsz_nr_interp = interp1d(tsz_nr_freqs, tsz_nr_vals, kind='cubic', 
                             bounds_error=False, fill_value='extrapolate')
    
    # Relativistic tSZ (T_e = 8 keV) (blue, solid)
    tsz_rel_freqs = np.array([0, 100, 150, 200, 260, 310, 350, 390, 500, 600, 700])
    tsz_rel_vals = np.array([-2.1, -1.4, -0.8, 0.0, 1.0, 2.0, 3.0, 3.9, 4.6, 5.1, 5.6])
    tsz_rel_interp = interp1d(tsz_rel_freqs, tsz_rel_vals, kind='cubic',
                              bounds_error=False, fill_value='extrapolate')
    
    # CIB (green, dash-dot) - exponential growth
    # From description: 0 at 0 GHz, 1 at 150 GHz, 2 at 200 GHz, 4 at 250 GHz, continues rapidly
    cib_freqs = np.array([0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700])
    cib_vals = np.array([0.0, 0.05, 0.3, 1.0, 2.0, 4.0, 6.5, 9.5, 13.5, 18.5, 25.0, 45.0, 80.0])
    cib_interp = interp1d(cib_freqs, cib_vals, kind='cubic',
                          bounds_error=False, fill_value='extrapolate')
    
    # CIB-δβ (orange-brown, dash-dot) - negative, decreasing
    # From description: 0 at 0 GHz, -1 at 150 GHz, -1.5 at 250 GHz, -2.5 at 350 GHz, -3.9 at 400 GHz
    cib_dbeta_freqs = np.array([0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700])
    cib_dbeta_vals = np.array([0.0, -0.2, -0.6, -1.0, -1.3, -1.5, -2.0, -2.5, -3.9, -4.8, -5.8, -7.5, -9.5])
    cib_dbeta_interp = interp1d(cib_dbeta_freqs, cib_dbeta_vals, kind='cubic',
                                bounds_error=False, fill_value='extrapolate')
    
    # CMB (pink/magenta, constant at 1)
    cmb_vals = np.ones_like(freqs)
    
    return {
        'tsz_nr': tsz_nr_interp(freqs),
        'tsz_rel': tsz_rel_interp(freqs),
        'cib': cib_interp(freqs),
        'cib_dbeta': cib_dbeta_interp(freqs),
        'cmb': cmb_vals
    }


def plot_frequency_responses(output_file='frequency_responses_150GHz.png'):
    """Create and save frequency response plot normalized at 150 GHz.

    Parameters
    ----------
    output_file : str
        Output filename for the plot
    """
    # Frequency range from 0 to 700 GHz to match target plot
    freqs = np.linspace(0, 700, 1000)

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

    # Estimate target curves from Fig 6
    target_curves = estimate_target_curves(freqs)

    # ACT frequency bands to highlight
    act_freqs = [90, 150]

    # Create figure with subplots - main plot on top, residual plot on bottom
    fig = plt.figure(figsize=(5, 8))
    gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.35, left=0.15, right=0.95, top=0.95, bottom=0.1)
    ax = fig.add_subplot(gs[0])
    ax_res = fig.add_subplot(gs[1])

    # Plot target curves first (in background) with fatter linewidth and lighter alpha
    # Order: Non-relativistic tSZ, CIB, CMB (matching paper figure order)
    ax.plot(freqs, target_curves['tsz_nr'], '-', color='gold', linewidth=4.0,
            alpha=0.4, label='Non-relativistic tSZ (target)', zorder=0)
    ax.plot(freqs, target_curves['cib'], '-.', color='lightgreen', linewidth=4.0,
            alpha=0.4, label='CIB (target)', zorder=0)
    ax.plot(freqs, target_curves['cmb'], '-', color='pink', linewidth=4.0,
            alpha=0.4, label='CMB (target)', zorder=0)

    # Plot our computed frequency responses (replication) on top
    # Order: Non-relativistic tSZ, CIB, CMB (matching paper figure order)
    ax.plot(freqs, tsz_plot, '--', color='darkorange', linewidth=2.5,
            label='Non-relativistic tSZ (replication)', alpha=0.8)
    ax.plot(freqs, cib_plot, '-.', color='darkgreen', linewidth=2.5,
            label='CIB (replication)', alpha=0.8)
    ax.plot(freqs, cmb_plot, '-', color='purple', linewidth=2.5,
            label='CMB (replication)', alpha=0.8)

    # Add vertical bands for ACT frequencies
    for freq in act_freqs:
        ax.axvspan(freq - 10, freq + 10, alpha=0.15, color='green', zorder=0)

    # Add horizontal line at y=1 (normalization reference)
    ax.axhline(1, color='purple', linestyle='-', linewidth=1.5, alpha=0.5)

    # Formatting - match target plot axes
    ax.set_xlabel(r'$\nu$ (GHz)', fontsize=16)
    ax.set_ylabel(r'$f_X(\nu) / |f_X(150)|$', fontsize=16)
    ax.set_xlim(0, 700)
    ax.set_ylim(-3, 4)

    # Grid
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)

    # Legend - order matches figure: tSZ, CIB, CMB (each with target then replication)
    # Get handles and labels to reorder
    handles, labels = ax.get_legend_handles_labels()
    # Reorder to: tSZ (target), tSZ (replication), CIB (target), CIB (replication), CMB (target), CMB (replication)
    desired_order = [0, 3, 1, 4, 2, 5]  # Indices after plotting in order above
    handles_ordered = [handles[i] for i in desired_order]
    labels_ordered = [labels[i] for i in desired_order]
    ax.legend(handles_ordered, labels_ordered, fontsize=11, loc='lower right', framealpha=0.9)

    # Tick label sizes
    ax.tick_params(labelsize=12)

    # Calculate fractional residuals (replication - target) / target
    # Avoid division by zero
    residual_tsz = np.where(np.abs(target_curves['tsz_nr']) > 1e-6,
                           (tsz_plot - target_curves['tsz_nr']) / target_curves['tsz_nr'],
                           np.nan)
    residual_cib = np.where(np.abs(target_curves['cib']) > 1e-6,
                           (cib_plot - target_curves['cib']) / target_curves['cib'],
                           np.nan)
    residual_cmb = np.where(np.abs(target_curves['cmb']) > 1e-6,
                           (cmb_plot - target_curves['cmb']) / target_curves['cmb'],
                           np.nan)

    # Plot fractional residuals in bottom subpanel
    ax_res.plot(freqs, residual_tsz, '--', color='darkorange', linewidth=2.0,
                label='tSZ residual', alpha=0.8)
    ax_res.plot(freqs, residual_cib, '-.', color='darkgreen', linewidth=2.0,
                label='CIB residual', alpha=0.8)
    ax_res.plot(freqs, residual_cmb, '-', color='purple', linewidth=2.0,
                label='CMB residual', alpha=0.8)
    ax_res.axhline(0, color='black', linestyle='-', linewidth=1.0, alpha=0.3)

    # Format residual subpanel
    ax_res.set_xlabel(r'$\nu$ (GHz)', fontsize=14)
    ax_res.set_ylabel('Fractional Residual', fontsize=12)
    ax_res.set_xlim(0, 700)
    ax_res.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
    ax_res.tick_params(labelsize=10)
    ax_res.legend(fontsize=9, loc='lower right', framealpha=0.9)

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
    output_file = sys.argv[1] if len(sys.argv) > 1 else 'frequency_responses_150GHz_validation.png'
    plot_frequency_responses(output_file)
    print("\nDone!")
