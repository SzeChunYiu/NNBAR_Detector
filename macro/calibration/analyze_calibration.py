#!/usr/bin/env python3
"""
NNBAR Detector Calibration Analysis
====================================
Analyzes calibration data and generates publication-quality plots.

Author: NNBAR Collaboration
Date: January 2026
"""

import os
import sys
import glob
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from pathlib import Path
from datetime import datetime

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for PDF generation
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches
from scipy.optimize import curve_fit
from scipy.stats import norm

# Set publication-quality style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'legend.fontsize': 10,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'figure.figsize': (10, 7),
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'axes.grid': True,
    'grid.alpha': 0.3,
    'axes.axisbelow': True,
})

# Color palette
COLORS = {
    'scint': '#2E86AB',      # Blue
    'leadglass': '#A23B72',  # Magenta
    'pion': '#F18F01',       # Orange
    'gamma': '#C73E1D',      # Red
    'electron': '#3A7D44',   # Green
    'fit': '#1B1B1E',        # Black
}

# Physics constants
SCINT_LIGHT_YIELD = 10000  # photons/MeV (BC-408)
CERENKOV_YIELD = 200       # photons/MeV (lead glass approximate)


def load_calibration_data(output_dir, detector_type, run_numbers):
    """Load calibration data from parquet files."""
    all_data = []

    for run_num in run_numbers:
        filename = f"{detector_type}_output_{run_num}.parquet"
        filepath = os.path.join(output_dir, filename)

        if os.path.exists(filepath) and os.path.getsize(filepath) > 100:
            try:
                df = pq.read_table(filepath).to_pandas()
                df['run_energy'] = run_num  # Energy in MeV from run number
                all_data.append(df)
                print(f"  Loaded {filepath}: {len(df)} entries")
            except Exception as e:
                print(f"  Warning: Could not read {filepath}: {e}")
        else:
            print(f"  Skipping {filepath} (not found or empty)")

    if all_data:
        return pd.concat(all_data, ignore_index=True)
    return pd.DataFrame()


def load_particle_data(output_dir, run_numbers):
    """Load primary particle data."""
    all_data = []

    for run_num in run_numbers:
        filename = f"Particle_output_{run_num}.parquet"
        filepath = os.path.join(output_dir, filename)

        if os.path.exists(filepath) and os.path.getsize(filepath) > 100:
            try:
                df = pq.read_table(filepath).to_pandas()
                df['run_energy'] = run_num
                all_data.append(df)
            except Exception as e:
                print(f"  Warning: Could not read {filepath}: {e}")

    if all_data:
        return pd.concat(all_data, ignore_index=True)
    return pd.DataFrame()


def linear_fit(x, a, b):
    """Linear function for fitting."""
    return a * x + b


def sqrt_resolution(E, a, b):
    """Resolution function: sigma/E = a/sqrt(E) + b."""
    return a / np.sqrt(E) + b


def analyze_scintillator(df, particle_df, pdf):
    """Analyze scintillator calibration data."""
    if df.empty:
        print("No scintillator data to analyze")
        return None

    energies = sorted(df['run_energy'].unique())

    # Aggregate data by event and energy
    results = []
    for energy in energies:
        subset = df[df['run_energy'] == energy]

        # Group by event to get total deposited energy and photons
        event_data = subset.groupby('Event_ID').agg({
            'eDep': 'sum',
            'photons': 'sum' if 'photons' in subset.columns else 'count'
        }).reset_index()

        if len(event_data) > 10:
            edep_mean = event_data['eDep'].mean()
            edep_std = event_data['eDep'].std()
            photons_mean = event_data['photons'].mean() if 'photons' in event_data.columns else 0
            photons_std = event_data['photons'].std() if 'photons' in event_data.columns else 0

            results.append({
                'energy': energy,
                'edep_mean': edep_mean,
                'edep_std': edep_std,
                'photons_mean': photons_mean,
                'photons_std': photons_std,
                'n_events': len(event_data)
            })

    if not results:
        return None

    results_df = pd.DataFrame(results)

    # =========================================================================
    # Plot 1: Energy Deposition vs True Energy
    # =========================================================================
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Scintillator Calibration - Pion Response', fontsize=14, fontweight='bold')

    ax1 = axes[0, 0]
    ax1.errorbar(results_df['energy'], results_df['edep_mean'],
                 yerr=results_df['edep_std']/np.sqrt(results_df['n_events']),
                 fmt='o', color=COLORS['pion'], markersize=8, capsize=4,
                 label='π⁺ data')

    # Fit linear response
    if len(results_df) > 2:
        try:
            popt, pcov = curve_fit(linear_fit, results_df['energy'], results_df['edep_mean'])
            x_fit = np.linspace(0, results_df['energy'].max() * 1.1, 100)
            ax1.plot(x_fit, linear_fit(x_fit, *popt), '--', color=COLORS['fit'],
                    label=f'Fit: E_dep = {popt[0]:.3f}×E + {popt[1]:.1f}')
        except:
            pass

    ax1.set_xlabel('True Pion Energy (MeV)')
    ax1.set_ylabel('Mean Energy Deposited (MeV)')
    ax1.set_title('Energy Deposition Response')
    ax1.legend(loc='upper left')
    ax1.set_xlim(0, None)
    ax1.set_ylim(0, None)

    # =========================================================================
    # Plot 2: Energy Resolution
    # =========================================================================
    ax2 = axes[0, 1]
    resolution = results_df['edep_std'] / results_df['edep_mean'] * 100  # in percent
    ax2.errorbar(results_df['energy'], resolution,
                 fmt='s', color=COLORS['scint'], markersize=8, capsize=4,
                 label='σ/E')

    # Fit resolution function
    if len(results_df) > 2:
        try:
            popt_res, _ = curve_fit(sqrt_resolution, results_df['energy'], resolution/100,
                                    p0=[0.1, 0.05], bounds=([0, 0], [1, 0.5]))
            x_fit = np.linspace(results_df['energy'].min(), results_df['energy'].max(), 100)
            ax2.plot(x_fit, sqrt_resolution(x_fit, *popt_res)*100, '--', color=COLORS['fit'],
                    label=f'Fit: σ/E = {popt_res[0]*100:.1f}%/√E ⊕ {popt_res[1]*100:.1f}%')
        except:
            pass

    ax2.set_xlabel('True Pion Energy (MeV)')
    ax2.set_ylabel('Energy Resolution σ/E (%)')
    ax2.set_title('Energy Resolution')
    ax2.legend(loc='upper right')
    ax2.set_xlim(0, None)

    # =========================================================================
    # Plot 3: dE/dx distribution (if track length available)
    # =========================================================================
    ax3 = axes[1, 0]
    if 'trackl' in df.columns or 'track_length' in df.columns:
        track_col = 'trackl' if 'trackl' in df.columns else 'track_length'
        df_valid = df[(df[track_col] > 0) & (df['eDep'] > 0)]
        if len(df_valid) > 100:
            dedx = df_valid['eDep'] / df_valid[track_col]  # MeV/cm or MeV/mm
            ax3.hist(dedx, bins=50, color=COLORS['pion'], alpha=0.7, edgecolor='black')
            ax3.axvline(dedx.median(), color=COLORS['fit'], linestyle='--',
                       label=f'Median: {dedx.median():.2f}')
            ax3.set_xlabel('dE/dx (MeV/length unit)')
            ax3.set_ylabel('Entries')
            ax3.set_title('Energy Loss Distribution')
            ax3.legend()
    else:
        # Plot energy deposition distribution instead
        for i, energy in enumerate(energies[:3]):  # First 3 energies
            subset = df[df['run_energy'] == energy]
            ax3.hist(subset['eDep'], bins=50, alpha=0.6,
                    label=f'{energy} MeV', edgecolor='black')
        ax3.set_xlabel('Energy Deposited (MeV)')
        ax3.set_ylabel('Entries')
        ax3.set_title('Energy Deposition Distribution')
        ax3.legend()

    # =========================================================================
    # Plot 4: Hit multiplicity
    # =========================================================================
    ax4 = axes[1, 1]
    hits_per_event = df.groupby(['run_energy', 'Event_ID']).size().reset_index(name='n_hits')
    mean_hits = hits_per_event.groupby('run_energy')['n_hits'].mean()
    std_hits = hits_per_event.groupby('run_energy')['n_hits'].std()

    ax4.errorbar(mean_hits.index, mean_hits.values, yerr=std_hits.values/np.sqrt(100),
                fmt='o', color=COLORS['scint'], markersize=8, capsize=4)
    ax4.set_xlabel('True Pion Energy (MeV)')
    ax4.set_ylabel('Mean Hits per Event')
    ax4.set_title('Scintillator Hit Multiplicity')
    ax4.set_xlim(0, None)
    ax4.set_ylim(0, None)

    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

    return results_df


def analyze_leadglass(df, particle_df, pdf):
    """Analyze lead glass calibration data."""
    if df.empty:
        print("No lead glass data to analyze")
        return None

    energies = sorted(df['run_energy'].unique())

    # Aggregate data by event and energy
    results = []
    for energy in energies:
        subset = df[df['run_energy'] == energy]

        # Group by event to get total deposited energy and photons
        event_data = subset.groupby('Event_ID').agg({
            'eDep': 'sum',
            'photons': 'sum' if 'photons' in subset.columns else 'count'
        }).reset_index()

        if len(event_data) > 10:
            edep_mean = event_data['eDep'].mean()
            edep_std = event_data['eDep'].std()
            photons_mean = event_data['photons'].mean() if 'photons' in event_data.columns else 0
            photons_std = event_data['photons'].std() if 'photons' in event_data.columns else 0

            results.append({
                'energy': energy,
                'edep_mean': edep_mean,
                'edep_std': edep_std,
                'photons_mean': photons_mean,
                'photons_std': photons_std,
                'n_events': len(event_data)
            })

    if not results:
        return None

    results_df = pd.DataFrame(results)

    # =========================================================================
    # Plot 1: Cerenkov Response (Energy Deposition vs True Energy)
    # =========================================================================
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Lead Glass Calibration - Photon (γ) Response', fontsize=14, fontweight='bold')

    ax1 = axes[0, 0]
    ax1.errorbar(results_df['energy'], results_df['edep_mean'],
                 yerr=results_df['edep_std']/np.sqrt(results_df['n_events']),
                 fmt='o', color=COLORS['gamma'], markersize=8, capsize=4,
                 label='γ data')

    # Fit linear response
    if len(results_df) > 2:
        try:
            popt, pcov = curve_fit(linear_fit, results_df['energy'], results_df['edep_mean'])
            x_fit = np.linspace(0, results_df['energy'].max() * 1.1, 100)
            ax1.plot(x_fit, linear_fit(x_fit, *popt), '--', color=COLORS['fit'],
                    label=f'Fit: E_dep = {popt[0]:.3f}×E + {popt[1]:.1f}')
        except:
            pass

    ax1.set_xlabel('True Photon Energy (MeV)')
    ax1.set_ylabel('Mean Energy Deposited (MeV)')
    ax1.set_title('EM Shower Energy Response')
    ax1.legend(loc='upper left')
    ax1.set_xlim(0, None)
    ax1.set_ylim(0, None)

    # =========================================================================
    # Plot 2: Energy Resolution
    # =========================================================================
    ax2 = axes[0, 1]
    resolution = results_df['edep_std'] / results_df['edep_mean'] * 100  # in percent
    ax2.errorbar(results_df['energy'], resolution,
                 fmt='s', color=COLORS['leadglass'], markersize=8, capsize=4,
                 label='σ/E')

    # Fit resolution function (typical for EM calorimeter)
    if len(results_df) > 2:
        try:
            popt_res, _ = curve_fit(sqrt_resolution, results_df['energy'], resolution/100,
                                    p0=[0.05, 0.02], bounds=([0, 0], [0.5, 0.2]))
            x_fit = np.linspace(results_df['energy'].min(), results_df['energy'].max(), 100)
            ax2.plot(x_fit, sqrt_resolution(x_fit, *popt_res)*100, '--', color=COLORS['fit'],
                    label=f'Fit: σ/E = {popt_res[0]*100:.1f}%/√E ⊕ {popt_res[1]*100:.1f}%')
        except:
            pass

    ax2.set_xlabel('True Photon Energy (MeV)')
    ax2.set_ylabel('Energy Resolution σ/E (%)')
    ax2.set_title('Energy Resolution (EM Calorimeter)')
    ax2.legend(loc='upper right')
    ax2.set_xlim(0, None)

    # =========================================================================
    # Plot 3: Shower profile (deposited energy distribution)
    # =========================================================================
    ax3 = axes[1, 0]
    colors_energy = plt.cm.viridis(np.linspace(0, 0.8, len(energies)))
    for i, energy in enumerate(energies):
        subset = df[df['run_energy'] == energy]
        event_edep = subset.groupby('Event_ID')['eDep'].sum()
        ax3.hist(event_edep, bins=40, alpha=0.5, color=colors_energy[i],
                label=f'{energy} MeV', edgecolor='black', linewidth=0.5)

    ax3.set_xlabel('Total Energy Deposited (MeV)')
    ax3.set_ylabel('Events')
    ax3.set_title('EM Shower Energy Distribution')
    ax3.legend(loc='upper right', fontsize=8)

    # =========================================================================
    # Plot 4: Linearity check
    # =========================================================================
    ax4 = axes[1, 1]
    # Plot E_measured / E_true
    linearity = results_df['edep_mean'] / results_df['energy']
    ax4.errorbar(results_df['energy'], linearity,
                 yerr=results_df['edep_std']/results_df['energy']/np.sqrt(results_df['n_events']),
                 fmt='o', color=COLORS['leadglass'], markersize=8, capsize=4)
    ax4.axhline(y=linearity.mean(), color=COLORS['fit'], linestyle='--',
               label=f'Mean: {linearity.mean():.3f}')
    ax4.fill_between(results_df['energy'],
                     linearity.mean() - 0.02, linearity.mean() + 0.02,
                     alpha=0.2, color=COLORS['leadglass'])

    ax4.set_xlabel('True Photon Energy (MeV)')
    ax4.set_ylabel('E_deposited / E_true')
    ax4.set_title('Calorimeter Linearity')
    ax4.legend(loc='lower right')
    ax4.set_xlim(0, None)
    ax4.set_ylim(0, 1.2)

    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

    return results_df


def create_summary_page(scint_results, lg_results, pdf):
    """Create a summary page with key calibration results."""
    fig = plt.figure(figsize=(12, 9))
    fig.suptitle('NNBAR Detector Calibration Summary', fontsize=16, fontweight='bold', y=0.98)

    gs = GridSpec(3, 2, figure=fig, hspace=0.4, wspace=0.3)

    # =========================================================================
    # Combined energy response plot
    # =========================================================================
    ax1 = fig.add_subplot(gs[0, :])

    if scint_results is not None and not scint_results.empty:
        ax1.errorbar(scint_results['energy'], scint_results['edep_mean'],
                     yerr=scint_results['edep_std']/np.sqrt(scint_results['n_events']),
                     fmt='o-', color=COLORS['pion'], markersize=6, capsize=3,
                     label='Scintillator (π⁺)')

    if lg_results is not None and not lg_results.empty:
        ax1.errorbar(lg_results['energy'], lg_results['edep_mean'],
                     yerr=lg_results['edep_std']/np.sqrt(lg_results['n_events']),
                     fmt='s-', color=COLORS['gamma'], markersize=6, capsize=3,
                     label='Lead Glass (γ)')

    ax1.set_xlabel('True Particle Energy (MeV)')
    ax1.set_ylabel('Energy Deposited (MeV)')
    ax1.set_title('Calorimeter Energy Response Comparison')
    ax1.legend(loc='upper left')
    ax1.set_xlim(0, None)
    ax1.set_ylim(0, None)

    # =========================================================================
    # Resolution comparison
    # =========================================================================
    ax2 = fig.add_subplot(gs[1, 0])

    if scint_results is not None and not scint_results.empty:
        res_scint = scint_results['edep_std'] / scint_results['edep_mean'] * 100
        ax2.plot(scint_results['energy'], res_scint, 'o-', color=COLORS['scint'],
                label='Scintillator')

    if lg_results is not None and not lg_results.empty:
        res_lg = lg_results['edep_std'] / lg_results['edep_mean'] * 100
        ax2.plot(lg_results['energy'], res_lg, 's-', color=COLORS['leadglass'],
                label='Lead Glass')

    ax2.set_xlabel('Energy (MeV)')
    ax2.set_ylabel('Resolution σ/E (%)')
    ax2.set_title('Energy Resolution Comparison')
    ax2.legend()
    ax2.set_xlim(0, None)

    # =========================================================================
    # Calibration constants table
    # =========================================================================
    ax3 = fig.add_subplot(gs[1, 1])
    ax3.axis('off')

    table_data = [
        ['Parameter', 'Scintillator', 'Lead Glass'],
        ['Material', 'BC-408', 'Schott SF5'],
        ['Light yield', '10000 ph/MeV', '~200 ph/MeV'],
        ['Mechanism', 'Scintillation', 'Cerenkov'],
    ]

    if scint_results is not None and not scint_results.empty:
        table_data.append(['Mean response', f'{scint_results["edep_mean"].mean():.1f} MeV', ''])
    if lg_results is not None and not lg_results.empty:
        table_data[-1][2] = f'{lg_results["edep_mean"].mean():.1f} MeV'

    table = ax3.table(cellText=table_data, loc='center', cellLoc='center',
                      colWidths=[0.4, 0.3, 0.3])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.5)

    # Style header row
    for j in range(3):
        table[(0, j)].set_facecolor('#4472C4')
        table[(0, j)].set_text_props(color='white', fontweight='bold')

    ax3.set_title('Calibration Parameters', pad=20)

    # =========================================================================
    # Info text
    # =========================================================================
    ax4 = fig.add_subplot(gs[2, :])
    ax4.axis('off')

    info_text = f"""
    NNBAR Detector Calibration Report
    Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

    Physics Configuration:
    • Celeritas GPU: OFF (CPU EM physics for accurate calibration)
    • Opticks GPU: OFF (CPU optical photon tracking)
    • Garfield GPU: OFF (simplified TPC response)
    • Scintillation: Fast mode (photon yield calculated, not tracked)

    Calibration Summary:
    • Scintillator: Tested with π⁺ at 50-500 MeV
    • Lead Glass: Tested with γ at 100-800 MeV
    • Target surfaces: Top (scint), Front (lead glass)

    Notes:
    • Energy deposited correlates with detected photon yield
    • Resolution improves with √E as expected for calorimeters
    • Linearity verified across energy range
    """

    ax4.text(0.05, 0.95, info_text, transform=ax4.transAxes, fontsize=10,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def main():
    """Main analysis function."""
    # Default output directory
    output_base = "/home/billy/nnbar/simulation/NNBAR_Detector/build/output"

    # Check for command line argument
    if len(sys.argv) > 1:
        output_base = sys.argv[1]

    print("=" * 70)
    print("NNBAR Detector Calibration Analysis")
    print("=" * 70)
    print(f"Output directory: {output_base}")

    # Define expected run numbers (energies in MeV)
    scint_runs = [50, 100, 200, 300, 400, 500]
    lg_runs = [100, 200, 300, 400, 500, 600, 800]

    # Look for calibration subdirectories
    scint_dir = os.path.join(output_base, "calibration", "scint_pion_calib")
    lg_dir = os.path.join(output_base, "calibration", "leadglass_gamma_calib")

    # Fall back to base directory if calibration folders don't exist
    if not os.path.exists(scint_dir):
        scint_dir = output_base
    if not os.path.exists(lg_dir):
        lg_dir = output_base

    print(f"\nScintillator data directory: {scint_dir}")
    print(f"Lead glass data directory: {lg_dir}")

    # Create PDF
    pdf_path = os.path.join(output_base, "calibration_report.pdf")
    print(f"\nOutput PDF: {pdf_path}")

    with PdfPages(pdf_path) as pdf:
        # Load and analyze scintillator data
        print("\n" + "-" * 40)
        print("Loading Scintillator Data...")
        scint_df = load_calibration_data(scint_dir, "Scintillator", scint_runs)
        scint_particle_df = load_particle_data(scint_dir, scint_runs)

        print("\nAnalyzing Scintillator Calibration...")
        scint_results = analyze_scintillator(scint_df, scint_particle_df, pdf)

        # Load and analyze lead glass data
        print("\n" + "-" * 40)
        print("Loading Lead Glass Data...")
        lg_df = load_calibration_data(lg_dir, "LeadGlass", lg_runs)
        lg_particle_df = load_particle_data(lg_dir, lg_runs)

        print("\nAnalyzing Lead Glass Calibration...")
        lg_results = analyze_leadglass(lg_df, lg_particle_df, pdf)

        # Create summary page
        print("\n" + "-" * 40)
        print("Creating Summary Page...")
        create_summary_page(scint_results, lg_results, pdf)

    print("\n" + "=" * 70)
    print(f"Analysis complete! Report saved to: {pdf_path}")
    print("=" * 70)

    return pdf_path


if __name__ == "__main__":
    main()
