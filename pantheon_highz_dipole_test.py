"""
EFT-Lapse Theory: Empirical High-z Dipole Test on Pantheon+ Data
File: pantheon_highz_dipole_test.py
Author: LLM-EFT-Lapse Collaboration (February 2026)

This script applies the JWST "Residual-mu" methodology (Mode E Dipole)
to the existing Pantheon+SH0ES dataset for supernovae at z > 0.5.
It tests whether the spatial ∇Ψ gradient observed in the local universe
persists into the deep Hubble flow.
"""

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import sys
import os

# --- CONFIGURATION ---
DATA_FILE = "PantheonSH0ES.dat"
COLD_SPOT_RA = 55.0   # Eridanus Supervoid (Minimum H0 expected)
COLD_SPOT_DEC = -25.0
M_ABS = -19.2435      # SH0ES calibration baseline
OMEGA_M_REF = 0.315   # Reference background cosmology

def angular_separation(ra1, dec1, ra2, dec2):
    """Haversine formula for angular distance on celestial sphere."""
    ra1, dec1 = np.radians(ra1), np.radians(dec1)
    ra2, dec2 = np.radians(ra2), np.radians(dec2)
    dlon = ra2 - ra1
    dlat = dec2 - dec1
    a = np.sin(dlat/2)**2 + np.cos(dec1) * np.cos(dec2) * np.sin(dlon/2)**2
    return np.degrees(2 * np.arcsin(np.sqrt(np.clip(a, 0, 1))))

def theoretical_distance_modulus(z, h0, omega_m):
    """Calculates theoretical distance modulus for a flat LCDM background."""
    from scipy.integrate import quad
    def integrand(x):
        return 1.0 / np.sqrt(omega_m * (1+x)**3 + (1-omega_m))

    mu_theoretical = []
    for z_i in z:
        integral, _ = quad(integrand, 0, z_i)
        # Luminosity distance in Mpc
        dL = (299792.458 * (1 + z_i) / h0) * integral
        mu_theoretical.append(5 * np.log10(dL) + 25)
    return np.array(mu_theoretical)

def dipole_model(theta, h_iso, a_dip):
    """Cosine dipole model: H_eff(theta) = H_iso + A_dip * cos(theta)"""
    return h_iso + a_dip * np.cos(np.radians(theta))

def load_pantheon_data(filepath):
    """Loads and standardizes Pantheon+ column names."""
    if not os.path.exists(filepath):
        print(f"ERROR: Dataset '{filepath}' not found. Please ensure it is in the same directory.")
        sys.exit(1)

    try:
        df = pd.read_csv(filepath, sep=r'\s+')
        if 'zHD' not in df.columns:
            df = pd.read_csv(filepath, sep=',')
    except Exception as e:
        print(f"ERROR reading file: {e}")
        sys.exit(1)

    # Standardize column names
    col_map = {
        'ZHD': 'zHD', 'Z_HD': 'zHD', 
        'MBCORR': 'mbcorr', 'MB_CORR': 'mbcorr', 'mb': 'mbcorr',
        'RA': 'RA', 'DEC': 'DEC', 'DECL': 'DEC'
    }
    df.columns = df.columns.str.upper() # make upper for matching
    # re-map back to exact names needed
    df = df.rename(columns={k: v for k, v in col_map.items() if k in df.columns})

    # Ensure required columns exist
    required = ['zHD', 'mbcorr', 'RA', 'DEC']
    missing = [col for col in required if col not in df.columns]
    if missing:
        print(f"ERROR: Missing required columns in dataset: {missing}")
        print(f"Available columns: {list(df.columns)}")
        sys.exit(1)

    return df

def analyze_dipole_in_z_bin(df, z_min, z_max, bin_name, h0_ref=73.0):
    """
    Extracts effective H0 from mu-residuals and fits the dipole amplitude.
    Uses h0_ref as the baseline for linearization.
    """
    subset = df[(df['zHD'] >= z_min) & (df['zHD'] < z_max)].copy()
    if len(subset) < 15:
        print(f"--- Epoch: {bin_name} (z={z_min} to {z_max}) ---")
        print(f"  Insufficient data (N={len(subset)}). Skipping.\n")
        return None

    # 1. Observed distance modulus from Pantheon
    subset['MU_OBS'] = subset['mbcorr'] - M_ABS

    # 2. Compute angular distance to Cold Spot (Dipole minimum)
    theta = angular_separation(subset['RA'].values, subset['DEC'].values, 
                               COLD_SPOT_RA, COLD_SPOT_DEC)

    # 3. Extract Effective H0 from residuals
    # theoretical mu assuming isotropic baseline h0_ref
    mu_bg = theoretical_distance_modulus(subset['zHD'].values, h0=h0_ref, omega_m=OMEGA_M_REF)
    delta_mu = subset['MU_OBS'].values - mu_bg

    # Linearization: Delta mu ≈ - (5 / ln(10)) * (Delta H0 / H0)
    # => H_eff = H0_ref * (1 - (ln(10) / 5) * delta_mu)
    h0_effective = h0_ref * (1 - (np.log(10) / 5) * delta_mu)

    # 4. Fit Dipole Model
    try:
        popt, pcov = curve_fit(dipole_model, theta, h0_effective, p0=[h0_ref, -2.0])
        h_iso, a_dip = popt
        err_h_iso, err_a_dip = np.sqrt(np.diag(pcov))
    except RuntimeError:
        print(f"--- Epoch: {bin_name} ---")
        print("  Curve fit failed. Data might be too sparse/noisy.\n")
        return None

    # 5. Compute Pearson correlation
    corr = np.corrcoef(np.cos(np.radians(theta)), h0_effective)[0, 1]

    print(f"--- Epoch: {bin_name} (z={z_min} to {z_max}) | N={len(subset)} ---")
    print(f"  Isotropic Baseline (H_iso) : {h_iso:.2f} ± {err_h_iso:.2f} km/s/Mpc")
    print(f"  Dipole Amplitude (A_dip)   : {abs(a_dip):.2f} ± {err_a_dip:.2f} km/s/Mpc")
    print(f"  Correlation r(H0, cos(th)) : {corr:.3f}")

    # Significance Check
    significance = abs(a_dip) / err_a_dip if err_a_dip > 0 else 0
    if significance > 2.0 and corr > 0.2:
        print(f"  -> RESULT: Mode E Dipole DETECTED ({significance:.1f}σ).")
    else:
        print("  -> RESULT: Dipole NOT statistically significant in this bin.")
    print("")
    return a_dip

def main():
    print("===================================================================")
    print(" EFT-LAPSE V5.1: HIGH-Z DIPOLE TEST ON PANTHEON+ (Residual Method) ")
    print("===================================================================\n")

    df = load_pantheon_data(DATA_FILE)

    print(f"Dataset loaded: {len(df)} total supernovae.")
    print(f"Testing Dipole aligned with Eridanus Cold Spot (RA={COLD_SPOT_RA}, DEC={COLD_SPOT_DEC}).\n")

    # Test 1: Low-z Anchor (Sanity Check vs anisotropytest.py)
    # Uses Residual method on low-z to verify it matches direct H0=cz/d calculation
    analyze_dipole_in_z_bin(df, 0.015, 0.15, "Low-z Anchor (Local Voids)", h0_ref=73.0)

    # Test 2: Intermediate-z (Transition into Hubble Flow)
    analyze_dipole_in_z_bin(df, 0.15, 0.5, "Intermediate-z", h0_ref=70.0)

    # Test 3: High-z (Deep Hubble Flow / Dark Energy dominated)
    analyze_dipole_in_z_bin(df, 0.5, 2.3, "High-z (Pantheon Deep)", h0_ref=70.0)

if __name__ == "__main__":
    main()
