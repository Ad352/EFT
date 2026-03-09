import numpy as np
import pandas as pd
import json

# ================================================================================
# J.8.1 ANHANG: VOLLSTÄNDIGE TEST-SCRIPTS FÜR EFT-LAPSE ORTHOGONALITÄT
# ================================================================================
# Version: 5.1 (März 2026)
# Zweck: Reproduzierbare Validierung der Hardy-Projektor-Orthogonalität
# ================================================================================

print("="*80)
print("EFT-LAPSE THEORY V5.1 - ADDENDUM J.8.1")
print("Vollständige Test-Scripts für Orthogonalitäts-Validierung")
print("="*80 + "\n")

# ================================================================================
# SCRIPT 1: BULK-FLOW-GRADIENT BERECHNUNG (Mode E Parameter)
# ================================================================================
print("SCRIPT 1: Bulk-Flow-Gradient aus CosmicFlows-4 Daten")
print("-"*80)

# Datenquelle: CosmicFlows-4 (Tully et al. 2023, MNRAS 524, 1857)
# URL: https://academic.oup.com/mnras/article/524/2/1857/7218572
bulkflow_data = {
    'Radius_h1Mpc': [150, 200],
    'Velocity_kms': [390, 419],
    'Sigma_kms': [78, 36],
    'Reference': ['Tully Fig.4 (2023)', 'Watkins Fig.4 (2007)']
}
df_bulk = pd.DataFrame(bulkflow_data)

# EFT-Lapse Parameter
f_growth = 0.45  # Wachstumsfaktor Ωm^0.55
H0_bare = 67.4   # km/s/Mpc (Planck 2018)
alpha_cosmo = 1/(4*np.pi)  # 0.0796

# Gradient des Gravitationspotentials
nabla_Phi = (4 * np.pi * df_bulk['Velocity_kms']) / (f_growth * H0_bare * df_bulk['Radius_h1Mpc'])

# Mode E Amplitude Ψ₀
Psi0_mean = alpha_cosmo * np.mean(nabla_Phi)
Psi0_relerr = np.sqrt(np.mean((df_bulk['Sigma_kms']/df_bulk['Velocity_kms'])**2))
Psi0_err = Psi0_mean * Psi0_relerr

# Hubble-Modulation ΔH_E
Delta_HE_max = Psi0_mean * 100  # Max bei θ=0°
Delta_HE_err = Psi0_err * 100
Delta_HE_MW = Delta_HE_max * np.cos(np.radians(71))  # MW @ 71° zu Shapley

print("\nCF4 Bulk Flow Daten:")
print(df_bulk.to_string(index=False))
print(f"\nGradient ∇Φ [1/Mpc]: {nabla_Phi.values}")
print(f"Ψ₀ = {Psi0_mean:.4f} ± {Psi0_err:.4f}")
print(f"\nMode E Vorhersagen:")
print(f"  ΔH_E(max) = {Delta_HE_max:.2f} ± {Delta_HE_err:.2f} km/s/Mpc")
print(f"  ΔH_E(MW@71°) = {Delta_HE_MW:.2f} km/s/Mpc")
print(f"  Falsifikation: JWST sollte 71° ± 5° messen\n")

results_script1 = {
    'Psi0': round(Psi0_mean, 4),
    'Psi0_err': round(Psi0_err, 4),
    'Delta_HE_max_kms_Mpc': round(Delta_HE_max, 2),
    'Delta_HE_MW_kms_Mpc': round(Delta_HE_MW, 2),
    'MW_angle_deg': 71
}

print("="*80 + "\n")