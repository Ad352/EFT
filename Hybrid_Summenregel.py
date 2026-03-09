# Mode E Parameter (aus Script 1: Bulk-Flow-Gradient)
Psi0 = 0.0774  # Gradient-Amplitude (dimensionslos)
Psi0_err = 0.0119
theta_MW = 71  # Grad (Milchstraßen-Position relativ zu Shapley)

# ================================================================================
# BERECHNUNG: ISOTROPER BLOCK (MODEN A-D)
# ================================================================================
print("\n1. ISOTROPER BLOCK (Monopol L=0, richtungsgemittelt)")
print("-"*80)

# Deterministische Beiträge (linear addiert, H⁺-Sektor)
H_det_iso = modes['A']['value'] + modes['D']['value']
print(f"Deterministisch (A+D):")
print(f"  Mode A (Lapse): {modes['A']['value']:.2f} km/s/Mpc")
print(f"  Mode D (Wisconsin): {modes['D']['value']:.2f} km/s/Mpc")
print(f"  → H_det = {H_det_iso:.2f} km/s/Mpc (linear)")

# Stochastische Beiträge (quadratisch addiert wegen Jensen-Bias, H⁻-Sektor)
H_stoch_variance = modes['B']['value']**2 + modes['C']['value']**2
H_stoch = np.sqrt(H_stoch_variance)
print(f"\nStochastisch (√(B²+C²)):")
print(f"  Mode B (Backreaction): {modes['B']['value']:.2f} km/s/Mpc")
print(f"  Mode C (Void-Selektion): {modes['C']['value']:.2f} km/s/Mpc")
print(f"  → Varianz = {modes['B']['value']:.2f}² + {modes['C']['value']:.2f}² = {H_stoch_variance:.4f}")
print(f"  → H_stoch = √{H_stoch_variance:.4f} = {H_stoch:.2f} km/s/Mpc (pythagoreisch)")

# Isotrope Gesamtverschiebung
H_iso = H_det_iso + H_stoch
print(f"\nIsotrope Gesamtverschiebung:")
print(f"  ΔH_iso = {H_det_iso:.2f} + {H_stoch:.2f} = {H_iso:.2f} km/s/Mpc")

# Isotrope Hubble-Konstante
H0_iso = H0_Planck + H_iso
print(f"  H₀(iso) = {H0_Planck:.1f} + {H_iso:.2f} = {H0_iso:.2f} km/s/Mpc")
print(f"\n  ✓ Vergleich mit TRGB (Freedman et al. 2021): 69.8 ± 1.7 km/s/Mpc")
print(f"    Differenz: {abs(H0_iso - 69.8):.2f} km/s/Mpc (innerhalb 1σ)")

# ================================================================================
# BERECHNUNG: ANISOTROPER BLOCK (MODE E)
# ================================================================================
print("\n" + "="*80)
print("2. ANISOTROPER BLOCK (Dipol L=1, richtungsabhängig)")
print("-"*80)

# Mode E Maximum bei θ=0° (Shapley-Attraktor-Zentrum)
Delta_HE_max = Psi0 * H0_Planck
Delta_HE_max_err = Psi0_err * H0_Planck

print(f"Dipol-Gradient (Mode E):")
print(f"  Ψ₀ = {Psi0:.4f} ± {Psi0_err:.4f} (aus CF4 Bulk Flow)")
print(f"  ΔH_E(max, θ=0°) = Ψ₀ × H₀_Planck")
print(f"                  = {Psi0:.4f} × {H0_Planck:.1f}")
print(f"                  = {Delta_HE_max:.2f} ± {Delta_HE_max_err:.2f} km/s/Mpc")

# Projektion auf Milchstraßen-Position
cos_MW = np.cos(np.radians(theta_MW))
Delta_HE_MW = Delta_HE_max * cos_MW
Delta_HE_MW_err = Delta_HE_max_err * cos_MW

print(f"\nProjektion auf MW-Position (θ_MW = {theta_MW}°):")
print(f"  cos({theta_MW}°) = {cos_MW:.4f}")
print(f"  ΔH_E(MW) = {Delta_HE_max:.2f} × {cos_MW:.4f}")
print(f"           = {Delta_HE_MW:.2f} ± {Delta_HE_MW_err:.2f} km/s/Mpc")

# Update Mode E mit berechnetem Wert
modes['E']['value'] = Delta_HE_MW

print(f"\n  Physikalische Bedeutung:")
print(f"  • MW liegt nicht IM Shapley-Attraktor (θ=0°)")
print(f"  • Sondern 71° entfernt → Kosinus-Dämpfung")
print(f"  • Maximum bei θ=0°: {Delta_HE_max:.2f} km/s/Mpc")
print(f"  • Lokal bei MW: {Delta_HE_MW:.2f} km/s/Mpc")

# ================================================================================
# FINALE BEOBACHTETE HUBBLE-KONSTANTE
# ================================================================================
print("\n" + "="*80)
print("3. FINALE BEOBACHTETE HUBBLE-KONSTANTE (H₀_obs,MW)")
print("-"*80)

H0_obs_MW = H0_iso + Delta_HE_MW
H0_obs_MW_err = np.sqrt(Delta_HE_MW_err**2 + 0.5**2)  # Systematische Unsicherheit

print(f"Vollständige Summenregel:")
print(f"  H₀(obs,MW) = H₀(bare) + ΔH_iso + ΔH_E(MW)")
print(f"             = H₀_Planck + [A+D + √(B²+C²)] + [Ψ₀·H₀·cos(θ_MW)]")
print(f"             = {H0_Planck:.1f} + {H_iso:.2f} + {Delta_HE_MW:.2f}")
print(f"             = {H0_obs_MW:.2f} ± {H0_obs_MW_err:.2f} km/s/Mpc")

# ================================================================================
# VALIDIERUNG GEGEN BEOBACHTUNGEN
# ================================================================================
print("\n" + "="*80)
print("4. VALIDIERUNG GEGEN BEOBACHTUNGEN")
print("-"*80)

observations = {
    'Planck 2018': {'value': 67.4, 'error': 0.5, 'type': 'CMB (volumengemittelt)'},
    'TRGB 2021': {'value': 69.8, 'error': 1.7, 'type': 'Lokale SNIa (isotrop gemittelt)'},
    'SH0ES 2022': {'value': 73.04, 'error': 1.04, 'type': 'Cepheids (Shapley-korreliert)'}
}

print(f"{'Messung':<15} {'Wert [km/s/Mpc]':<20} {'EFT-Vorhersage':<20} {'Δ [km/s/Mpc]':<15} {'Signifikanz'}")
print("-"*95)

for name, obs in observations.items():
    if 'Planck' in name:
        eft_prediction = H0_Planck
        description = "H₀(bare)"
    elif 'TRGB' in name:
        eft_prediction = H0_iso
        description = "H₀(iso)"
    else:  # SH0ES
        eft_prediction = H0_obs_MW
        description = "H₀(obs,MW)"
    
    delta = abs(obs['value'] - eft_prediction)
    sigma = delta / obs['error']
    
    status = "✓" if sigma < 1.0 else ("~" if sigma < 2.0 else "✗")
    print(f"{name:<15} {obs['value']:.2f} ± {obs['error']:.2f}        {eft_prediction:.2f} ({description:<12}) {delta:>6.2f}          {sigma:.2f}σ {status}")

print("\n" + "-"*95)
print("Legende: ✓ < 1σ  |  ~ 1-2σ  |  ✗ > 2σ")

# ================================================================================
# MODEN-KLASSIFIKATIONS-TABELLE
# ================================================================================
print("\n" + "="*80)
print("5. VOLLSTÄNDIGE MODEN-KLASSIFIKATION")
print("-"*80)

print(f"{'Mode':<6} {'Beschreibung':<25} {'ΔH [km/s/Mpc]':<15} {'Typ':<15} {'Hardy':<6} {'L':<3}")
print("-"*95)
for mode_name, params in modes.items():
    value_str = f"{params['value']:.2f}" if params['value'] is not None else "berechnet"
    print(f"{mode_name:<6} {params['description']:<25} {value_str:<15} {params['type']:<15} {params['Hardy']:<6} {params['L']:<3}")

# ================================================================================
# FALSIFIKATIONS-KRITERIEN
# ================================================================================
print("\n" + "="*80)
print("6. FALSIFIKATIONS-KRITERIEN")
print("-"*80)
print(f"Die EFT-Lapse-Theorie wird widerlegt, falls:")
print(f"")
print(f"1. JWST/Euclid High-z SNIa (2027+):")
print(f"   → Keine Dipol-Signatur in H₀(θ) mit A < 2.0σ")
print(f"   → Erwartung: A = {Delta_HE_max/2:.2f} ± 0.5 km/s/Mpc bei z=0.5-1.0")
print(f"")
print(f"2. Dipol-Achsen-Inkonsistenz:")
print(f"   → θ_max ≠ 71° ± 10° (Shapley/CMB-Dipol)")
print(f"   → Erwartung: RA ≈ 170°, Dec ≈ 7° (Shapley-Attraktor)")
print(f"")
print(f"3. TRGB vs. SH0ES Tension verschwindet:")
print(f"   → Falls neue Kalibrierungen beide zu H₀ ≈ 67-68 km/s/Mpc führen")
print(f"   → Mode E würde dann nicht mehr benötigt")

print("\n" + "="*80)
print("SCRIPT 4 ABGESCHLOSSEN")
print("Alle Berechnungen folgen EFT-Lapse V5.1 (Februar 2026)")
print("="*80)

# ================================================================================
# EXPORT FÜR WEITERE VERWENDUNG
# ================================================================================
results_script4 = {
    'H_iso_kms_Mpc': round(H_iso, 2),
    'H0_iso_kms_Mpc': round(H0_iso, 2),
    'Delta_HE_max_kms_Mpc': round(Delta_HE_max, 2),
    'Delta_HE_MW_kms_Mpc': round(Delta_HE_MW, 2),
    'H0_obs_MW_kms_Mpc': round(H0_obs_MW, 2),
    'theta_MW_deg': theta_MW,
    'Psi0': Psi0,
    'residual_SH0ES_kms_Mpc': round(abs(H0_obs_MW - 73.04), 2),
    'residual_SH0ES_sigma': round(abs(H0_obs_MW - 73.04) / 1.04, 2)
}

print(f"\nJSON-Export:")
import json
print(json.dumps(results_script4, indent=2))