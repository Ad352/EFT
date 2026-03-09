
import pandas as pd

# ================================================================================
# Script-Übersicht für J.8.1
# ================================================================================

print("KORRIGIERTE SCRIPT-ÜBERSICHT FÜR J.8.1 ANHANG")
print("="*100)

# Daten für die Tabelle
data = {
    'Script': [
        '1. Bulk-Flow-Gradient',
        '2. Hardy-Orthogonalität',
        '3. Angular-Orthogonalität',
        '4. Hybrid-Summenregel',
        '5. JWST Mock-Test'
    ],
    'Zweck': [
        'Ψ₀-Berechnung aus CosmicFlows-4',
        'Temporal H⁺⊥H⁻ Test',
        'Spherical Y₀⁰⊥Y₁⁰ Test',
        'Moden A-E Addition',
        'High-z Dipol-Detektion'
    ],
    'Kernresultat': [
        'Ψ₀=0.0774±0.0119, ΔH_E(max)=5.22 km/s/Mpc, ΔH_E(MW)=1.70 km/s/Mpc',
        '⟨H⁺|H⁻⟩=0.0 (exakt orthogonal) ✓',
        '∫Y₀⁰Y₁⁰ dΩ=-7.4×10⁻¹⁷ (orthogonal) ✓',
        'H₀(obs,MW)=72.52 km/s/Mpc | SH0ES: 73.04±1.04 → Δ=0.52 km/s/Mpc (0.5σ) ✓',
        'A_fit=2.59±0.13 km/s/Mpc, Signifikanz: 20.2σ'
    ]
}

df = pd.DataFrame(data)

# Ausgabe als formatierte Tabelle
print("\n" + df.to_string(index=False))

print("\n" + "="*100)
print("\nÄNDERUNGEN IN SCRIPT 4:")
print("-"*100)
print("ALT (fehlerhaft):")
print("  H₀(obs)=71.47 km/s/Mpc (SH0ES: 73.04)")
print("  → Differenz: 1.57 km/s/Mpc")
print("")
print("NEU (korrigiert):")
print("  H₀(obs,MW)=72.52 km/s/Mpc | SH0ES: 73.04±1.04 → Δ=0.52 km/s/Mpc (0.5σ) ✓")
print("  → Differenz: 0.52 km/s/Mpc (innerhalb 1σ)")

print("\n" + "="*100)
print("\nMARKDOWN-FORMAT FÜR COPY-PASTE:")
print("-"*100)

# Markdown-Tabelle erstellen
print("\n| Script | Zweck | Kernresultat |")
print("|--------|-------|--------------|")
for idx, row in df.iterrows():
    print(f"| **{row['Script']}** | {row['Zweck']} | {row['Kernresultat']} |")

print("\n" + "="*100)