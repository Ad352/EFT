# ================================================================================
# SCRIPT 2: HARDY-PROJEKTOR ORTHOGONALITÄT (Temporal)
# ================================================================================
print("SCRIPT 2: Hardy-Raum Orthogonalität (H⁺ ⊥ H⁻)")
print("-"*80)

# Simuliere Moden im Frequenzraum
N = 1000
omega = np.linspace(-10, 10, N)
domega = omega[1] - omega[0]

# Mode A+D (retardiert, H⁺): nur ω>0
H_forward = np.where(omega > 0, np.exp(-omega**2/10), 0)

# Mode B+C (avanciert, H⁻): nur ω<0
H_backward = np.where(omega < 0, np.exp(-omega**2/10), 0)

# Skalarprodukt im Frequenzraum
inner_product = np.sum(H_forward * H_backward) * domega

print(f"H⁺(ω>0) Support: {np.sum(H_forward>0)} Punkte")
print(f"H⁻(ω<0) Support: {np.sum(H_backward>0)} Punkte")
print(f"Skalarprodukt ⟨H⁺|H⁻⟩_F = {inner_product:.2e}")
print(f"  → Orthogonal: {abs(inner_product) < 1e-10}")

# Norm der Einzelmoden
norm_forward = np.sqrt(np.sum(H_forward**2) * domega)
norm_backward = np.sqrt(np.sum(H_backward**2) * domega)

print(f"\nNormen:")
print(f"  ||H⁺|| = {norm_forward:.4f}")
print(f"  ||H⁻|| = {norm_backward:.4f}")
print(f"  Kosinus-Winkel: {inner_product/(norm_forward*norm_backward):.2e}")

results_script2 = {
    'inner_product': float(f"{inner_product:.2e}"),
    'orthogonal': abs(inner_product) < 1e-10,
    'norm_forward': round(norm_forward, 4),
    'norm_backward': round(norm_backward, 4)
}

print("\n" + "="*80 + "\n")