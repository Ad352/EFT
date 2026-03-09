# ================================================================================
# SCRIPT 3: ANGULAR PROJEKTOR ORTHOGONALITÄT (L=0 ⊥ L=1)
# ================================================================================
print("SCRIPT 3: Spherical Harmonic Orthogonalität (Y₀⁰ ⊥ Y₁⁰)")
print("-"*80)

# Diskretisierung der Himmelssphäre
n_theta = 200
n_phi = 400
theta = np.linspace(0, np.pi, n_theta)
phi = np.linspace(0, 2*np.pi, n_phi)
THETA, PHI = np.meshgrid(theta, phi, indexing='ij')

# Kugelflächenfunktionen
Y00 = 1/np.sqrt(4*np.pi) * np.ones_like(THETA)  # Monopol (L=0)
Y10 = np.sqrt(3/(4*np.pi)) * np.cos(THETA)      # Dipol (L=1)

# Integration über Sphäre: ∫ Y₀⁰ Y₁⁰ sin(θ) dθ dφ
integrand = Y00 * Y10 * np.sin(THETA)
dtheta = theta[1] - theta[0]
dphi = phi[1] - phi[0]
integral_L0_L1 = np.sum(integrand) * dtheta * dphi

print(f"Diskretisierung: {n_theta}×{n_phi} = {n_theta*n_phi} Punkte")
print(f"∫_S² Y₀⁰ Y₁⁰ dΩ = {integral_L0_L1:.2e}")
print(f"  → Orthogonal: {abs(integral_L0_L1) < 1e-10}")

# Normierungschecks
norm_Y00 = np.sqrt(np.sum(Y00**2 * np.sin(THETA)) * dtheta * dphi)
norm_Y10 = np.sqrt(np.sum(Y10**2 * np.sin(THETA)) * dtheta * dphi)

print(f"\nNormierungen:")
print(f"  ⟨Y₀⁰|Y₀⁰⟩ = {norm_Y00**2:.4f} (soll: 1.0000)")
print(f"  ⟨Y₁⁰|Y₁⁰⟩ = {norm_Y10**2:.4f} (soll: 1.0000)")

results_script3 = {
    'integral_L0_L1': float(f"{integral_L0_L1:.2e}"),
    'orthogonal': abs(integral_L0_L1) < 1e-10,
    'norm_Y00_squared': round(norm_Y00**2, 4),
    'norm_Y10_squared': round(norm_Y10**2, 4)
}

print("\n" + "="*80 + "\n")