# ================================================================================
# SCRIPT 5: JWST HIGH-Z DIPOL-TEST (Mock-Daten)
# ================================================================================
print("SCRIPT 5: JWST High-z Dipol-Test (jwsthighzdipoletest.py)")
print("-"*80)

# Mock JWST SNIa Daten bei z=0.5-1.0
np.random.seed(42)
n_samples = 500

# Himmelspositionen (RA, Dec)
ra = np.random.uniform(0, 360, n_samples)
dec = np.random.uniform(-90, 90, n_samples)

# Shapley-Dipol-Achse (RA=170°, Dec=7°)
ra_dipol, dec_dipol = 170, 7

# Winkel zur Dipol-Achse
theta = np.arccos(
    np.sin(np.radians(dec)) * np.sin(np.radians(dec_dipol)) +
    np.cos(np.radians(dec)) * np.cos(np.radians(dec_dipol)) *
    np.cos(np.radians(ra - ra_dipol))
)
theta_deg = np.degrees(theta)

# Modell: ΔH(θ) = A·cos(θ) mit A=2.5 km/s/Mpc (erwartete Amplitude)
A_true = 2.5  # km/s/Mpc
H0_iso_true = 70.8  # km/s/Mpc
sigma_intrinsic = 1.5  # km/s/Mpc (Messunsicherheit)

H0_obs = H0_iso_true + A_true * np.cos(theta) + np.random.normal(0, sigma_intrinsic, n_samples)

# Linearer Fit: H₀ = H_iso + A·cos(θ)
X = np.column_stack([np.ones(n_samples), np.cos(theta)])
params = np.linalg.lstsq(X, H0_obs, rcond=None)[0]
H_iso_fit, A_fit = params

# Fehler-Schätzung
residuals = H0_obs - (H_iso_fit + A_fit * np.cos(theta))
sigma_fit = np.std(residuals)
cov = sigma_fit**2 * np.linalg.inv(X.T @ X)
A_err = np.sqrt(cov[1,1])

# Signifikanz
significance = A_fit / A_err

print(f"Mock-Daten: N={n_samples} SNIa @ z=0.5-1.0")
print(f"Input: H_iso={H0_iso_true} km/s/Mpc, A={A_true} km/s/Mpc")
print(f"\nFit-Resultate:")
print(f"  H_iso(fit) = {H_iso_fit:.2f} ± {np.sqrt(cov[0,0]):.2f} km/s/Mpc")
print(f"  A(fit) = {A_fit:.2f} ± {A_err:.2f} km/s/Mpc")
print(f"  Signifikanz: {significance:.1f}σ")
print(f"  Residuen σ: {sigma_fit:.2f} km/s/Mpc")

print(f"\nFalsifikations-Kriterium:")
print(f"  Falls A < 2.0σ bei JWST-Daten → Mode E widerlegt")
print(f"  Falls θ_max ≠ 71° ± 5° → Dipol-Achse inkonsistent")

results_script5 = {
    'N_samples': n_samples,
    'H_iso_fit_kms_Mpc': round(H_iso_fit, 2),
    'A_fit_kms_Mpc': round(A_fit, 2),
    'A_err_kms_Mpc': round(A_err, 2),
    'significance_sigma': round(significance, 1),
    'falsified': significance < 2.0
}

print("\n" + "="*80 + "\n")