import numpy as np

# 1. Daten laden (falls nicht noch im Speicher)
delta = np.load('twompp_density.npy')
N = delta.shape[0]        # 257
L_box = 400.0             # Kantenlänge der Box in Mpc/h (BITTE PRÜFEN!)

print("Starte FFT...")

# 2. Dichte in den Fourier-Raum transformieren
# Wir nutzen rfftn (real-input FFT), das spart 50% Speicher
delta_k = np.fft.rfftn(delta)

# 3. Wellenzahlen k für jede Dimension erstellen
# k = 2*pi * n / L
k_per_step = 2 * np.pi / L_box
kx = np.fft.fftfreq(N) * N * k_per_step
ky = np.fft.fftfreq(N) * N * k_per_step
kz = np.fft.rfftfreq(N) * N * k_per_step  # rfftfreq für die letzte Dimension

# 3D-Gitter der k-Werte erstellen (Broadcasting für Speicher-Effizienz)
# k^2 = kx^2 + ky^2 + kz^2
ky_grid, kx_grid, kz_grid = np.meshgrid(ky, kx, kz, indexing='ij') 
k_sq = kx_grid**2 + ky_grid**2 + kz_grid**2

# 4. Poisson-Gleichung im k-Raum lösen: Phi_k = - delta_k / k^2
# Achtung: Bei k=0 ist k^2=0 (Division durch Null).
# Wir setzen den DC-Mode (Mittelwert) manuell auf 0.
k_sq[0,0,0] = 1.0  # Temporär, um Warnung zu vermeiden
potential_k = - delta_k / k_sq
potential_k[0,0,0] = 0.0 # Den "Null-Punkt" des Potentials fixieren

# 5. Zurück in den Realraum transformieren
phi = np.fft.irfftn(potential_k, s=delta.shape)

# Optional: Skalierung auf physikalische Einheiten (ohne G/c^2 Faktoren)
# Das Ergebnis ist jetzt proportional zu Phi. Für Korrelation reicht das.
# Physikalisch wäre da noch ein Faktor: 1.5 * Omega_m * (H0/c)^2 ...

print("Fertig!")
print(f"Phi Min: {phi.min():.4e}")
print(f"Phi Max: {phi.max():.4e}")

# 6. Speichern für den nächsten Schritt (Interpolation)
np.save('phi_field.npy', phi)
print("Potential als 'phi_field.npy' gespeichert.")
