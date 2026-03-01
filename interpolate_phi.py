import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
from astropy.coordinates import SkyCoord
import astropy.units as u

# 1. Daten laden
print("Lade Felder und SNe-Daten...")
phi = np.load('phi_field.npy')
sne = pd.read_csv('my_sne.csv')

# 2. Gitter-Koordinaten definieren (wie im README: -200 bis +200 Mpc/h)
N = phi.shape[0]          # 257
L_box = 400.0             # Gesamte Kantenlänge
# Erzeugt 257 Punkte von -200 bis +200
coords = np.linspace(-L_box/2, L_box/2, N)

# Interpolator erstellen (ermöglicht Abfrage an exakten Koordinaten)
# bounds_error=False setzt Werte außerhalb der Box auf NaN
interp = RegularGridInterpolator((coords, coords, coords), phi, 
                                 bounds_error=False, fill_value=np.nan)

print("Starte Umrechnung der Koordinaten (RA/Dec/z -> SGX/SGY/SGZ)...")

def get_phi_for_sn(row):
    # a) Distanzschätzung (H0=70 km/s/Mpc)
    dist_mpc_h = (299792.458 * row['z']) / 70.0 
    
    # b) Transformation in Supergalaktische Koordinaten
    c = SkyCoord(ra=row['ra']*u.degree, 
                 dec=row['dec']*u.degree, 
                 distance=dist_mpc_h*u.Mpc, 
                 frame='icrs')
    
    # Umwandlung in das Supergalaktische System
    sg = c.supergalactic
    
    # WICHTIG: Explizit die kartesischen Werte (x, y, z) abgreifen
    # Diese entsprechen SGX, SGY, SGZ
    point = [sg.cartesian.x.value, 
             sg.cartesian.y.value, 
             sg.cartesian.z.value]
    
    # c) Wert aus dem Potentialfeld abfragen
    val = interp(point)[0]
    return val

# 3. Für alle SNe berechnen
sne['phi_lokal'] = sne.apply(get_phi_for_sn, axis=1)

# 4. Statistik anzeigen
valid_sne = sne.dropna(subset=['phi_lokal'])
print("-" * 30)
print(f"Verarbeitung abgeschlossen!")
print(f"SNe im Datensatz: {len(sne)}")
print(f"SNe innerhalb der Box gefunden: {len(valid_sne)}")
print(f"Phi-Bereich der SNe: {valid_sne['phi_lokal'].min():.2f} bis {valid_sne['phi_lokal'].max():.2f}")

# 5. Speichern
sne.to_csv('sne_mit_phi.csv', index=False)
print("Ergebnis in 'sne_mit_phi.csv' gespeichert.")
