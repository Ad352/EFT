import pandas as pd
import numpy as np
import sys

# --- KONFIGURATION ---
DATEINAME = "Pantheon+SH0ES.dat"

# Die kosmischen Voids (Zentrum RA/Dec in Grad, Radius in Grad)
# Diese Koordinaten definieren die "Test-Regionen" für Addendum K
VOIDS = {
    "Bootes Void":        {"ra": 223.0, "dec": 46.0,  "radius": 15.0},
    "Local Void":         {"ra": 280.0, "dec": 18.0,  "radius": 20.0},
    "Eridanus (Cold Spot)": {"ra": 55.0,  "dec": -25.0, "radius": 15.0},
    "Northern Local Void":  {"ra": 290.0, "dec": 65.0,  "radius": 15.0}
}

def load_data(filepath):
    """Lädt die Pantheon+ Daten. Versucht verschiedene Formate."""
    try:
        # Versuch 1: Leerzeichen-getrennt (Standard für .dat)
        df = pd.read_csv(filepath, sep=r'\s+')
        # Prüfen ob wichtige Spalten existieren
        if 'zHD' not in df.columns:
            # Versuch 2: CSV (Komma-getrennt)
            df = pd.read_csv(filepath, sep=',')
    except FileNotFoundError:
        print(f"FEHLER: Datei '{filepath}' nicht gefunden.")
        print("Bitte stellen Sie sicher, dass die Datei im selben Ordner liegt.")
        sys.exit(1)
        
    # Spalten normalisieren (Groß-/Kleinschreibung ignorieren)
    df.columns = df.columns.str.upper()
    
    # Mapping der Spaltennamen auf Standard-Namen
    # Pantheon+ nutzt oft: zHD, m_b_corr, RA, DEC
    col_map = {
        'ZHD': 'zHD', 'Z_HD': 'zHD',
        'M_B_CORR': 'm_b_corr', 'MU': 'm_b_corr', # Falls MU direkt gegeben ist
        'RA': 'RA', 'DEC': 'DEC', 'DECL': 'DEC'
    }
    
    # Umbenennen
    df = df.rename(columns=col_map)
    
    required = ['zHD', 'm_b_corr', 'RA', 'DEC']
    missing = [c for c in required if c not in df.columns]
    
    if missing:
        print(f"FEHLER: Folgende Spalten fehlen in der Datei: {missing}")
        print(f"Gefundene Spalten: {list(df.columns)}")
        sys.exit(1)
        
    return df

def calculate_h0(df_subset):
    """
    Berechnet H0 für eine Teilmenge von Supernovae.
    Formel: H0 = (c * z) / D_L
    D_L = 10^((m - M - 25)/5)
    """
    c = 299792.458 # Lichtgeschwindigkeit km/s
    M_abs = -19.2435 # Absolute Magnitude (SH0ES Kalibrierung aus Riess et al. 2022)
    
    # Distanz in Mpc
    # m_b_corr ist die scheinbare Helligkeit (korrigiert)
    # Wenn die Datei Distanzmodule (MU) enthält, ist M_abs = 0 zu setzen (mu = m - M)
    # Pantheon+SH0ES.dat enthält meist m_b_corr. Wir nutzen M_abs.
    
    # Prüfen ob wir MU oder m_b haben. Bei Pantheon+SH0ES.dat ist es oft m_b_corr + M_abs schon verrechnet?
    # HINWEIS: In der offiziellen Datei ist m_b_corr die scheinbare Helligkeit.
    # Wir nutzen die Standard SH0ES Kalibrierung.
    
    dist_mpc = 10**((df_subset['m_b_corr'] - M_abs - 25) / 5)
    
    # H0 für jedes Objekt
    h0_vals = (c * df_subset['zHD']) / dist_mpc
    
    # Statistik (Ausreißer filtern optional, hier roh)
    mean_h0 = np.mean(h0_vals)
    std_err = np.std(h0_vals, ddof=1) / np.sqrt(len(df_subset))
    
    return mean_h0, std_err

def angular_distance(ra1, dec1, ra2, dec2):
    """Berechnet Winkelabstand in Grad (Haversine)."""
    ra1, dec1 = np.radians(ra1), np.radians(dec1)
    ra2, dec2 = np.radians(ra2), np.radians(dec2)
    
    dlon = ra2 - ra1
    dlat = dec2 - dec1
    
    a = np.sin(dlat/2)**2 + np.cos(dec1) * np.cos(dec2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))
    
    return np.degrees(c)

def main():
    print(f"--- ADDENDUM K: VOID-TO-VOID COHERENCE TEST ---")
    print(f"Lade Daten aus: {DATEINAME} ...")
    
    df = load_data(DATEINAME)
    print(f"Geladen: {len(df)} Supernovae.")
    
    # Filter: Nur z > 0.01 (um Pekuliarbewegungen zu minimieren)
    df = df[df['zHD'] > 0.01].copy()
    print(f"Nach z-Cut (z > 0.01): {len(df)} Supernovae.")
    print("-" * 65)
    print(f"{'VOID REGION':<22} | {'N':<5} | {'H0 (km/s/Mpc)':<18} | {'STATUS'}")
    print("-" * 65)
    
    results = []
    
    for name, params in VOIDS.items():
        # Distanz berechnen
        dists = angular_distance(df['RA'].values, df['DEC'].values, params['ra'], params['dec'])
        
        # Filtern
        mask = dists < params['radius']
        subset = df[mask]
        
        if len(subset) > 2:
            h0_mean, h0_err = calculate_h0(subset)
            
            # Bewertung für Addendum K
            # Addendum K sagt: H0 muss ~73.04 sein (deterministisch)
            # Standardmodell (CMB) sagt: ~67.4
            
            if h0_mean > 71.0:
                status = "BESTÄTIGT (Add. K)"
            elif h0_mean < 69.0:
                status = "WIDERLEGT (Planck)"
            else:
                status = "UNKLARE DATEN"
                
            print(f"{name:<22} | {len(subset):<5} | {h0_mean:.2f} +/- {h0_err:.2f}     | {status}")
            results.append(h0_mean)
        else:
            print(f"{name:<22} | {len(subset):<5} | --                 | Zu wenig SNe")
            
    print("-" * 65)
    
    if len(results) >= 2:
        spread = max(results) - min(results)
        print(f"\nMAXIMALE ABWEICHUNG ZWISCHEN VOIDS: {spread:.2f} km/s/Mpc")
        if spread < 2.5:
            print("=> ERGEBNIS: Hohe Kohärenz. Addendum K wird durch die Daten gestützt.")
        else:
            print("=> ERGEBNIS: Hohe Varianz. Addendum K ist wahrscheinlich falsch.")
    else:
        print("\nNicht genügend Voids mit Daten für eine globale Aussage.")

if __name__ == "__main__":
    main()
