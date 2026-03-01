import pandas as pd
import numpy as np
import sys

# --- KONFIGURATION ---
DATEINAME = "Pantheon+SH0ES.dat"

# Wir definieren Bootes zweimal: Einmal "Full", einmal "Core"
VOIDS = {
    "Bootes (Full 15deg)": {"ra": 223.0, "dec": 46.0, "radius": 15.0},
    "Bootes (Core 10deg)": {"ra": 223.0, "dec": 46.0, "radius": 10.0},
    "Local Void":          {"ra": 280.0, "dec": 18.0, "radius": 20.0},
}

def load_data(filepath):
    try:
        # Versuch 1: Whitespace-separated
        df = pd.read_csv(filepath, sep=r'\s+')
        if 'zHD' not in df.columns:
             # Versuch 2: Comma-separated
            df = pd.read_csv(filepath, sep=',')
    except FileNotFoundError:
        print(f"FEHLER: Datei '{filepath}' nicht gefunden.")
        sys.exit(1)
        
    df.columns = df.columns.str.upper()
    col_map = {'ZHD': 'zHD', 'Z_HD': 'zHD', 'M_B_CORR': 'm_b_corr', 'RA': 'RA', 'DEC': 'DEC', 'DECL': 'DEC'}
    df = df.rename(columns=col_map)
    return df

def angular_distance(ra1, dec1, ra2, dec2):
    ra1, dec1 = np.radians(ra1), np.radians(dec1)
    ra2, dec2 = np.radians(ra2), np.radians(dec2)
    dlon, dlat = ra2 - ra1, dec2 - dec1
    a = np.sin(dlat/2)**2 + np.cos(dec1) * np.cos(dec2) * np.sin(dlon/2)**2
    return np.degrees(2 * np.arcsin(np.sqrt(np.clip(a, 0, 1))))

def main():
    print(f"--- BOOTES VOID ANOMALY TEST ---")
    df = load_data(DATEINAME)
    
    # Filter: Nur z > 0.01 (Pekuliarbewegungen minimieren)
    df = df[df['zHD'] > 0.01].copy()
    
    print(f"{'REGION':<20} | {'N':<5} | {'H0 (km/s/Mpc)':<18} | {'INTERPRETATION'}")
    print("-" * 65)
    
    results = {}
    
    for name, params in VOIDS.items():
        dists = angular_distance(df['RA'].values, df['DEC'].values, params['ra'], params['dec'])
        subset = df[dists < params['radius']]
        
        if len(subset) > 2:
            # H0 Berechnung (SH0ES Kalibrierung M = -19.2435)
            M_abs = -19.2435
            dist_mpc = 10**((subset['m_b_corr'] - M_abs - 25) / 5)
            h0_vals = (299792.458 * subset['zHD']) / dist_mpc
            
            mean_h0 = np.mean(h0_vals)
            std_err = np.std(h0_vals, ddof=1) / np.sqrt(len(subset))
            
            # Interpretation
            if mean_h0 < 62:
                interp = "ANOMALIE (Sehr tief)"
            elif mean_h0 > 70:
                interp = "SH0ES / Addendum K"
            else:
                interp = "Planck-Level"
                
            print(f"{name:<20} | {len(subset):<5} | {mean_h0:.2f} +/- {std_err:.2f}     | {interp}")
            results[name] = mean_h0
        else:
            print(f"{name:<20} | {len(subset):<5} | --                 | Zu wenig Daten")

    print("-" * 65)
    
    # Analyse der Veränderung
    if "Bootes (Full 15deg)" in results and "Bootes (Core 10deg)" in results:
        diff = results["Bootes (Core 10deg)"] - results["Bootes (Full 15deg)"]
        print(f"\nVeränderung Full -> Core: {diff:+.2f} km/s/Mpc")
        
        if abs(diff) < 1.5:
             print("=> ERGEBNIS: Die Anomalie ist REAL und stabil (kein Randeffekt).")
             print("   Dies widerlegt Addendum K (H0 sollte ~73 sein).")
        elif diff > 3.0:
             print("=> ERGEBNIS: Der Wert steigt im Kern an. Der niedrige Wert war ein Randeffekt.")

if __name__ == "__main__":
    main()
