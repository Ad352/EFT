import pandas as pd
import numpy as np
import sys

# --- KONFIGURATION ---
DATEINAME = "Pantheon+SH0ES.dat"

# Referenzpunkt: CMB Cold Spot (Eridanus)
# Koordinaten ca. RA 55°, Dec -25° (Chilton et al.)
COLD_SPOT = {"ra": 55.0, "dec": -25.0}

# Zu testende Regionen (Voids & Supercluster)
REGIONS = {
    "Eridanus (Cold Spot)": {"ra": 55.0,  "dec": -25.0, "radius": 20.0}, # Nullpunkt
    "Bootes Void":          {"ra": 223.0, "dec": 46.0,  "radius": 15.0},
    "Local Void":           {"ra": 280.0, "dec": 18.0,  "radius": 20.0},
    "Northern Local Void":  {"ra": 290.0, "dec": 65.0,  "radius": 15.0},
    "Southern Supervoid":   {"ra": 15.0,  "dec": -45.0, "radius": 20.0},
    # Gegenrichtung (Antipode zum Cold Spot -> RA ~235, Dec ~+25)
    "Antipode Region":      {"ra": 235.0, "dec": 25.0,  "radius": 20.0}
}

def load_data(filepath):
    try:
        df = pd.read_csv(filepath, sep=r'\s+')
        if 'zHD' not in df.columns:
            df = pd.read_csv(filepath, sep=',')
    except FileNotFoundError:
        print(f"FEHLER: Datei '{filepath}' nicht gefunden.")
        sys.exit(1)
        
    df.columns = df.columns.str.upper()
    col_map = {'ZHD': 'zHD', 'Z_HD': 'zHD', 'M_B_CORR': 'm_b_corr', 'RA': 'RA', 'DEC': 'DEC', 'DECL': 'DEC'}
    df = df.rename(columns=col_map)
    return df

def angular_separation(ra1, dec1, ra2, dec2):
    """Berechnet Winkelabstand in Grad zwischen zwei Punkten"""
    ra1, dec1 = np.radians(ra1), np.radians(dec1)
    ra2, dec2 = np.radians(ra2), np.radians(dec2)
    dlon, dlat = ra2 - ra1, dec2 - dec1
    a = np.sin(dlat/2)**2 + np.cos(dec1) * np.cos(dec2) * np.sin(dlon/2)**2
    return np.degrees(2 * np.arcsin(np.sqrt(np.clip(a, 0, 1))))

def main():
    print(f"--- EFT-LAPSE ANISOTROPY TEST (COLD SPOT DIPOLE) ---")
    df = load_data(DATEINAME)
    # Filter: Nur z > 0.015 (sehr lokale SNe ausschließen)
    df = df[df['zHD'] > 0.015].copy()
    
    print(f"{'REGION':<22} | {'ANGLE':<5} | {'N':<5} | {'H0 (km/s/Mpc)':<18} | {'TREND'}")
    print("-" * 75)
    
    results = []
    
    for name, params in REGIONS.items():
        # 1. Winkelabstand zum Cold Spot berechnen
        angle = angular_separation(COLD_SPOT['ra'], COLD_SPOT['dec'], params['ra'], params['dec'])
        
        # 2. SNe in dieser Region finden
        dists = angular_separation(df['RA'].values, df['DEC'].values, params['ra'], params['dec'])
        subset = df[dists < params['radius']]
        
        if len(subset) > 5:
            # 3. H0 berechnen
            M_abs = -19.2435
            dist_mpc = 10**((subset['m_b_corr'] - M_abs - 25) / 5)
            h0_vals = (299792.458 * subset['zHD']) / dist_mpc
            
            mean_h0 = np.mean(h0_vals)
            std_err = np.std(h0_vals, ddof=1) / np.sqrt(len(subset))
            
            # Trend bewerten
            if mean_h0 < 65: trend = "LOW (Cold)"
            elif mean_h0 > 71: trend = "HIGH (Warm)"
            else: trend = "Medium"
            
            print(f"{name:<22} | {angle:>3.0f}°  | {len(subset):<5} | {mean_h0:.2f} +/- {std_err:.2f}     | {trend}")
            results.append({"angle": angle, "h0": mean_h0})
        else:
            print(f"{name:<22} | {angle:>3.0f}°  | {len(subset):<5} | --                 | Zu wenig Daten")

    print("-" * 75)
    
    # Korrelation prüfen
    if len(results) >= 3:
        angles = [r['angle'] for r in results]
        h0s = [r['h0'] for r in results]
        corr = np.corrcoef(angles, h0s)[0, 1]
        
        print(f"\nKorrelation (Winkel vs. H0): {corr:.2f}")
        if corr > 0.5:
            print("=> POSITIVE KORRELATION: Je weiter vom Cold Spot, desto höher H0.")
            print("   Dies stützt die Hypothese eines kosmischen Dipols (EFT-Lapse Addendum J).")
        elif corr < -0.5:
            print("=> NEGATIVE KORRELATION: H0 sinkt mit Abstand zum Cold Spot.")
        else:
            print("=> KEINE KLARE KORRELATION: Die Verteilung scheint zufällig.")

if __name__ == "__main__":
    main()
