import pandas as pd

# 1. Pantheon Datei laden
filename = 'Pantheon+SH0ES.dat' 

print(f"Lese {filename}...")
df = pd.read_csv(filename, sep='\s+', comment='#')

# 2. Relevante Spalten auswählen
# m_b_corr ist die für Bias und Form korrigierte Helligkeit
mapping = {
    'CID': 'name',
    'zCMB': 'z',
    'RA': 'ra',
    'DEC': 'dec',
    'm_b_corr': 'mag',
    'm_b_corr_err_DIAG': 'mag_err'
}

df_final = df[list(mapping.keys())].rename(columns=mapping)

# 3. Speichern
df_final.to_csv('my_sne_with_mags.csv', index=False)

print("-" * 30)
print(f"Erfolg! {len(df_final)} SNe mit Helligkeiten extrahiert.")
print("Datei 'my_sne_with_mags.csv' erstellt.")
