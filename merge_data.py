import pandas as pd

# Bestehende Phi-Daten und neue Helligkeits-Daten laden
df_phi = pd.read_csv('sne_mit_phi.csv')
df_mags = pd.read_csv('my_sne_with_mags.csv')

# Zusammenführen über den Namen
# (Wir nehmen nur die Spalten 'name', 'mag' und 'mag_err' aus der neuen Datei)
df_final = pd.merge(df_phi, df_mags[['name', 'mag', 'mag_err']], on='name', how='left')

# Speichern
df_final.to_csv('sne_final_analysis.csv', index=False)
print("Finale Datei 'sne_final_analysis.csv' gespeichert.")
