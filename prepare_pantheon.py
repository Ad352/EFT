import pandas as pd

# 1. Pantheon Datei laden
# Stellen Sie sicher, dass die Datei im selben Ordner liegt!
filename = 'Pantheon+SH0ES.dat' # Passen Sie den Namen ggf. an

print(f"Lese {filename}...")
# sep='\s+' bedeutet: Trennung durch beliebig viele Leerzeichen
df = pd.read_csv(filename, sep='\s+', comment='#')

# 2. Relevante Spalten auswählen und umbenennen
mapping = {
    'CID': 'name',
    'zCMB': 'z',
    'RA': 'ra',
    'DEC': 'dec'
}

# Nur die Spalten nehmen, die existieren
df_final = df[list(mapping.keys())].rename(columns=mapping)

# 3. Speichern
df_final.to_csv('my_sne.csv', index=False)

print("-" * 30)
print(f"Erfolg! {len(df_final)} Supernovae gefunden.")
print("Datei 'my_sne.csv' wurde erstellt.")
print("Nächster Schritt: python interpolate_phi.py")
