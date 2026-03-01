import numpy as np

delta = np.load('twompp_density.npy')
print(f"Array-Form: {delta.shape}")
print(f"Datentyp: {delta.dtype}")
print(f"Min/Max Dichte: {delta.min():.4f} / {delta.max():.4f}")
print(f"Mittlere Dichte: {delta.mean():.4f}")