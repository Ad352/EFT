import numpy as np
from scipy.linalg import expm
from scipy.sparse import diags, kron as spkron, csr_matrix, eye as sp_eye
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt
import pandas as pd

# --- Parameters ---
L = 4          # Number of sites
dloc = 4       # Local dimension per site (Casimir levels)
N_full = dloc**L # 256
trunc = 128    # Truncation size

g = 0.5        # Coupling
d_angle = 30 * np.pi / 180
cos_d, sin_d = np.cos(d_angle), np.sin(d_angle)

# --- Build Full Sparse Hamiltonian ---
def site_op(op, i):
    left = sp_eye(dloc**i, format='csr') if i > 0 else sp_eye(1, format='csr')
    right = sp_eye(dloc**(L - i - 1), format='csr') if i < L-1 else sp_eye(1, format='csr')
    if i == 0:
        return spkron(op, right, format='csr')
    elif i == L-1:
        return spkron(left, op, format='csr')
    else:
        return spkron(spkron(left, op, format='csr'), right, format='csr')

def bond_op(op_bond, i):
    left = sp_eye(dloc**i, format='csr') if i > 0 else sp_eye(1, format='csr')
    right = sp_eye(dloc**(L - i - 2), format='csr') if i < L-2 else sp_eye(1, format='csr')
    if L == 2:
        return op_bond
    if i == 0:
        return spkron(op_bond, right, format='csr')
    elif i == L-2:
        return spkron(left, op_bond, format='csr')
    else:
        return spkron(spkron(left, op_bond, format='csr'), right, format='csr')

# Local operators
e_loc = np.linspace(1.0, 2.0, dloc)
H_iso_loc = diags([e_loc], [0], shape=(dloc, dloc), format='csr')

# Anisotropic local (Dipole bias)
H_aniso_loc = g * sin_d * diags([np.ones(dloc)], [0], shape=(dloc, dloc), format='csr')

# Hopping (Dipole gradient off-diagonal)
H_hop_loc = g * cos_d * diags([np.ones(dloc-1)], [-1], shape=(dloc, dloc), format='csr')
H_hop_loc += H_hop_loc.T
O_hop = spkron(H_hop_loc, sp_eye(dloc, format='csr')) + spkron(sp_eye(dloc, format='csr'), H_hop_loc)

H_iso = csr_matrix((N_full, N_full))
for i in range(L):
    H_iso += site_op(H_iso_loc, i)

H_aniso = csr_matrix((N_full, N_full))
for i in range(L):
    H_aniso += 0.1 * site_op(H_aniso_loc, i)
for i in range(L-1):
    H_aniso += 0.3 * bond_op(O_hop, i) 

H_full = H_iso + H_aniso

# --- Truncation ---
evals, evecs = eigsh(H_full, k=trunc, which='SM')
H_trunc = np.real(evecs.T @ H_full @ evecs)
H_iso_trunc = np.real(evecs.T @ H_iso @ evecs)
V = np.real(evecs.T @ H_aniso @ evecs)

pre_off = np.linalg.norm(H_trunc - np.diag(np.diag(H_trunc)))

# --- Magnus-SW Transformation (Order 12) ---
def ad(A, B):
    return A @ B - B @ A

def magnus_S(V, H_iso, order=12):
    S = np.zeros_like(V, dtype=complex)
    adV_cum = V.copy().astype(complex)
    for k in range(order):
        coeff = ((-1j / 2) ** (k+1)) / np.math.factorial(k+1)
        S += coeff * adV_cum
        adV_cum = ad(V, adV_cum @ H_iso)
    return S

S = magnus_S(V, H_iso_trunc, order=12)
U = expm(S)
H_eff = U @ H_trunc @ U.T.conj()

post_off = np.linalg.norm(H_eff - np.diag(np.diag(H_eff)))
bd_percent = 100 * (1 - post_off / np.linalg.norm(H_eff))

e_eff_sort = np.sort(np.real(np.diag(H_eff)))
deg_count = np.sum(np.diff(e_eff_sort) < 1e-3)

# --- Lindblad Dissipation ---
gamma = 1 / (4 * np.pi)  # KSS bound ~0.0796
dt, steps = 0.05, 100
psi0 = np.zeros(trunc, dtype=complex); psi0[0] = 1.0
rho = np.outer(psi0, psi0.conj())

np.random.seed(42)
Ls = [np.diag(np.random.uniform(0.05, 0.2, trunc)) for _ in range(L)]

purities = []
times = []

for t in range(steps):
    dH = -1j * (H_eff @ rho - rho @ H_eff)
    dL = np.zeros_like(rho, dtype=complex)
    for Lj in Ls:
        dL += Lj @ rho @ Lj.conj().T - 0.5 * (Lj.conj().T @ Lj @ rho + rho @ Lj.conj().T @ Lj)
    rho += dt * (dH + gamma * dL)
    purities.append(np.real(np.trace(rho @ rho)))
    times.append(t * dt)

final_purity = purities[-1]

# --- Output & Plots ---
plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.plot(np.sort(np.real(np.diag(H_iso_trunc))), 'b--', alpha=0.7, label='Isotrope Basis')
plt.plot(np.sort(np.real(np.diag(H_trunc))), 'orange', alpha=0.7, label='Full (Pre-SW)')
plt.plot(e_eff_sort, 'r-', linewidth=1.5, label='Effektiv (Post-SW Order=12)')
plt.xlabel('Krylov Mode Index')
plt.ylabel('Energie Eigenwert')
plt.title('Spektrales Splitting (Mode E Emergenz)')
plt.legend()
plt.grid(True, alpha=0.3)

plt.subplot(1, 2, 2)
plt.plot(times, purities, 'g-', linewidth=2, label=r'Purity $Tr(\rho^2)$')
plt.axhline(1.0, color='k', linestyle='--', alpha=0.5)
plt.ylim(0.8, 1.05)
plt.xlabel('Zeit (a.u.)')
plt.ylabel('Purity')
plt.title(r'Lindblad Dissipation ($\eta/s=1/4\pi$)')
plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('eft_lapse_v52_validation.png', dpi=300)

metrics_df = pd.DataFrame({
    'N_full': [N_full],
    'Truncation': [trunc],
    'Magnus_Order': [12],
    'Pre_OffDiag_Norm': [pre_off],
    'Post_OffDiag_Norm': [post_off],
    'Block_Diagonal_Pct': [bd_percent],
    'Pd0_Sectors': [deg_count],
    'Final_Purity': [final_purity]
})
metrics_df.to_csv('eft_lapse_v52_metrics.csv', index=False)

print(f"Metrics: Off-Diag reduced {pre_off:.2e} -> {post_off:.2e}")
print(f"Pd=0 Sectors: {deg_count}")
print(f"Final Purity: {final_purity:.4f}")