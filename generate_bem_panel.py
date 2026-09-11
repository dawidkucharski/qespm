#!/usr/bin/env python3
"""Compact main-text panel: independent solver cross-checks of the linear ITF.

Relative error of the linear ITF (first-order boundary perturbation) and of
the independent boundary-element solution against the exact Rayleigh-Floquet
solution, as functions of the perturbation parameter kA.

Prints a verification table of the numbers quoted in main.tex:
  - ITF relative error at A = 100 nm and A = 2 um (expected ~0.14 (kA)^2)
  - BEM-vs-Rayleigh relative difference (expected <= 0.6%)

Saves manuscript/figB1_bem_panel.pdf
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.linalg import solve

from audit_fixes import rayleigh_map, E_BG, E_CHARGE, M_CA40, OMEGA_X

LAM = 40e-6
k_s = 2 * np.pi / LAM
h = 40e-6
C_ITF = E_CHARGE * E_BG / (2 * M_CA40 * OMEGA_X)   # = 1.92e3 m/s


def itf_shift(A):
    """Linear ITF secular shift at the crest, rad/s."""
    return -C_ITF * k_s ** 2 * A * np.exp(-k_s * h)


# ---------------- exact Rayleigh vs ITF (error vs kA) ----------------
A_itf = np.logspace(np.log10(20e-9), np.log10(4e-6), 8)
err_itf, kA_list = [], []
for A in A_itf:
    dw_ex = rayleigh_map(A, k_s, h, N=28)
    dw_lin = itf_shift(A)
    err_itf.append(abs(dw_lin - dw_ex) / abs(dw_ex))
    kA_list.append(k_s * A)
err_itf = np.array(err_itf)
kA_list = np.array(kA_list)

# coefficient of the quadratic law at the largest amplitudes
kA_hi = kA_list[kA_list > 0.2]
coef = np.mean(err_itf[kA_list > 0.2] / kA_hi ** 2)

# ---------------- independent BEM vs exact Rayleigh ----------------
A_bem = np.array([100e-9, 500e-9, 1e-6, 2e-6])
err_bem, kA_bem = [], []
for A in A_bem:
    Nb = 600
    L = 2 * LAM
    s = np.linspace(0, L, Nb, endpoint=False)
    ds = s[1] - s[0]
    xs = s
    zs = A * np.sin(k_s * s)
    r = np.stack([xs, zs], axis=1)
    G = np.zeros((Nb, Nb))
    for i in range(Nb):
        dxv = (r[:, 0] - r[i, 0]) * 2 * np.pi / L
        dzv = (r[:, 1] - r[i, 1]) * 2 * np.pi / L
        arg = 2 * np.cosh(dzv) - 2 * np.cos(dxv)
        G[i, :] = -np.log(arg) * ds / (4 * np.pi)
        G[i, i] = -ds / (2 * np.pi) * (np.log(2 * np.pi * ds / L) - 1.0)
    sigma = np.linalg.solve(G, E_BG * zs)
    x0, z0 = LAM / 4, h
    dxv = (r[:, 0] - x0) * 2 * np.pi / L
    dzv = (r[:, 1] - z0) * 2 * np.pi / L
    csh = np.cosh(dzv)
    cs = np.cos(dxv)
    D = csh - cs
    kp = 2 * np.pi / L
    d2G = -kp ** 2 / (4 * np.pi) * (csh * cs - 1.0) / D ** 2
    d2Phi = np.sum(sigma * d2G) * ds
    dw_b = (E_CHARGE / (2 * M_CA40 * OMEGA_X)) * d2Phi
    dw_ex = rayleigh_map(A, k_s, h, N=28)
    err_bem.append(abs(dw_b - dw_ex) / abs(dw_ex))
    kA_bem.append(k_s * A)
err_bem = np.array(err_bem)
kA_bem = np.array(kA_bem)

# ---------------- verification table ----------------
print("ITF vs exact Rayleigh:")
for A, kA, e in zip(A_itf, kA_list, err_itf):
    print(f"  A = {A*1e9:7.1f} nm   kA = {kA:.4f}   rel err = {e*100:.4f}%")
print(f"  quadratic coefficient (kA>0.2): {coef:.4f}")
print(f"  predicted at A=2 um (kA={k_s*2e-6:.4f}): {coef*(k_s*2e-6)**2*100:.2f}%")
print("BEM vs exact Rayleigh:")
for A, kA, e in zip(A_bem, kA_bem, err_bem):
    print(f"  A = {A*1e9:7.1f} nm   kA = {kA:.4f}   rel diff = {e*100:.3f}%")
print(f"  max BEM diff: {err_bem.max()*100:.2f}%")

# ---------------- figure ----------------
fig, ax = plt.subplots(figsize=(8.2, 4.6))
ax.semilogy(kA_list, err_itf * 100, 'o-', color='#C62828', ms=7, lw=2,
            label='Linear ITF vs exact Rayleigh')
ax.semilogy(kA_list, coef * kA_list ** 2 * 100, '--', color='black', lw=1.6,
            label=r'$0.14\,(kA)^2$ law')
ax.semilogy(kA_bem, err_bem * 100, 's', color='#1f4e79', ms=9, mec='white',
            label='BEM vs exact Rayleigh')
ax.axhline(0.6, color='#1f4e79', ls=':', lw=1.2)
ax.text(0.42, 0.62, '0.6%', color='#1f4e79', fontsize=9, va='bottom')
ax.set_xlabel(r'Perturbation parameter $kA$', fontsize=12)
ax.set_ylabel(r'Relative error [%]', fontsize=12)
ax.set_xlim(0.0, 0.66)
ax.set_ylim(1e-4, 30)
ax.grid(True, alpha=0.3, which='both')
ax.legend(fontsize=9.5, loc='lower right')
ax.set_title(r'Independent solver cross-checks of the linear ITF '
             r'($\lambda=40\,\mu$m, $h=40\,\mu$m)', fontsize=11.5)

plt.tight_layout()
plt.savefig('manuscript/figB1_bem_panel.pdf', dpi=150, bbox_inches='tight')
plt.close()
print('saved manuscript/figB1_bem_panel.pdf')
