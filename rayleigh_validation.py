#!/usr/bin/env python3
"""
Exact (Rayleigh/Floquet) validation of the linear ITF.

Solves the Laplace equation exactly for a grounded conducting surface with
the sinusoidal corrugation z_s(x) = A sin(kx) in a uniform background field
E_bg z-hat, using the Rayleigh mode expansion:

    Phi(x, z) = -E_bg z + sum_n c_n sin(n k x) exp(-n k z)    (z >= max z_s)

The coefficients c_n are determined by the exact Dirichlet boundary
condition Phi(x, A sin(kx)) = 0, giving the linear system

    sum_n c_n M_{mn} = E_bg A delta_{m1},
    M_{mn} = (1/pi) integral_0^{2pi} sin(mu) sin(nu) exp(-n k A sin u) du.

The exact secular frequency shift at the crest (x = lambda/4) is
Delta omega_x = (e/2m omega_x) sum_{n odd} -(nk)^2 c_n exp(-n k h),
which is compared with the first-order ITF prediction.

Output: figV7_rayleigh_validation.pdf
"""

import numpy as np
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt

# ============================================================
# Physical parameters
# ============================================================
E_CHARGE = 1.602176634e-19
M_CA40 = 40 * 1.66053906660e-27
E_BG = 1e4
OMEGA_X = 2 * np.pi * 1e6

lambda_s = 40e-6          # surface period
k_s = 2 * np.pi / lambda_s
h_nom = 40e-6             # ion height


def rayleigh_coeffs(A, N=60):
    """Solve the Rayleigh system for the exact corrugated-boundary solution."""
    u = np.linspace(0, 2 * np.pi, 8000, endpoint=False)
    M = np.zeros((N, N))
    rhs = np.zeros(N)
    for m in range(N):
        sin_m = np.sin((m + 1) * u)
        rhs[m] = E_BG * A * np.mean(sin_m * np.sin(u)) * 2.0
        for n in range(N):
            f = np.sin((n + 1) * u) * sin_m * np.exp(-(n + 1) * k_s * A * np.sin(u))
            M[m, n] = np.mean(f) * 2.0
    return np.linalg.solve(M, rhs)


def exact_dw(A, N=60):
    """Exact secular frequency shift at the crest (x = lambda/4)."""
    c = rayleigh_coeffs(A, N)
    dw = 0.0
    for n in range(0, N, 2):  # odd harmonics survive at x = lambda/4
        nn = n + 1
        dw += -(E_CHARGE / (2 * M_CA40 * OMEGA_X)) * c[n] * (nn * k_s)**2 \
            * np.exp(-nn * k_s * h_nom)
    return dw


def itf_dw(A):
    """First-order ITF prediction at the crest."""
    return -(E_CHARGE * E_BG / (2 * M_CA40 * OMEGA_X)) * A * k_s**2 \
        * np.exp(-k_s * h_nom)


# ============================================================
# Amplitude sweep
# ============================================================
print("Running Rayleigh exact validation for sinusoidal surfaces...")
A_list = np.logspace(np.log10(10e-9), np.log10(2e-6), 15)
dw_exact = []
dw_itf = []

for A in A_list:
    dw_exact.append(exact_dw(A))
    dw_itf.append(itf_dw(A))

dw_exact = np.array(dw_exact)
dw_itf = np.array(dw_itf)
rel_err = np.abs(dw_exact - dw_itf) / np.abs(dw_exact)

# fit leading quadratic coefficient: rel_err ~ alpha * (kA)^2
alpha_fit = np.mean(rel_err / (k_s * A_list)**2)

# ============================================================
# Figure V7: exact validation
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Panel (a): exact vs ITF frequency shift vs amplitude
axes[0].loglog(A_list * 1e9, np.abs(dw_exact) / (2 * np.pi), 'bo-', ms=6, lw=2,
               label='Exact (Rayleigh)')
axes[0].loglog(A_list * 1e9, np.abs(dw_itf) / (2 * np.pi), 'r--', lw=2,
               label='ITF (linear)')
axes[0].set_xlabel('Surface amplitude $A$ [nm]')
axes[0].set_ylabel(r'$|\Delta\omega_x|/2\pi$ [Hz]')
axes[0].set_title(f'(a) Frequency shift vs amplitude '
                  f'($h$={h_nom*1e6:.0f} µm, $\lambda$={lambda_s*1e6:.0f} µm)')
axes[0].legend(fontsize=9)
axes[0].grid(True, alpha=0.3)

# Panel (b): relative error exact vs ITF
axes[1].loglog(A_list * 1e9, rel_err * 100, 'ko-', ms=6, lw=2,
               label='Exact vs ITF')
axes[1].loglog(A_list * 1e9, alpha_fit * (k_s * A_list)**2 * 100, 'g--', lw=1.5,
               label=rf'$\propto (kA)^2$')
axes[1].set_xlabel('Surface amplitude $A$ [nm]')
axes[1].set_ylabel('Relative error [%]')
axes[1].set_title('(b) ITF accuracy vs exact solution')
axes[1].legend(fontsize=9)
axes[1].grid(True, alpha=0.3)

# Panel (c): A/h ratio vs error
z_over_h = A_list / h_nom
axes[2].loglog(z_over_h, rel_err * 100, 'ko-', ms=6, lw=2)
axes[2].axvline(0.1, color='gray', ls='--', lw=1, label='$z_s/h = 0.1$')
axes[2].set_xlabel('$A/h$')
axes[2].set_ylabel('Relative error [%]')
axes[2].set_title('(c) Error vs $z_s/h$ ratio')
axes[2].legend(fontsize=9)
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('manuscript/figV7_rayleigh_validation.pdf', dpi=150,
            bbox_inches='tight')
plt.close()
print('Saved: figV7_rayleigh_validation.pdf')

print(f"\nRayleigh validation summary:")
for A_i in [10e-9, 100e-9, 300e-9, 1e-6, 2e-6]:
    idx = np.argmin(np.abs(A_list - A_i))
    print(f"  A = {A_i*1e9:6.0f} nm: |dw|/2pi (exact) = "
          f"{np.abs(dw_exact[idx])/(2*np.pi):10.1f} Hz, "
          f"relative error = {rel_err[idx]*100:6.2f}%")
print(f"  Quadratic coefficient: rel_err = {alpha_fit:.3f} x (kA)^2")
