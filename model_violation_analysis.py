#!/usr/bin/env python3
"""
Robustness analysis: test QESPM reconstruction quality when model assumptions
are systematically violated.

Three violation scenarios:
  1. Non-perfect conductor (thin dielectric layer on conductor)
  2. Non-uniform background field E_bg (gradient across FOV)
  3. Finite surface slope (Ak -> 1, beyond small-slope limit)

Output: figV5_model_violation.pdf
"""

import numpy as np
from numpy.fft import fft2, ifft2, fftfreq
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt

# Parameters
e = 1.602176634e-19
m = 40 * 1.66053906660e-27
epsilon0 = 8.8541878128e-12
omega_x = 2 * np.pi * 1e6
E_bg_nom = 1e4
delta_omega = 2 * np.pi * 5
h = 40e-6

N = 128
FOV = 128e-6
dx = FOV / N
x = np.linspace(0, FOV, N, endpoint=False)
y = np.linspace(0, FOV, N, endpoint=False)
X, Y = np.meshgrid(x, y)

k_vals = 2 * np.pi * fftfreq(N, dx)
KX, KY = np.meshgrid(k_vals, k_vals)
K = np.sqrt(KX**2 + KY**2)
K[0, 0] = 1e-12

# ============================================================
# Helper: ITF forward + inverse
# ============================================================
def forward(z_s, h_val, Ebg=E_bg_nom):
    phi_s = Ebg * z_s
    phi_fft = fft2(phi_s)
    phi_at_ion = phi_fft * np.exp(-K * h_val)
    dw_fft = (e/(2*m*omega_x)) * (-KX**2) * phi_at_ion
    dw_fft[0, 0] = 0
    return np.real(ifft2(dw_fft))

def invert(dw, h_val, lamb=1e-4):
    dw_fft = fft2(dw)
    H = -(e*E_bg_nom/(2*m*omega_x)) * KX**2 * np.exp(-K * h_val)
    lam = lamb * np.max(np.abs(H))
    H_inv = np.conj(H) / (np.abs(H)**2 + lam**2)
    H_inv[0, 0] = 0
    return np.real(ifft2(H_inv * dw_fft))

# ============================================================
# Generate test surface (sinusoidal grating, variable slope)
# ============================================================
np.random.seed(123)

lambda_s = np.logspace(np.log10(1e-6), np.log10(50e-6), 20)
A_s = np.logspace(np.log10(10e-9), np.log10(500e-9), 20)

# ============================================================
# Test 1: Dielectric layer (non-perfect conductor)
# ============================================================
# Exact linear 3-layer solution (verified analytically, Sec. Methods):
# Phi^(1)(k,h) = E_bg z_s(k) F(k) e^{-kh},  F(k) = eps_r e^{kd}/(eps_r cosh(kd)+sinh(kd)),
# obtained by matching the potential and D-field at the oxide-vacuum interface.
# Thin-film limit: F = 1 + kd(eps_r-1)/eps_r  (a slight ENHANCEMENT: the effective
# image plane is displaced toward the ion by the dielectric).

def dielectric_factor(k, d, eps_r):
    kd = k * d
    return eps_r * np.exp(kd) / (eps_r * np.cosh(kd) + np.sinh(kd))

epsilon_r_list = [2.0, 3.9, 10.0]  # polymer, SiO2, high-k
d_dielectric = 2e-6  # 2 um thick overlayer (k-dependence visible in passband)
h_test = 10e-6
amp_err_vals = []

# 1D band-limited random profile with power inside the ITF passband
rng = np.random.default_rng(42)
kxv = 2 * np.pi * fftfreq(N, dx)
band1d = (np.abs(kxv) >= 2*np.pi/60e-6) & (np.abs(kxv) <= 2*np.pi/25e-6)
prof = np.real(np.fft.ifft(rng.normal(0, 1, N) * band1d))
prof *= 200e-9 / np.std(prof)
z_diel = np.tile(prof, (N, 1))  # x-varying surface

for eps_r in epsilon_r_list:
    f_corr = dielectric_factor(K, d_dielectric, eps_r)
    z_eff = np.real(ifft2(fft2(z_diel) * f_corr))

    dw_dielectric = forward(z_eff, h_test)
    dw_dielectric += np.random.normal(0, delta_omega, (N, N))

    # Invert assuming PERFECT CONDUCTOR (incorrect model)
    z_recon_wrong = invert(dw_dielectric, h_test)
    # Dielectric rescales the reconstructed heights -> amplitude error
    amp_err = abs(1.0 - np.std(z_recon_wrong)/np.std(z_diel))
    amp_err_vals.append(amp_err)
    print(f"  dielectric eps_r={eps_r}: F(k) over passband = "
          f"{dielectric_factor(2*np.pi/60e-6, d_dielectric, eps_r):.4f} .. "
          f"{dielectric_factor(2*np.pi/25e-6, d_dielectric, eps_r):.4f}, "
          f"amplitude error = {amp_err:.3f}")

# ============================================================
# Test 2: Non-uniform E_bg
# ============================================================
gradient_pct = np.array([0, 5.0, 10.0, 20.0, 50.0])
corr_gradient = []
lam_grad = 50e-6
grad_A = 200e-9
z_grad = grad_A * np.sin(2*np.pi*X/lam_grad)

for g in gradient_pct:
    # E_bg varies across FOV: +-g%
    E_bg_spatial = E_bg_nom * (1 + g/100 * (X - FOV/2) / (FOV/2))
    dw_grad = forward(z_grad, h_test, Ebg=E_bg_spatial)
    dw_grad += np.random.normal(0, delta_omega, (N, N))
    z_recon_grad = invert(dw_grad, h_test)
    corr_gradient.append(np.corrcoef(z_grad.ravel(), z_recon_grad.ravel())[0, 1])

# ============================================================
# Test 3: Slope violation (Ak parameter sweep)
# ============================================================
# Fixed wavelength inside the ITF passband; sweep amplitude A = Ak/k so
# that the slope parameter Ak is the independent variable.
lam_slope = 50e-6
k_slope = 2 * np.pi / lam_slope
h_slope = 10e-6  # smaller height makes higher-order harmonics observable
Ak_values = np.logspace(np.log10(0.01), np.log10(1.0), 20)
corr_slope = []

for ak in Ak_values:
    A_val = ak / k_slope
    z_slope = A_val * np.sin(k_slope * X)

    # "True" data to second order in the boundary perturbation
    # Phi^(2)(x,y,0) = -z_s dPhi^(1)/dz|_0  ->  z_eff = A sin(kx)
    # + (k A^2/2)(1 - cos 2kx): DC shift (filtered) + 2k harmonic.
    z_eff = z_slope + (k_slope * A_val**2 / 2.0) * (1.0 - np.cos(2 * k_slope * X))

    dw_slope = forward(z_eff, h_slope)
    dw_slope += np.random.normal(0, delta_omega, (N, N))
    z_recon_slope = invert(dw_slope, h_slope)
    c = np.corrcoef(z_slope.ravel(), z_recon_slope.ravel())[0, 1]
    corr_slope.append(c)

# ============================================================
# Figure V5: Model violation analysis
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Panel (a): Dielectric effect
colors_d = ['#2196F3', '#4CAF50', '#FF9800']
for eps_r, c, col in zip(epsilon_r_list, amp_err_vals, colors_d):
    axes[0].bar(epsilon_r_list.index(eps_r), c, color=col, alpha=0.7,
                label=f'$\\varepsilon_r = {eps_r}$')
axes[0].set_xticks(range(len(epsilon_r_list)))
axes[0].set_xticklabels([f'{r}' for r in epsilon_r_list])
axes[0].set_ylabel('Amplitude error $|1-R_q^{\\rm rec}/R_q|$')
axes[0].set_title(f'(a) Dielectric overlayer ($d={d_dielectric*1e6:.0f}$ µm)\n'
                  f'$h={h_test*1e6:.0f}$ µm, band-limited profile')
axes[0].legend(fontsize=9)
axes[0].grid(True, alpha=0.3, axis='y')
axes[0].set_ylim(0, 0.5)

# Panel (b): E_bg non-uniformity
axes[1].plot(gradient_pct, 1 - np.array(corr_gradient), 'o-', color='teal', lw=2, ms=8)
axes[1].set_xlabel('$E_{\\rm bg}$ gradient [% across FOV]')
axes[1].set_ylabel('Reconstruction error $1-r$')
axes[1].set_title(f'(b) Non-uniform $E_{{\\rm bg}}$\n'
                  f'$h={h_test*1e6:.0f}$ µm, $\\lambda={lam_grad*1e6:.0f}$ µm')
axes[1].grid(True, alpha=0.3)
axes[1].set_ylim(0, 0.05)

# Panel (c): Slope violation
valid = ~np.isnan(corr_slope)
Ak_valid = Ak_values[valid]
corr_valid = np.array(corr_slope)[valid]
axes[2].semilogx(Ak_valid, 1 - corr_valid, 'o-', color='crimson', lw=2, ms=6)
axes[2].axvline(1.0, color='gray', ls='--', lw=1, label='$Ak=1$ (slope limit)')
axes[2].set_xlabel('Surface slope parameter $Ak$')
axes[2].set_ylabel('Reconstruction error $1-r$')
axes[2].set_title(f'(c) Large-slope breakdown\n$h = {h_slope*1e6:.0f}$ µm, '
                  f'$\\lambda = {lam_slope*1e6:.0f}$ µm')
axes[2].legend(fontsize=9, loc='upper left')
axes[2].grid(True, alpha=0.3)
axes[2].set_ylim(0, 0.12)

# Annotate regimes
axes[2].annotate('Valid', xy=(0.06, 0.08), xycoords='axes fraction',
                 fontsize=10, color='green')
axes[2].annotate('Breakdown', xy=(0.72, 0.85), xycoords='axes fraction',
                 fontsize=10, color='red')

plt.tight_layout()
plt.savefig('manuscript/figV5_model_violation.pdf', dpi=150, bbox_inches='tight')
plt.close()
print('Saved: figV5_model_violation.pdf')

# Print summary
print(f"\nModel violation summary:")
print(f"  Dielectric (d=2um) amplitude error: "
      f"{[f'{a:.3f}' for a in amp_err_vals]}")
print(f"  E_bg gradient 10%: 1-r = {1-corr_gradient[2]:.4f}")
print(f"  E_bg gradient 50%: 1-r = {1-corr_gradient[4]:.4f}")
print(f"  Slope Ak=0.1: 1-r = {1-corr_valid[np.argmin(np.abs(Ak_valid-0.1))]:.4f}")
print(f"  Slope Ak=1.0: 1-r = {1-corr_valid[np.argmin(np.abs(Ak_valid-1.0))]:.4f}")
