#!/usr/bin/env python3
"""
Monte Carlo uncertainty propagation (JCGM 101:2008 Supplement 1).
Propagates input uncertainties through the full forward + inversion chain.

Output: figV6_monte_carlo.pdf
"""

import numpy as np
from numpy.fft import fft2, ifft2, fftfreq
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt

np.random.seed(789)

# Parameters
e = 1.602176634e-19
m = 40 * 1.66053906660e-27
epsilon0 = 8.8541878128e-12
omega_x_nom = 2 * np.pi * 1e6
E_bg_nom = 1e4
h_nom = 40e-6
delta_omega_nom = 2 * np.pi * 5

N = 128
FOV = 128e-6
dx = FOV / N
x = np.linspace(0, FOV, N, endpoint=False)
X, Y = np.meshgrid(x, x)

k_vals = 2 * np.pi * fftfreq(N, dx)
KX, KY = np.meshgrid(k_vals, k_vals)
K = np.sqrt(KX**2 + KY**2)
K[0, 0] = 1e-12

# Test surface: sinusoidal grating at the GUM evaluation wavenumber k = 2e5
lambda_s = 2*np.pi/2e5  # 31.4 um, inside the ITF passband at h = 40 um
A_true = 100e-9
z_true = A_true * np.sin(2*np.pi*X/lambda_s)

# Forward
def forward(z_s, h_val, omega_val, Ebg):
    phi_s = Ebg * z_s
    phi_fft = fft2(phi_s)
    phi_at_ion = phi_fft * np.exp(-K * h_val)
    dw_fft = (e/(2*m*omega_val)) * (-KX**2) * phi_at_ion
    dw_fft[0, 0] = 0
    return np.real(ifft2(dw_fft))

def invert(dw, h_val, omega_val, Ebg, lamb=1e-4):
    dw_fft = fft2(dw)
    H = -(e*Ebg/(2*m*omega_val)) * KX**2 * np.exp(-K * h_val)
    lam = lamb * np.max(np.abs(H))
    H_inv = np.conj(H) / (np.abs(H)**2 + lam**2)
    H_inv[0, 0] = 0
    return np.real(ifft2(H_inv * dw_fft))

# ============================================================
# Monte Carlo propagation
# ============================================================
N_MC = 5000
print(f"Running {N_MC} Monte Carlo iterations...")

# Input uncertainty distributions (from GUM budget)
u_h = 100e-9           # 100 nm height uncertainty
u_Ebg = 10.0           # 10 V/m
u_omega = 2*np.pi*10   # 10 Hz secular frequency
u_sigma = 1.13e-8       # surface charge (m equivalent)

# Storage
z_mc = np.zeros((N_MC, N))

for i in range(N_MC):
    if i % 1000 == 0:
        print(f"  Iteration {i}/{N_MC}...")
    
    # Draw random inputs
    h_i = np.random.normal(h_nom, u_h)
    Ebg_i = np.random.normal(E_bg_nom, u_Ebg)
    omega_i = np.random.normal(omega_x_nom, u_omega)
    
    # Add charge contribution (coherent at the grating wavenumber k,
    # matching the single-k GUM evaluation)
    sigma_i = np.random.normal(0, u_sigma)   # scalar: charge density at k
    z_eff = z_true + sigma_i / (2 * epsilon0 * (2*np.pi/lambda_s) * Ebg_i) \
        * np.sin(2*np.pi*X/lambda_s)
    
    # Forward
    dw_i = forward(z_eff, h_i, omega_i, Ebg_i)
    # Add measurement noise
    dw_i += np.random.normal(0, delta_omega_nom, (N, N))
    
    # Invert with NOMINAL parameters (as done in practice)
    z_recon_i = invert(dw_i, h_nom, omega_x_nom, E_bg_nom)
    
    # Store cross-section at centre
    z_mc[i, :] = z_recon_i[N//2, :]

print(f"Done. Computing statistics...")

# Statistics
z_mean = np.mean(z_mc, axis=0)
z_std = np.std(z_mc, axis=0)
z_lo = np.percentile(z_mc, 2.5, axis=0)
z_hi = np.percentile(z_mc, 97.5, axis=0)

# True cross-section
z_true_cross = z_true[N//2, :]

# ============================================================
# Figure V6: Monte Carlo uncertainty
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

x_um = x * 1e6

# Panel (a): Cross-section with 95% coverage interval
axes[0].fill_between(x_um, z_lo*1e9, z_hi*1e9, alpha=0.3, color='blue',
                      label='95% coverage')
axes[0].plot(x_um, z_true_cross*1e9, 'k-', lw=1.5, label='True $z$')
axes[0].plot(x_um, z_mean*1e9, 'r--', lw=1, label='MC mean')
axes[0].set_xlabel('$x$ [µm]')
axes[0].set_ylabel('$z$ [nm]')
axes[0].set_title(f'(a) Reconstruction ±95% CI  ($N_{{\\rm MC}}$={N_MC})')
axes[0].legend(fontsize=8)
axes[0].grid(True, alpha=0.3)

# Panel (b): Uncertainty budget decomposition (spatial)
axes[1].plot(x_um, z_std*1e9, 'b-', lw=1.5, label='Combined $u_c(z)$')
axes[1].set_xlabel('$x$ [µm]')
axes[1].set_ylabel('Standard uncertainty [nm]')
axes[1].set_title('(b) Spatial uncertainty map')
axes[1].grid(True, alpha=0.3)
axes[1].legend(fontsize=9)

# Panel (c): Histogram at the grating crest (x = lambda_s/4)
idx_peak = int(round(lambda_s / dx / 4))  # crest where |sin| ~ 1
z_at_peak = z_mc[:, idx_peak] * 1e9
axes[2].hist(z_at_peak, bins=60, density=True, color='steelblue', alpha=0.7, edgecolor='white')
axes[2].axvline(z_true_cross[idx_peak]*1e9, color='k', ls='--', lw=2, label='True')
axes[2].axvline(np.mean(z_at_peak), color='red', ls='-', lw=1.5, label=f'MC mean ({np.mean(z_at_peak):.1f} nm)')
axes[2].set_xlabel('$z$ [nm]')
axes[2].set_ylabel('Probability density')
axes[2].set_title(f'(c) Distribution at $x$={x_um[idx_peak]:.1f} µm')
axes[2].legend(fontsize=8)

# Annotate stats
axes[2].annotate(f'$u_c$ = {np.std(z_at_peak):.1f} nm\n95% CI: [{np.percentile(z_at_peak,2.5):.1f}, {np.percentile(z_at_peak,97.5):.1f}] nm',
                 xy=(0.95, 0.95), xycoords='axes fraction',
                 ha='right', va='top', fontsize=9,
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
plt.savefig('manuscript/figV6_monte_carlo.pdf', dpi=150, bbox_inches='tight')
plt.close()
print('Saved: figV6_monte_carlo.pdf')

# Print results
print(f"\nMonte Carlo results (N_MC = {N_MC}):")
print(f"  Mean u_c(z) across FOV: {np.mean(z_std)*1e9:.2f} nm")
print(f"  Max u_c(z): {np.max(z_std)*1e9:.2f} nm")
print(f"  At peak: z = {np.mean(z_at_peak):.1f} ± {np.std(z_at_peak):.1f} nm (1σ)")
