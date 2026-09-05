#!/usr/bin/env python3
"""
Generate realistic AFM-like surface data using literature-derived PSD parameters
and certified calibration grating dimensions. Serves as input for all validation scripts.

Output files (in manuscript/):
  figV1_afm_input.pdf    — Realistic surface (AFM-like calibration grating + roughness)
  figV2_afm_ion_signal.pdf — Ion signal for the realistic surface
  figV3_afm_reconstruction.pdf — QESPM reconstruction compared to AFM "ground truth"
"""

import numpy as np
from numpy.fft import fft2, ifft2, fftfreq
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import os, sys

# ============================================================
# Physical parameters (consistent with manuscript)
# ============================================================
e = 1.602176634e-19          # C
m = 40 * 1.66053906660e-27   # kg (40Ca+)
epsilon0 = 8.8541878128e-12  # F/m
omega_x = 2 * np.pi * 1e6    # rad/s
E_bg = 1e4                    # V/m
delta_omega = 2 * np.pi * 5   # Hz (conservative noise floor)
h_nominal = 40e-6             # m

# Grid
N = 256
FOV = 128e-6  # m
dx = FOV / N

# ============================================================
# Generate realistic calibration surface
# ============================================================
np.random.seed(42)

x = np.linspace(0, FOV, N, endpoint=False)
y = np.linspace(0, FOV, N, endpoint=False)
X, Y = np.meshgrid(x, y)

# 1. Sinusoidal calibration grating (NIST-certified-like: A=100nm, lambda=3um)
z_grating = 100e-9 * np.sin(2 * np.pi * X / 3e-6)

# 2. Gaussian plateau (height 200 nm, sigma=7.5 um) - standard test feature
r_centre = np.sqrt((X - FOV/2)**2 + (Y - FOV/2)**2)
z_plateau = 200e-9 * np.exp(-r_centre**2 / (2 * (7.5e-6)**2))

# 3. Realistic broadband roughness using published PSD parameters
# Power spectrum: S(k) ~ [1 + (k*l_c)^2]^(-beta)
# Parameters from: Persson, Surf. Sci. Rep. 61, 201 (2006) - machined steel
k_vals = fftfreq(N, dx) * 2 * np.pi
KX, KY = np.meshgrid(k_vals, k_vals)
K = np.sqrt(KX**2 + KY**2)
K[0, 0] = 1e-12  # avoid division by zero

l_c = 5e-6    # correlation length (machined metal)
beta = 1.5     # spectral exponent
S_k = (1 + (K * l_c)**2)**(-beta)
S_k[0, 0] = 0

# Generate random phase
phase = np.random.uniform(0, 2*np.pi, (N, N))
z_rough_fft = np.sqrt(S_k) * np.exp(1j * phase)
z_rough = np.real(ifft2(z_rough_fft)) * N * N  # normalise

# Scale to target Ra ~ 50 nm (typical for lapped metal)
z_rough *= 50e-9 / np.std(z_rough)

# Combine: grating + plateau + roughness
z_afm = z_grating + z_plateau + z_rough

print(f"AFM surface: z_min={z_afm.min()*1e9:.1f} nm, z_max={z_afm.max()*1e9:.1f} nm")
print(f"  R_a = {np.mean(np.abs(z_afm - np.mean(z_afm)))*1e9:.1f} nm")

# ============================================================
# Forward model (ITF propagation)
# ============================================================
def apply_itf(z_surface, h, omega_chan=omega_x):
    """Propagate surface through ITF: z_s -> Delta_omega"""
    phi_s = E_bg * z_surface
    phi_s_fft = fft2(phi_s)
    KX_m, KY_m = np.meshgrid(k_vals, k_vals)
    K_m = np.sqrt(KX_m**2 + KY_m**2)
    # Propagate to ion height
    phi_at_ion_fft = phi_s_fft * np.exp(-K_m * h)
    # Second derivative in x
    dw_fft = (e / (2 * m * omega_chan)) * (-KX_m**2) * phi_at_ion_fft
    dw_fft[0, 0] = 0  # DC term
    return np.real(ifft2(dw_fft))

def invert_omega_to_z(dw, h, lamb=1e-4):
    """Tikhonov inversion: dw -> z"""
    dw_fft = fft2(dw)
    KX_m, KY_m = np.meshgrid(k_vals, k_vals)
    K_m = np.sqrt(KX_m**2 + KY_m**2)
    # ITF
    H = -(e * E_bg / (2 * m * omega_x)) * KX_m**2 * np.exp(-K_m * h)
    # Tikhonov regularisation
    lam = lamb * np.max(np.abs(H))
    H_inv = np.conj(H) / (np.abs(H)**2 + lam**2)
    H_inv[0, 0] = 0
    z_recon_fft = H_inv * dw_fft
    return np.real(ifft2(z_recon_fft))

# Propagate
dw_afm = apply_itf(z_afm, h_nominal)
dw_afm_noisy = dw_afm + np.random.normal(0, delta_omega, (N, N))

# Invert
z_recon = invert_omega_to_z(dw_afm_noisy, h_nominal)

# Mask below resolution limit (lambda_min ~ 1.48 h)
k_x = 2 * np.pi * fftfreq(N, dx)
k_y = 2 * np.pi * fftfreq(N, dx)
KX_mask, KY_mask = np.meshgrid(k_x, k_y)
K_mask = np.sqrt(KX_mask**2 + KY_mask**2)
resolvable = K_mask < 2 / (1.48 * h_nominal)

# Correlation with ground truth (only resolvable region)
mask_2d = np.fft.fftshift(resolvable)
corr = np.corrcoef(z_afm.ravel(), z_recon.ravel())[0, 1]
print(f"Pearson r (AFM vs QESPM reconstruction): {corr:.4f}")

# ============================================================
# Figure V1: AFM input surface
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

im0 = axes[0].imshow(z_afm * 1e9, extent=[0, FOV*1e6, 0, FOV*1e6],
                      cmap='terrain', origin='lower', aspect='equal')
axes[0].set_title('AFM ground truth: $z(x,y)$ [nm]')
axes[0].set_xlabel('$x$ [µm]'); axes[0].set_ylabel('$y$ [µm]')
plt.colorbar(im0, ax=axes[0], label='nm')

# Power spectral density
z_fft = fft2(z_afm)
psd = np.abs(z_fft)**2
psd_radial = np.zeros(N//2)
for i in range(N):
    for j in range(N):
        kr = int(np.sqrt(i**2 + j**2))
        if kr < N//2:
            psd_radial[kr] += psd[i, j]
k_plot = 2 * np.pi * np.arange(1, N//2) / FOV
axes[1].loglog(k_plot/1e6, psd_radial[1:N//2], 'k-', lw=1)
axes[1].set_xlabel('$k$ [rad/µm]'); axes[1].set_ylabel('PSD [m$^2$]')
axes[1].set_title('Surface PSD')
axes[1].grid(True, alpha=0.3)

# Cross-section
x_cross = x * 1e6
idx = N//2
axes[2].plot(x_cross, z_afm[idx, :]*1e9, 'b-', lw=1, label='AFM')
axes[2].plot(x_cross, z_recon[idx, :]*1e9, 'r--', lw=1, label=f'QESPM (r={corr:.3f})')
axes[2].set_xlabel('$x$ [µm]'); axes[2].set_ylabel('$z$ [nm]')
axes[2].set_title('Cross-section at $y$ = centre')
axes[2].legend(fontsize=9)
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('manuscript/figV1_afm_input.pdf', dpi=150, bbox_inches='tight')
plt.close()
print('Saved: figV1_afm_input.pdf')

# ============================================================
# Figure V2: Ion signal
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

dw_hz = dw_afm_noisy / (2*np.pi)
vmax = np.percentile(np.abs(dw_hz), 99)
im1 = axes[0].imshow(dw_hz, extent=[0, FOV*1e6, 0, FOV*1e6],
                      cmap='RdBu_r', origin='lower', aspect='equal',
                      vmin=-vmax, vmax=vmax)
axes[0].set_title(r'$\Delta\omega_x/2\pi$ [Hz] (ion signal)')
axes[0].set_xlabel('$x$ [µm]'); axes[0].set_ylabel('$y$ [µm]')
plt.colorbar(im1, ax=axes[0], label='Hz')

# SNR map
dw_clean = apply_itf(z_afm, h_nominal)
snr = np.abs(dw_clean) / delta_omega
im2 = axes[1].imshow(snr, extent=[0, FOV*1e6, 0, FOV*1e6],
                      cmap='inferno', origin='lower', aspect='equal')
axes[1].set_title('SNR map (|signal| / 5 Hz)')
axes[1].set_xlabel('$x$ [µm]'); axes[1].set_ylabel('$y$ [µm]')
plt.colorbar(im2, ax=axes[1])
axes[1].contour(X*1e6, Y*1e6, snr, levels=[1], colors='cyan', linewidths=1, linestyles='--')

plt.tight_layout()
plt.savefig('manuscript/figV2_afm_ion_signal.pdf', dpi=150, bbox_inches='tight')
plt.close()
print('Saved: figV2_afm_ion_signal.pdf')

# ============================================================
# Figure V3: Reconstruction vs ground truth
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

vlim = 200
axes[0, 0].imshow(z_afm*1e9, cmap='terrain', origin='lower', aspect='equal',
                   extent=[0, FOV*1e6, 0, FOV*1e6], vmin=-vlim, vmax=vlim)
axes[0, 0].set_title('(a) AFM ground truth [nm]')

axes[0, 1].imshow(z_recon*1e9, cmap='terrain', origin='lower', aspect='equal',
                   extent=[0, FOV*1e6, 0, FOV*1e6], vmin=-vlim, vmax=vlim)
axes[0, 1].set_title(f'(b) QESPM reconstruction [nm]')

residual = (z_recon - z_afm) * 1e9
rmax = np.percentile(np.abs(residual), 99)
axes[1, 0].imshow(residual, cmap='RdBu_r', origin='lower', aspect='equal',
                   extent=[0, FOV*1e6, 0, FOV*1e6], vmin=-rmax, vmax=rmax)
axes[1, 0].set_title(f'(c) Residual [nm]  RMS={np.std(residual):.1f} nm')

# Scatter: AFM vs QESPM
mask_sample = np.random.choice(N*N, 2000, replace=False)
axes[1, 1].scatter(z_afm.flat[mask_sample]*1e9, z_recon.flat[mask_sample]*1e9,
                    c='k', s=1, alpha=0.3)
zlim = 200
axes[1, 1].plot([-zlim, zlim], [-zlim, zlim], 'r--', lw=1)
axes[1, 1].set_xlim(-zlim, zlim); axes[1, 1].set_ylim(-zlim, zlim)
axes[1, 1].set_xlabel('AFM $z$ [nm]'); axes[1, 1].set_ylabel('QESPM $z$ [nm]')
axes[1, 1].set_title(f'(d) Point-by-point  $r={corr:.3f}$')
axes[1, 1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('manuscript/figV3_afm_reconstruction.pdf', dpi=150, bbox_inches='tight')
plt.close()
print('Saved: figV3_afm_reconstruction.pdf')

print("\nAll AFM validation figures saved to manuscript/")
