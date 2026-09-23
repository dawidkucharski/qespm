#!/usr/bin/env python3
"""
Supplementary figures V1-V3 from the REAL Camargo AFM data (20x20 um).

Consistent with real_afm_benchmark.py (main-text Fig. V9).

Output files (in manuscript/):
  figV1_afm_input.pdf          — real AFM topography + PSD
  figV2_afm_ion_signal.pdf     — simulated ion signal at h = 2, 5, 10, 40 um
  figV3_afm_reconstruction.pdf — reconstructions + residuals at 4 heights
"""

import numpy as np
from numpy.fft import fft2, ifft2, fftfreq
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import os, struct

# ============================================================
# Load & parse NanoScope AFM data (Camargo Jr. et al. 2019)
# ============================================================
nc_file = os.path.join('data', 'Raw AFM NanoScope images', 'SS 20um.002')
with open(nc_file, 'rb') as f:
    raw = f.read()

marker = b'\\*File list end'
idx = raw.find(marker)
data_start = idx + len(marker)
while data_start < len(raw) and raw[data_start:data_start+1] in (b'\r', b'\n', b'\x00', b'\x1a'):
    data_start += 1

N_px = 512
data_bytes = raw[data_start:data_start + N_px*N_px*2]
raw_vals = struct.unpack('<' + 'h'*(len(data_bytes)//2), data_bytes[:N_px*N_px*2])
z_raw = np.array(raw_vals).reshape(N_px, N_px).astype(np.float64)
z_range_nm = 500.0
z_afm = z_raw * (z_range_nm / (z_raw.max() - z_raw.min()))
z_afm -= z_afm.mean()
print(f"Camargo stainless steel: {N_px}x{N_px}, "
      f"{z_range_nm:.0f} nm p-v, R_a={np.std(z_afm):.0f} nm")

# ============================================================
# QESPM parameters
# ============================================================
e = 1.602176634e-19
m = 40 * 1.66053906660e-27
omega_x = 2 * np.pi * 1e6
E_bg = 1e4
delta_omega = 2 * np.pi * 5

scan_um = 20.0
FOV = scan_um * 1e-6
N = N_px
dx = FOV / N

k_vals = 2*np.pi*fftfreq(N, dx)
KX, KY = np.meshgrid(k_vals, k_vals)
K = np.sqrt(KX**2+KY**2)
K[0, 0] = 1e-12


def forward(z_s, h_val):
    phi_fft = fft2(E_bg*z_s)*np.exp(-K*h_val)
    dw_fft = (e/(2*m*omega_x))*(-KX**2)*phi_fft
    dw_fft[0, 0] = 0
    return np.real(ifft2(dw_fft))


def invert(dw, h_val, lamb=1e-4):
    dw_fft = fft2(dw)
    H = -(e*E_bg/(2*m*omega_x))*KX**2*np.exp(-K*h_val)
    lam = lamb*np.max(np.abs(H))
    H_inv = np.conj(H)/(np.abs(H)**2+lam**2)
    H_inv[0, 0] = 0
    return np.real(ifft2(H_inv*dw_fft))


# ============================================================
# Figure V1: real AFM input + PSD
# ============================================================
fig1, axes1 = plt.subplots(1, 2, figsize=(13, 5.5))
vlim = np.percentile(np.abs(z_afm), 99)

im0 = axes1[0].imshow(z_afm, extent=[0, scan_um, 0, scan_um], cmap='terrain',
                      origin='lower', aspect='equal', vmin=-vlim, vmax=vlim)
axes1[0].set_title('(a) AFM topography — stainless steel [nm]')
axes1[0].set_xlabel('x [µm]'); axes1[0].set_ylabel('y [µm]')
plt.colorbar(im0, ax=axes1[0], fraction=0.046)

z_fft = np.fft.fftshift(fft2(z_afm))
psd_r = np.zeros(N//2)
for i in range(N):
    for j in range(N):
        kr = int(np.sqrt((i-N//2)**2+(j-N//2)**2))
        if 0 < kr < N//2:
            psd_r[kr] += np.abs(z_fft[i, j])**2
k_plot = 2*np.pi*np.arange(1, N//2)/FOV
axes1[1].loglog(k_plot/1e6, psd_r[1:], 'k-', lw=1.5)
axes1[1].set_xlabel('$k$ [rad/µm]')
axes1[1].set_ylabel('PSD [nm²]')
axes1[1].set_title('(b) Surface power spectral density')
axes1[1].grid(True, alpha=0.3)

fig1.suptitle('Validation Input: Camargo Jr. et al. (2019), Mendeley Data, '
              'doi:10.17632/6dzmrjngcg.3 (CC BY 4.0)',
              fontsize=11, fontweight='bold')
plt.tight_layout()
fig1.savefig('manuscript/figV1_afm_input.pdf', dpi=150, bbox_inches='tight')
plt.close()
print('Saved: figV1_afm_input.pdf')

# ============================================================
# Figure V2: simulated ion signal at 4 heights
# ============================================================
heights = [2, 5, 10, 40]
fig2, axes2 = plt.subplots(2, 2, figsize=(11, 10))
sig_rms = {}
for h_um, ax in zip(heights, axes2.ravel()):
    dw = forward(z_afm, h_um*1e-6)
    sig_rms[h_um] = np.std(dw)
    vmax_d = np.percentile(np.abs(dw), 99)
    im = ax.imshow(dw/(2*np.pi), extent=[0, scan_um, 0, scan_um], cmap='coolwarm',
                   origin='lower', aspect='equal', vmin=-vmax_d/(2*np.pi),
                   vmax=vmax_d/(2*np.pi))
    s_val = sig_rms[h_um] / (2*np.pi)
    exp10 = int(np.floor(np.log10(s_val)))
    mant = s_val / 10**exp10
    ax.set_title(f'$h$ = {h_um} µm  ($\sigma_{{\Delta\omega}}/2\pi$ = '
                 f'${mant:.1f}\times10^{{{exp10}}}$ Hz)', fontsize=10)
    ax.set_xlabel('x [µm]', fontsize=8); ax.set_ylabel('y [µm]', fontsize=8)
    if h_um <= 5:
        ax.text(0.5, 0.03, 'asymptotic ITF prediction --- outside validity domain',
                transform=ax.transAxes, ha='center', va='bottom', fontsize=8,
                color='white', bbox=dict(facecolor='black', alpha=0.6))
    plt.colorbar(im, ax=ax, fraction=0.046)

fig2.suptitle('Simulated ion signal $\\Delta\\omega_x(x,y)$ from AFM topography',
              fontsize=13, fontweight='bold')
plt.tight_layout()
fig2.savefig('manuscript/figV2_afm_ion_signal.pdf', dpi=150, bbox_inches='tight')
plt.close()
print('Saved: figV2_afm_ion_signal.pdf')
print(f"Signal RMS: h=2um {sig_rms[2]/(2*np.pi):.0f} Hz, "
      f"h=5um {sig_rms[5]/(2*np.pi):.0f} Hz, h=10um {sig_rms[10]/(2*np.pi):.0f} Hz, "
      f"h=40um {sig_rms[40]/(2*np.pi):.0f} Hz, "
      f"factor 2->40um: {sig_rms[2]/sig_rms[40]:.0f}x")

# ============================================================
# Figure V3: reconstructions at 4 heights + residuals
# ============================================================
fig3, axes3 = plt.subplots(2, 4, figsize=(18, 9))
r_vals = {}
for col, h_um in enumerate([40, 10, 5, 2]):
    dw = forward(z_afm, h_um*1e-6) + np.random.normal(0, delta_omega, (N, N))
    zr = invert(dw, h_um*1e-6)
    r = np.corrcoef(z_afm.ravel(), zr.ravel())[0, 1]
    r_vals[h_um] = r
    im_t = axes3[0, col].imshow(zr, extent=[0, scan_um, 0, scan_um], cmap='terrain',
                                origin='lower', aspect='equal', vmin=-vlim, vmax=vlim)
    axes3[0, col].set_title(f'({chr(97+col)}) Recon., $h$={h_um} µm ($r$={r:.2f})',
                            fontsize=10)
    axes3[0, col].set_xlabel('x [µm]', fontsize=8)
    if col == 0:
        axes3[0, col].set_ylabel('y [µm]', fontsize=8)
    resid = zr - z_afm
    vmax_r = np.percentile(np.abs(resid), 99)
    im_r = axes3[1, col].imshow(resid, extent=[0, scan_um, 0, scan_um], cmap='RdBu_r',
                                origin='lower', aspect='equal',
                                vmin=-vmax_r, vmax=vmax_r)
    axes3[1, col].set_title(f'({chr(101+col)}) Residual, $h$={h_um} µm', fontsize=10)
    axes3[1, col].set_xlabel('x [µm]', fontsize=8)
    if col == 0:
        axes3[1, col].set_ylabel('y [µm]', fontsize=8)

fig3.suptitle('QESPM reconstruction from AFM topography (Camargo data)',
              fontsize=13, fontweight='bold')
plt.tight_layout()
fig3.savefig('manuscript/figV3_afm_reconstruction.pdf', dpi=150, bbox_inches='tight')
plt.close()
print('Saved: figV3_afm_reconstruction.pdf')
print(f"Correlations: h=40 {r_vals[40]:.3f}, h=10 {r_vals[10]:.3f}, "
      f"h=5 {r_vals[5]:.3f}, h=2 {r_vals[2]:.3f}")
