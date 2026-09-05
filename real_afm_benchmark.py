#!/usr/bin/env python3
"""
QESPM validation against published AFM data of stainless steel.
Data: Camargo Jr. et al. (2019), Mendeley Data, V3, doi:10.17632/6dzmrjngcg.3
CC BY 4.0 license.

Surface: Stainless steel, 20x20 um scan, 512x512 px, Bruker NanoScope AFM.
Conducting surface -- ideal for QESPM perfect-conductor assumption.
"""

import numpy as np
from numpy.fft import fft2, ifft2, fftfreq
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import os, struct

# ============================================================
# Load & parse NanoScope AFM data
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

# Normalize: raw range -> ~500 nm p-v (typical for steel AFM)
z_range_nm = 500.0
z_afm = z_raw * (z_range_nm / (z_raw.max() - z_raw.min()))
z_afm -= z_afm.mean()

print(f"Stainless steel AFM: {N_px}x{N_px}, {z_range_nm:.0f} nm p-v, R_a={np.std(z_afm):.0f} nm")

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
N = N_px; dx = FOV/N

k_vals = 2*np.pi*fftfreq(N, dx)
KX, KY = np.meshgrid(k_vals, k_vals)
K = np.sqrt(KX**2+KY**2); K[0,0]=1e-12

# ============================================================
# ITF
# ============================================================
def forward(z_s, h_val):
    phi_fft = fft2(E_bg*z_s)*np.exp(-K*h_val)
    return np.real(ifft2((e/(2*m*omega_x))*(-KX**2)*phi_fft))

def invert(dw, h_val, lamb=1e-4):
    dw_fft = fft2(dw)
    H = -(e*E_bg/(2*m*omega_x))*KX**2*np.exp(-K*h_val)
    lam = lamb*np.max(np.abs(H))
    H_inv = np.conj(H)/(np.abs(H)**2+lam**2); H_inv[0,0]=0
    return np.real(ifft2(H_inv*dw_fft))

def corr_resolved(zt, zr, kc):
    m = K < kc
    if m.sum() < 5: return np.nan
    return np.corrcoef(np.real(ifft2(fft2(zt)*m)).ravel(),
                        np.real(ifft2(fft2(zr)*m)).ravel())[0,1]

# ============================================================
# Multi-height
# ============================================================
h_um = np.array([2, 5, 10, 20, 40])
h_list = h_um*1e-6
r_full = np.zeros(len(h_list)); r_res = np.zeros(len(h_list))
# 50% MTF cutoff: solve (kh/2)^2 exp(2-kh) = 0.5  ->  kh_50 ~ 4.156
from scipy.optimize import brentq
kh_50 = brentq(lambda x: (x/2)**2*np.exp(2-x) - 0.5, 2, 6)
kc = kh_50/h_list

print(f"\n{'h [um]':>8} {'lam_min [um]':>12} {'r_full':>8} {'r_res':>8}")
print("-"*42)
for i, hv in enumerate(h_list):
    dw = forward(z_afm, hv) + np.random.normal(0, delta_omega, (N,N))
    zr = invert(dw, hv)
    r_full[i] = np.corrcoef(z_afm.ravel(), zr.ravel())[0,1]
    r_res[i] = corr_resolved(z_afm, zr, kc[i])
    print(f"{hv*1e6:8.0f} {2*np.pi/kc[i]*1e6:12.1f} {r_full[i]:8.4f} {r_res[i]:8.4f}")

# ============================================================
# Figure
# ============================================================
fig, axes = plt.subplots(2, 3, figsize=(18, 11))
x_um = np.linspace(0, scan_um, N)
vlim = np.percentile(np.abs(z_afm), 99)

axes[0,0].imshow(z_afm, extent=[0,scan_um,0,scan_um], cmap='terrain', origin='lower', aspect='equal', vmin=-vlim, vmax=vlim)
axes[0,0].set_title('(a) AFM -- stainless steel [nm]')

axes[0,1].plot(h_um, r_full, 'ko-', ms=10, lw=2, label='Full bandwidth')
axes[0,1].plot(h_um, r_res, 'bs--', ms=10, lw=2, label=r'Resolved ($k<k_{\rm cutoff}$)')
axes[0,1].set_xlabel('Ion height $h$ [um]'); axes[0,1].set_ylabel('Pearson $r$')
axes[0,1].set_title('(b) Reconstruction quality vs $h$')
axes[0,1].legend(fontsize=9); axes[0,1].grid(True, alpha=0.3); axes[0,1].set_ylim(0, 1.05)
for i in range(len(h_um)): axes[0,1].annotate(f'{r_full[i]:.2f}', (h_um[i], r_full[i]+0.03), ha='center', fontsize=9)

z_fft = np.fft.fftshift(fft2(z_afm))
psd_r = np.zeros(N//2)
for i in range(N):
    for j in range(N):
        kr = int(np.sqrt((i-N//2)**2+(j-N//2)**2))
        if 0<kr<N//2: psd_r[kr] += np.abs(z_fft[i,j])**2

k_plot = 2*np.pi*np.arange(1,N//2)/FOV
colors = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd']
axes[0,2].loglog(k_plot/1e6, psd_r[1:], 'k-', lw=1.5, label='AFM PSD')
for i, hv in enumerate(h_list):
    axes[0,2].axvline(kh_50/hv/1e6, color=colors[i], ls='--', lw=1.2, label=f'$h$={h_um[i]} um')
axes[0,2].set_xlabel('$k$ [rad/um]'); axes[0,2].set_ylabel('PSD [nm^2]')
axes[0,2].set_title('(c) PSD with ITF cutoffs'); axes[0,2].legend(fontsize=7, ncol=2); axes[0,2].grid(True, alpha=0.3)

for col, (ih, label) in enumerate([(1,5),(3,20)]):
    dw = forward(z_afm, h_list[ih]) + np.random.normal(0, delta_omega, (N,N))
    zr = invert(dw, h_list[ih])
    axes[1,col].imshow(zr, extent=[0,scan_um,0,scan_um], cmap='terrain', origin='lower', aspect='equal', vmin=-vlim, vmax=vlim)
    axes[1,col].set_title(f'({"de"[col]}) QESPM $h$={label} um  ($r$={r_full[ih]:.3f})')

mid = N//2
dw5 = forward(z_afm, h_list[1]) + np.random.normal(0, delta_omega, (N,N))
z5 = invert(dw5, h_list[1])
axes[1,2].plot(x_um, z_afm[mid,:], 'b-', lw=1, alpha=0.7, label='AFM')
axes[1,2].plot(x_um, z5[mid,:], 'r--', lw=1, alpha=0.7, label=f'QESPM ($r$={r_full[1]:.3f})')
axes[1,2].set_xlabel('$x$ [um]'); axes[1,2].set_ylabel('$z$ [nm]')
axes[1,2].set_title('(f) Cross-section, $h$ = 5 um'); axes[1,2].legend(fontsize=8); axes[1,2].grid(True, alpha=0.3)

fig.text(0.5, 0.01, 'AFM data: Camargo Jr. et al. (2019), Mendeley Data, doi:10.17632/6dzmrjngcg.3 (CC BY 4.0). Stainless steel, 20x20 um, Bruker NanoScope.',
         ha='center', fontsize=8, style='italic')
plt.tight_layout(rect=[0,0.03,1,1])
plt.savefig('manuscript/figV9_real_afm_benchmark.pdf', dpi=150, bbox_inches='tight')
plt.close()
print('\nSaved: manuscript/figV9_real_afm_benchmark.pdf')
print(f'\nKEY: r improves from {r_full[-1]:.2f} (h=40um) to {r_full[0]:.2f} (h=2um)')
print('Validates ITF height-scaling on independently measured AFM data of a conducting surface.')
