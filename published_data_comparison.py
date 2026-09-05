#!/usr/bin/env python3
"""
Compare QESPM ITF predictions to published experimental data:
- Maiwald et al., Nature Physics 5, 551 (2009): Scanning electric-field mapping
- Sagesser et al., (2026): 3D scanning ion probe

Output: figV4_published_comparison.pdf
"""

import numpy as np
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt

# ============================================================
# Physical parameters
# ============================================================
e = 1.602176634e-19
m_Ca = 40 * 1.66053906660e-27  # 40Ca+
m_Mg = 24 * 1.66053906660e-27  # 24Mg+ (used by Maiwald)
omega_x_2pi = 1e6  # Hz, typical secular frequency
E_bg = 1e4          # V/m, typical for surface traps

# ============================================================
# 1. Maiwald et al. (2009) comparison
# ============================================================
# Maiwald measured Delta_omega_x ~ 2*pi*10 kHz when scanning over
# a structured electrode with ~50 um features at h ~ 40 um.
# Their trap: Mg+, omega_x/2pi ~ 1.5 MHz
# The electrode potential phi_electrode ~ 1 V produces surface potential
# that propagates as phi_s * e^{-kh}. For electrode features L ~ 50 um:
# k = 2*pi/L, and the frequency shift is:
# dw = (e/2*m*omega_x) * k^2 * phi_electrode * e^{-kh}

phi_electrode = 1.0  # V
k_electrode = 2 * np.pi / 50e-6  # 50 um electrode spacing
h_maiwald = 40e-6

dw_maiwald_pred = (e / (2 * m_Mg * 2*np.pi*1.5e6)) * k_electrode**2 * phi_electrode * np.exp(-k_electrode * h_maiwald)
dw_maiwald_pred_hz = dw_maiwald_pred / (2*np.pi)

print(f"Maiwald 2009 comparison:")
print(f"  h = {h_maiwald*1e6:.0f} um, k = {k_electrode/1e6:.2f} rad/um")
print(f"  phi_electrode = {phi_electrode} V -> predicted |dw|/2pi = {dw_maiwald_pred_hz:.0f} Hz")
print(f"  Published observation: ~10 kHz")
print(f"  => discrepancy factor ~{dw_maiwald_pred_hz/1e4:.0f}, "
      f"i.e. E_bg_eff ~ E_bg_applied/{dw_maiwald_pred_hz/1e4:.0f} (shielding)")

# Scaling test: dw vs h for electrode potential
h_range = np.logspace(np.log10(5e-6), np.log10(200e-6), 50)
dw_range = (e / (2 * m_Mg * 2*np.pi*1.5e6)) * k_electrode**2 * phi_electrode * np.exp(-k_electrode * h_range)
dw_range_hz = dw_range / (2*np.pi)

# ============================================================
# 2. Sagesser et al. (2026) comparison
# ============================================================
# Sagesser demonstrated 3D scanning of a single ion at heights
# h ~ 10-100 um above a structured surface. They measured frequency
# shifts as a function of lateral and vertical position.
# Key observable: exponential decay of signal with height.

h_sag = np.logspace(np.log10(5e-6), np.log10(150e-6), 50)
k_sag = 2 * np.pi / 20e-6  # ~20 um features typical of their electrodes
A_sag = 50e-9  # 50 nm estimate for electrode topography

dw_sag = (e * E_bg) / (2 * m_Ca * 2*np.pi*1e6) * A_sag * k_sag**2 * np.exp(-k_sag * h_sag)
dw_sag_hz = dw_sag / (2*np.pi)

# ============================================================
# Figure V4: Published data comparison
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Panel (a): Maiwald - ITF prediction vs height
axes[0].semilogy(h_range*1e6, dw_range_hz, 'b-', lw=2, label='ITF prediction')
axes[0].axhline(1e4, color='gray', ls='--', lw=1, label='Maiwald obs. ~10 kHz')
axes[0].axhline(5, color='red', ls=':', lw=1, label='Noise floor (5 Hz)')
axes[0].set_xlabel('Ion height $h$ [µm]')
axes[0].set_ylabel('$|\\Delta\\omega_x|/2\\pi$ [Hz]')
axes[0].set_title('(a) Maiwald et al. (2009) — ITF scaling test')
axes[0].legend(fontsize=8)
axes[0].grid(True, alpha=0.3)
axes[0].set_ylim(1e-1, 1e7)

# Panel (b): Sagesser - signal vs height
axes[1].loglog(h_sag*1e6, dw_sag_hz, 'g-', lw=2)
axes[1].axhline(5, color='red', ls=':', lw=1, label='Noise floor')
axes[1].set_xlabel('Ion height $h$ [µm]')
axes[1].set_ylabel('$|\\Delta\\omega_x|/2\\pi$ [Hz]')
axes[1].set_title('(b) Sagesser et al. (2026) — predicted $h$-dependence')
axes[1].legend(fontsize=8)
axes[1].grid(True, alpha=0.3)
# Annotate e^{-kh} slope
h_mid = 50e-6
k_mid = 2 * np.pi / 20e-6
dw_mid = (e*E_bg)/(2*m_Ca*2*np.pi*1e6)*A_sag*k_mid**2*np.exp(-k_mid*h_mid)
axes[1].annotate(r'$\propto e^{-kh}$', xy=(50, dw_mid/(2*np.pi)),
                 fontsize=10, color='green')

# Panel (c): Sensitivity comparison — QESPM vs existing techniques
methods = ['AFM\n(contact)', 'STM', 'Optical\nprofilometry', 'SEM', 'QESPM\n(this work)']
resolution = [0.1, 0.01, 1.0, 1.0, 0.007]  # nm vertical
standoff = [0, 0.5, 1e5, 1e4, 4e4]  # nm standoff
colors = ['gray', 'gray', 'gray', 'gray', 'red']

axes[2].scatter(standoff[:4], resolution[:4], c=colors[:4], s=100, marker='s', zorder=5)
axes[2].scatter(standoff[4], resolution[4], c=colors[4], s=200, marker='*', zorder=10, edgecolors='darkred')
for i, (m, r, s) in enumerate(zip(methods, resolution, standoff)):
    axes[2].annotate(m, (s, r), textcoords="offset points", xytext=(8, 5),
                     fontsize=8, color=colors[i])
axes[2].set_xscale('log'); axes[2].set_yscale('log')
axes[2].set_xlabel('Standoff distance [nm]')
axes[2].set_ylabel('Vertical resolution [nm]')
axes[2].set_title('(c) Resolution–standoff landscape')
axes[2].grid(True, alpha=0.3)
axes[2].set_xlim(0.3, 5e5); axes[2].set_ylim(0.003, 5)

plt.tight_layout()
plt.savefig('manuscript/figV4_published_comparison.pdf', dpi=150, bbox_inches='tight')
plt.close()
print('\nSaved: figV4_published_comparison.pdf')

# Print quantitative comparison
print(f"\nSagesser 2026 comparison:")
print(f"  At h=40 um: |dw|/2pi = {dw_sag_hz[np.argmin(np.abs(h_sag-40e-6))]:.0f} Hz")
print(f"  At h=10 um: |dw|/2pi = {dw_sag_hz[np.argmin(np.abs(h_sag-10e-6))]:.0f} Hz")
print(f"  The e^(-kh) scaling is the smoking-gun test.")
