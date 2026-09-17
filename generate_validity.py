#!/usr/bin/env python3
"""Generate the QESPM forward-model domain-of-validity diagram (figD1_validity.pdf).

Regime map in the (lambda, h) plane, following Sec. 2.6 of the main text:
  - resolvable band: 1.51 h < lambda < L (FOV), above the large-slope limit
  - unresolvable:      lambda < 1.51 h
  - large-slope breakdown: k A >= 0.45  (A = 100 nm)  ->  lambda <= 2*pi*A/0.45
  - heating-limited:   h < 5 um
  - Casimir-Polder exclusion: h < 0.1 um (100 nm), conservative per Assumption 7
  - ITF peak:          lambda = pi*h
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

lam = np.logspace(-1, 2, 500)          # wavelength [um], 0.1 .. 100
h = np.logspace(-1, 2, 500)            # ion height [um], 0.1 .. 100
LAM, H = np.meshgrid(lam, h)

lam_hi = 64.0                          # field of view [um]
lam_slope = 2 * np.pi * 0.1 / 0.45     # large-slope limit for A = 100 nm (~1.4 um)
h_heat = 5.0                           # heating-dominated below 5 um
h_cp = 0.1                             # Casimir-Polder exclusion below 100 nm

fig, ax = plt.subplots(figsize=(8.6, 6.4))

# --- heating-limited band (h < 5 um), excluding the C-P strip below 0.1 um ---
ax.axhspan(h.min(), h_heat, xmin=0, xmax=1, color='#f4c7c3', alpha=0.55, zorder=1)

# --- large-slope breakdown (lambda < 1.4 um, all h) ---
ax.axvspan(lam.min(), lam_slope, color='#fde0b2', alpha=0.85, zorder=1)

# --- resolvable band: max(1.51 h, lam_slope) < lambda < lam_hi ---
lam_lo = np.maximum(1.51 * h, lam_slope)
ax.fill_between(h, lam_lo, lam_hi, where=(lam_lo < lam_hi),
                color='#c8e6c9', alpha=0.7, zorder=1)

# --- Casimir-Polder exclusion strip along the bottom edge (h < 0.1 um) ---
ax.axhspan(h.min(), h_cp, color='#b0bec5', alpha=0.85, zorder=3)

# --- FOV limit at lambda = L ---
ax.axvline(lam_hi, color='grey', ls=':', lw=2.0, zorder=2)

# --- ITF peak line lambda = pi h ---
ax.plot(lam, lam / np.pi, 'b--', lw=2.0, zorder=4, label=r'ITF peak $\lambda=\pi h$')
ax.legend(loc='lower right', fontsize=10, framealpha=0.9)

# --- labels ---
ax.text(0.13, 1.5, 'large-slope\nbreakdown\n($kA\geq0.45$, $A=100$\,nm)',
        fontsize=9, ha='left', va='top', color='#7a4b00')
ax.text(9.0, 40.0, 'resolvable band\n($\\lambda>1.51h$, $\\lambda<L$)',
        fontsize=10, ha='center', va='center', color='#1b5e20')
ax.text(0.35, 0.35, 'unresolvable\n($\\lambda<1.51h$)',
        fontsize=10, ha='center', va='center', color='#5d4037')
ax.text(28.0, 1.6, 'heating-limited ($h<5$ $\mu$m)',
        fontsize=10, ha='center', va='center', color='#8c1d18')
ax.text(70.0, 60.0, 'FOV-limited\n$\\lambda>L$', fontsize=9,
        ha='center', va='center', color='#37474f', rotation=90)
ax.text(0.13, 0.105, 'Casimir--Polder exclusion ($h<100$\,nm)',
        fontsize=9, ha='left', va='bottom', color='#263238')

ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(0.1, 100)
ax.set_ylim(0.1, 100)
ax.set_xlabel('surface feature wavelength $\lambda$ [$\mu$m]', fontsize=12)
ax.set_ylabel('ion--surface separation $h$ [$\mu$m]', fontsize=12)
ax.set_title('Domain of validity of the QESPM forward model', fontsize=13)
ax.grid(True, which='both', alpha=0.25)

plt.tight_layout()
plt.savefig('manuscript/figD1_validity.pdf', dpi=150, bbox_inches='tight')
plt.close()
print('saved manuscript/figD1_validity.pdf')
