#!/usr/bin/env python3
"""Generate the eight manuscript figures that lacked generator scripts.

Outputs (into manuscript/):
  fig1a_itf.pdf           ITF magnitude vs k for several heights
  fig6_qfi.pdf            QFI detection sensitivity and enhancement factor
  figI9_joint_bayes.pdf   joint Bayesian identifiability (319 nm floor)
  figI7_backus_gilbert.pdf Backus-Gilbert averaging kernels and trade-off
  fig3_amin.pdf           A_min(k;h) and reconstructable wavelength band
  fig5_uncertainty.pdf    uncertainty budget decomposition vs k
  fig8_micromotion.pdf    secular vs micromotion channel
  fig9_degeneracy.pdf     two-height ratio (rank-1 degeneracy)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

E = 1.602176634e-19
M = 40 * 1.66053906660e-27
W = 2 * np.pi * 1e6
EBG = 1e4
C = E * EBG / (2 * M * W)          # 1.92e3 m/s
EPS0 = 8.8541878128e-12
H_NOM = 40e-6
SIG_CHARGE = 1.13e-8               # C/m2 budget prior
K_BUDGET = 2e5                     # rad/m evaluation point

def amin(k, dw, h=H_NOM):
    with np.errstate(over='ignore'):
        return (2 * M * W / (E * EBG)) * dw * np.exp(k * h) / k**2

# ------------------------------------------------------------------
# fig1a_itf.pdf
# ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8.6, 5.6))
kum = np.logspace(-2, 1, 400)          # k in rad/um
for h_um in [5, 10, 20, 40, 80, 100]:
    h = h_um * 1e-6
    k = kum * 1e6
    ax.plot(kum, C * k**2 * np.exp(-k * h), lw=1.8,
            label=f'$h={h_um}$ $\\mu$m')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('spatial wavenumber $k$ [rad/$\\mu$m]', fontsize=11)
ax.set_ylabel('$|H_x(k;h)|$ [s$^{-1}$]', fontsize=11)
ax.set_title('Instrument transfer function magnitude $|H_x(k;h)| = C\\,k^2 e^{-kh}$',
             fontsize=12)
ax.axvline(2 / 40, color='k', ls='--', lw=1.2)
ax.text(2 / 40, 3e8, '  $k_{\\rm opt}=2/h$', fontsize=9, va='top')
ax.text(1.51 / 40 * 1e-2, 2e7,
        '$\\lambda_{\\min}^{\\rm(MTF)}\\approx1.51h=60$ $\\mu$m ($h=40$ $\\mu$m)',
        fontsize=9, rotation=0)
ax.legend(fontsize=9, ncol=2)
ax.grid(True, which='both', alpha=0.25)
plt.tight_layout()
plt.savefig('manuscript/fig1a_itf.pdf', dpi=150, bbox_inches='tight')
plt.close()

# ------------------------------------------------------------------
# fig6_qfi.pdf
# ------------------------------------------------------------------
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.4, 4.4))
k = np.logspace(-2, 0, 300) * 1e6      # 0.01 .. 1 rad/um (h=40 um passband)
floors = [('coherent (QPN, 0.7 Hz)', 2 * np.pi * 0.7, '#1f77b4'),
          ('squeezed ($r=1.5$, $\\times 4.5$)', 2 * np.pi * 0.7 / 4.5, '#d62728'),
          ('NOON ($N=100$, $\\times 10$)', 2 * np.pi * 0.7 / 10, '#2ca02c')]
for lab, dw, col in floors:
    a1.plot(k / 1e6, amin(k, dw) * 1e9, lw=2.0, color=col, label=lab)
k_opt = 2 / H_NOM / 1e6               # 2/h = 0.05 rad/um
a1.axvline(k_opt, color='k', ls='--', lw=1.2)
a1.text(k_opt * 1.02, 3e-3,
        '$k_{\\rm opt}=2/h=0.05$ rad/$\\mu$m\n($\\delta A_{\\min}^{\\rm coh}'
        '\\approx 7\\times10^{-3}$ nm)',
        fontsize=8, va='bottom')
a1.set_xscale('log'); a1.set_yscale('log')
a1.set_xlabel('spatial frequency $k$ [rad/$\\mu$m]', fontsize=10)
a1.set_ylabel('$\\delta A_{\\min}$ [nm]', fontsize=10)
a1.set_title('(a) Quantum-limited surface\namplitude sensitivity ($h=40$ $\\mu$m)',
             fontsize=10)
a1.legend(fontsize=8)
a1.grid(True, which='both', alpha=0.25)

a2.plot(k / 1e6, np.full_like(k, 4.5), lw=2.0, color='#d62728',
        label='squeezing ($\\times 4.5$)')
a2.plot(k / 1e6, np.full_like(k, 10.0), lw=2.0, color='#2ca02c',
        label='Heisenberg ($\\times\\sqrt{N}=10$)')
a2.set_xscale('log')
a2.set_xlabel('spatial frequency $k$ [rad/$\\mu$m]', fontsize=10)
a2.set_ylabel('$\\delta A_{\\rm coh}/\\delta A_{\\rm q}$', fontsize=10)
a2.set_title('(b) Quantum enhancement factor', fontsize=10)
a2.legend(fontsize=8)
a2.grid(True, which='both', alpha=0.25)
plt.tight_layout()
plt.savefig('manuscript/fig6_qfi.pdf', dpi=150, bbox_inches='tight')
plt.close()

# ------------------------------------------------------------------
# figI9_joint_bayes.pdf
# ------------------------------------------------------------------
phi_rms = SIG_CHARGE / (2 * EPS0 * K_BUDGET)      # 3.19 mV
p_phi = phi_rms**2
pz_list = [(0.32e-6, 'broad'), (3.2e-6, 'wider')]
s = np.logspace(-2, 8, 400)
fig, (b1, b2) = plt.subplots(1, 2, figsize=(10.4, 4.4))
for pz, lab in pz_list:
    var = (s * p_phi / EBG**2 + 1) / (s + s * p_phi / (EBG**2 * pz) + 1 / pz)
    b1.plot(s, np.sqrt(var) * 1e9, lw=1.8, label=f'$p_z^{{1/2}}={pz*1e6}$ $\\mu$m')
b1.axhline(319, color='k', ls='--', lw=1.4)
b1.text(3e-1, 370, 'charge floor $\\sqrt{p_\\varphi}/E_{\\rm bg}=319$ nm',
        fontsize=9)
b1.set_xscale('log'); b1.set_yscale('log')
b1.set_xlabel('data strength $|H|^2/\\sigma^2$', fontsize=10)
b1.set_ylabel('marginal std of $\\hat z_k$ [nm]', fontsize=10)
b1.set_title('(a) Marginal posterior std of the\ntopography mode (saturation at 319 nm)',
             fontsize=10)
b1.legend(fontsize=8)
b1.grid(True, which='both', alpha=0.25)

c_j = 1 / (2 * EPS0 * K_BUDGET * EBG)
J = np.array([1.0, c_j])
std1 = np.zeros_like(s); std2 = np.zeros_like(s)
for i, si in enumerate(s):
    P = si * np.outer(J, J) + np.diag([1 / (3.2e-6)**2, 1 / p_phi])
    ev = np.linalg.eigvalsh(P)
    std1[i] = 1 / np.sqrt(max(ev[1], 1e-300))
    std2[i] = 1 / np.sqrt(max(ev[0], 1e-300))
b2.plot(s, std1 * 1e9, lw=1.8, color='#1f77b4', label='constrained direction')
b2.plot(s, std2 * 1e9, lw=1.8, color='#d62728', label='degenerate direction')
b2.set_xscale('log'); b2.set_yscale('log')
b2.set_xlabel('data strength $|H|^2/\\sigma^2$', fontsize=10)
b2.set_ylabel('eigen-direction uncertainty [nm]', fontsize=10)
b2.set_title('(b) Posterior precision eigen-directions', fontsize=10)
b2.legend(fontsize=8)
b2.grid(True, which='both', alpha=0.25)
plt.tight_layout()
plt.savefig('manuscript/figI9_joint_bayes.pdf', dpi=150, bbox_inches='tight')
plt.close()

# ------------------------------------------------------------------
# figI7_backus_gilbert.pdf
# ------------------------------------------------------------------
# 2D computation on a well-resolved domain (L = 1024 um) so the
# passband (k_opt = 2/h = 0.05 rad/um) is sampled; spread and variance
# follow Eqs. (bg_spread) and (bg_variance) of the main text.
N, L = 512, 1024e-6
dx = L / N
k = 2 * np.pi * np.fft.fftfreq(N, d=dx)
KX, KY = np.meshgrid(k, k)
KR = np.sqrt(KX**2 + KY**2)
H = C * KX**2 * np.exp(-KR * H_NOM)
Hmax = np.max(np.abs(H))
delta = 2 * np.pi * 5.0
x = np.arange(N) * dx - L / 2
X, Y = np.meshgrid(x, x)
R2 = X**2 + Y**2

def bg_kernel(lam):
    R = H**2 / (H**2 + lam**2)
    A = np.fft.fftshift(np.real(np.fft.ifft2(np.fft.ifftshift(R))))
    return A / (np.sum(A) * dx * dx)

fig, (c1, c2) = plt.subplots(1, 2, figsize=(10.4, 4.4))
for lamr in [0.02, 0.1, 0.5]:
    A = bg_kernel(lamr * Hmax)
    row = A[N // 2]
    row = row / np.max(np.abs(row))
    c1.plot(x * 1e6, row, lw=1.6,
            label=f'$\\lambda={lamr:g}\\,H_{{\\rm max}}$')
c1.set_xlim(-150, 150)
c1.set_xlabel('$x - x_0$ [$\\mu$m]', fontsize=10)
c1.set_ylabel('normalised averaging kernel $A(x,0)$', fontsize=10)
c1.set_title('(a) Averaging kernels (cut through centre)', fontsize=10)
c1.legend(fontsize=8)
c1.grid(True, alpha=0.3)

lams = Hmax * np.logspace(-4, 0, 200)
spreads, varsz = [], []
for lam in lams:
    var = delta**2 * np.sum(H**2 / (H**2 + lam**2)**2)
    varsz.append(var)
    A = bg_kernel(lam)
    spreads.append(np.sqrt(np.sum(R2 * A**2) / np.sum(A**2)) * 1e6)
spreads = np.array(spreads); varsz = np.array(varsz)
sz = np.sqrt(varsz) * 1e9
c2.plot(sz, spreads, 'o-', lw=1.8, color='#1f77b4')
i1 = np.argmin(np.abs(sz - 1.0))
c2.plot(sz[i1], spreads[i1], 's', color='#d62728', ms=9)
c2.annotate(
    f'$\\sigma_z={sz[i1]:.1f}$ nm $\\leftrightarrow\\Delta x_{{\\rm res}}'
    f'\\approx{spreads[i1]:.0f}$ $\\mu$m',
    xy=(sz[i1], spreads[i1]), xytext=(0.42, 0.28),
    textcoords='axes fraction', fontsize=9,
    arrowprops=dict(arrowstyle='->', color='#d62728'))
c2.set_xscale('log')
c2.set_xlabel('height uncertainty $\\sigma_z$ [nm]', fontsize=10)
c2.set_ylabel('resolution spread $\\Delta x_{\\rm res}$ [$\\mu$m]', fontsize=10)
c2.set_title('(b) Resolution--variance trade-off', fontsize=10)
c2.grid(True, which='both', alpha=0.3)
plt.tight_layout()
plt.savefig('manuscript/figI7_backus_gilbert.pdf', dpi=150, bbox_inches='tight')
plt.close()

# ------------------------------------------------------------------
# fig3_amin.pdf
# ------------------------------------------------------------------
fig, (d1, d2) = plt.subplots(1, 2, figsize=(10.4, 4.6))
k = np.logspace(-2, 0, 300) * 1e6        # 0.01 .. 1 rad/um (passband window)
for h_um in [5, 10, 20, 40, 80, 100]:
    h = h_um * 1e-6
    d1.plot(k / 1e6, amin(k, 2 * np.pi * 5, h) * 1e9, lw=1.6,
            label=f'$h={h_um}$ $\\mu$m')
d1.set_xscale('log'); d1.set_yscale('log')
d1.set_xlabel('spatial frequency $k$ [rad/$\\mu$m]', fontsize=10)
d1.set_ylabel('$A_{\\min}(k;h)$ [nm]', fontsize=10)
d1.set_title('(a) Detection threshold ($5$ Hz floor)', fontsize=10)
d1.legend(fontsize=8, ncol=2)
d1.grid(True, which='both', alpha=0.25)

hs = np.logspace(0, 2, 120) * 1e-6
d2.fill_between(hs * 1e6, 64, 1e-3, color='#e8e8e8')
for A_nm in [1, 10, 100, 1000]:
    lmin = []
    for h in hs:
        rhs = (2 * M * W / (E * EBG)) * (2 * np.pi * 5) * h**2 / (A_nm * 1e-9)
        x = np.logspace(-3, 3, 4000)
        y = x**2 * np.exp(-x)
        ilo = np.where(y > rhs)[0]
        xhi = x[ilo[-1]] if len(ilo) else np.nan
        lmin.append(2 * np.pi * h / xhi * 1e6 if xhi else np.nan)
    d2.plot(hs * 1e6, lmin, lw=1.6, label=f'$A={A_nm}$ nm')
d2.axhline(64, color='grey', ls=':', lw=1.4)
d2.text(7, 66, 'FOV $L=64$ $\\mu$m', fontsize=9)
d2.set_xscale('log'); d2.set_yscale('log')
d2.set_ylim(1, 100)
d2.set_xlabel('ion height $h$ [$\\mu$m]', fontsize=10)
d2.set_ylabel('short-wavelength cutoff $\\lambda_{\\min}$ [$\\mu$m]', fontsize=10)
d2.set_title('(b) Reconstructable wavelength band', fontsize=10)
d2.legend(fontsize=8)
d2.grid(True, which='both', alpha=0.25)
plt.tight_layout()
plt.savefig('manuscript/fig3_amin.pdf', dpi=150, bbox_inches='tight')
plt.close()

# ------------------------------------------------------------------
# fig5_uncertainty.pdf
# ------------------------------------------------------------------
fig, e1 = plt.subplots(figsize=(8.6, 5.2))
k = np.logspace(4, 6.3, 400)
z = 100e-9
u_sigma = SIG_CHARGE / (2 * EPS0 * k * EBG)
u_h = k * z * 100e-9
u_ebg = z * 10 / EBG * np.ones_like(k)
u_A = (2 * M * W / (E * EBG)) * (2 * np.pi * 0.3) * np.exp(k * H_NOM) / k**2
u_cal = 1.0e-9 * np.ones_like(k)
uc = np.sqrt(u_sigma**2 + u_h**2 + u_ebg**2 + u_A**2 + u_cal**2)
e1.plot(k, u_sigma * 1e9, lw=2.0, color='#1f77b4', label='surface charge $\\sigma$')
e1.plot(k, u_h * 1e9, lw=1.6, color='#ff7f0e', label='ion height $h$')
e1.plot(k, u_ebg * 1e9, lw=1.4, color='#2ca02c', label='background field $E_{\\rm bg}$')
e1.plot(k, u_A * 1e9, lw=1.4, color='#d62728', label='Type A (frequency measurement)')
e1.plot(k, u_cal * 1e9, lw=1.4, color='#9467bd', label='ITF calibration')
e1.plot(k, uc * 1e9, lw=2.6, color='k', label='combined $u_c(z)$')
e1.axvline(K_BUDGET, color='grey', ls='--', lw=1.2)
e1.text(K_BUDGET, 0.02, '  $k=2\\times10^5$ rad/m\n  $u_c(z)=319$ nm', fontsize=9)
e1.set_xscale('log'); e1.set_yscale('log')
# Show the physically informative band: from the long-wavelength edge of the
# scan to the noise-limited short-wavelength cutoff of the reconstructable
# band (lambda in [19, 64] um at h=40 um), so that every contribution curve
# remains fully visible.
e1.set_xlim(1e4, 4e5)
e1.set_ylim(1e-4, 1e4)
e1.set_xticks([1e4, 2e4, 5e4, 1e5, 2e5, 4e5])
e1.set_xticklabels(['$10^4$', '$2\\times10^4$', '$5\\times10^4$',
                    '$10^5$', '$2\\times10^5$', '$4\\times10^5$'])
e1.axvspan(3.3e5, 4e5, color='grey', alpha=0.10)
e1.text(3.65e5, 2e3, 'beyond reconstructable band', fontsize=8,
        ha='center', va='top', color='0.25', rotation=90)
e1.set_xlabel('spatial frequency $k$ [rad/m]', fontsize=11)
e1.set_ylabel('uncertainty contribution [nm]', fontsize=11)
e1.set_title('GUM-based uncertainty budget vs spatial frequency ($z=100$ nm, $h=40$ $\\mu$m)',
             fontsize=11)
e1.legend(fontsize=8)
e1.grid(True, which='both', alpha=0.25)
plt.tight_layout()
plt.savefig('manuscript/fig5_uncertainty.pdf', dpi=150, bbox_inches='tight')
plt.close()

# ------------------------------------------------------------------
# fig8_micromotion.pdf
# ------------------------------------------------------------------
fig, (f1, f2) = plt.subplots(1, 2, figsize=(10.4, 4.4))
k = np.logspace(-1, 1.6, 400) * 1e6
h = 40e-6
f1.plot(k / 1e6, k**2 * np.exp(-k * h) / (k**2 * np.exp(-k * h)).max(),
        lw=2.0, color='#1f77b4', label='secular $\\Delta\\omega_x$ ($k^2 e^{-kh}$)')
f1.plot(k / 1e6, k * np.exp(-k * h) / (k * np.exp(-k * h)).max(),
        lw=2.0, color='#d62728', label='micromotion $\\beta_x$ ($k e^{-kh}$)')
f1.set_xscale('log')
f1.set_xlabel('spatial frequency $k$ [rad/$\\mu$m]', fontsize=10)
f1.set_ylabel('normalised transfer function', fontsize=10)
f1.set_title('(a) Complementary transfer-function channels ($h=40$ $\\mu$m)',
             fontsize=10)
f1.legend(fontsize=8)
f1.grid(True, which='both', alpha=0.25)

f2.plot(k / 1e6, k / 1e6, lw=2.0, color='#7f7f7f')
f2.set_xscale('log'); f2.set_yscale('log')
f2.set_xlabel('spatial frequency $k$ [rad/$\\mu$m]', fontsize=10)
f2.set_ylabel('$|H_{\\rm sec}|/|H_{\\rm mm}|\\propto k$', fontsize=10)
f2.set_title('(b) Channel ratio: independent height diagnostic', fontsize=10)
f2.grid(True, which='both', alpha=0.25)
plt.tight_layout()
plt.savefig('manuscript/fig8_micromotion.pdf', dpi=150, bbox_inches='tight')
plt.close()

# ------------------------------------------------------------------
# fig9_degeneracy.pdf
# ------------------------------------------------------------------
fig, g1 = plt.subplots(figsize=(8.6, 4.6))
k = np.logspace(4, 6, 400)
h1, h2 = 20e-6, 80e-6
R = np.exp(-k * (h1 - h2))          # = exp(+k*(h2-h1))
g1.plot(k, R, lw=2.2, color='#1f77b4')
g1.set_xscale('log'); g1.set_yscale('log')
g1.set_xlabel('spatial frequency $k$ [rad/m]', fontsize=11)
g1.set_ylabel('$\\Delta\\omega_x(h_1)/\\Delta\\omega_x(h_2)$', fontsize=11)
g1.set_title('Two-height ratio: $e^{-k(h_1-h_2)}$, $h_1=20$ $\\mu$m, $h_2=80$ $\\mu$m',
             fontsize=11)
g1.text(0.55, 0.18,
        'The ratio depends only on $k$ and the heights,\n'
        'not on the $(z_s,\\sigma)$ mixture: the two-height\n'
        'forward matrix is rank 1. Resolution strategies:\n'
        'in-situ charge elimination, a known reference scan,\n'
        'or an independent material-contrast channel (heating rate).',
        transform=g1.transAxes, fontsize=9, va='top',
        bbox=dict(boxstyle='round', facecolor='#fff9c4', alpha=0.9))
g1.grid(True, which='both', alpha=0.25)
plt.tight_layout()
plt.savefig('manuscript/fig9_degeneracy.pdf', dpi=150, bbox_inches='tight')
plt.close()

print('saved the eight orphan figures to manuscript/')
