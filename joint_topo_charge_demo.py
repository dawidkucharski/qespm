#!/usr/bin/env python3
"""Synthetic two-parameter (topography + charge) joint inversion demonstration.

Shows that: (i) single-channel inversion that ignores charge is badly
biased; (ii) joint MAP inversion with a smooth (GP) topography prior and a
sparse (Laplace) charge prior recovers both fields.
Output: manuscript/figV13_joint_tc.pdf
"""
import numpy as np
from numpy.fft import fft, ifft, fftfreq
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt

# ---- parameters (nominal, main text) ----
e = 1.602176634e-19
m = 40 * 1.66053906660e-27
omega_x = 2 * np.pi * 1e6
E_bg = 1e4
C = e * E_bg / (2 * m * omega_x)          # = 1.92e3 m/s
eps0 = 8.8541878128e-12
h = 40e-6
L = 64e-6
N = 256
dx = L / N
x = (np.arange(N) - N // 2) * dx
k = 2 * np.pi * fftfreq(N, dx)

# ---- true fields ----
A_z, w_z = 200e-9, 8e-6
z_true = A_z * np.exp(-x**2 / (2 * w_z**2))
A_s, w_s, x0 = 1.0e-8, 6e-6, 8e-6
sigma_true = A_s * np.exp(-(x - x0)**2 / (2 * w_s**2))

# ---- charge potential (method of images) ----
sig_t = fft(sigma_true)
kap = np.abs(k)
kap[0] = 1e30
phi_charge = np.real(ifft(sig_t / (2 * eps0 * kap)))
phi_charge = phi_charge - phi_charge.mean()
z_charge_eq = phi_charge / E_bg           # equivalent height

# ---- forward model with noise (5 Hz floor) ----
H = -C * k**2 * np.exp(-kap * h)
H[0] = 0.0
d = np.real(ifft(H * fft(z_true + z_charge_eq)))
rng = np.random.default_rng(7)
noise = rng.normal(0.0, 2 * np.pi * 5.0, N)   # 5 Hz floor in rad/s
d = d + noise

# ---- (i) naive single-channel inversion (assumes charge-free) ----
lam = 1e-4 * np.max(np.abs(H))
z_naive = np.real(ifft(np.conj(H) * fft(d) / (np.abs(H)**2 + lam**2)))

def corr(a, b):
    return np.corrcoef(a, b)[0, 1]

r_naive = corr(z_true, z_naive)

# ---- (ii) joint MAP: smooth prior on z, sparse prior on sigma ----
# P: linear map sigma -> equivalent height (phi/E_bg)
def P(sig):
    st = fft(sig)
    out = np.real(ifft(st / (2 * eps0 * kap)))
    return out - out.mean()

def PT(v):
    vt = fft(v)
    out = np.real(ifft(vt / (2 * eps0 * kap)))
    return out - out.mean()

def Hz(v):
    return np.real(ifft(H * fft(v)))

def Hstar(v):
    return np.real(ifft(np.conj(H) * fft(v)))

# ISTA for sigma (fixed z)
def solve_sigma(z, sig0, lam_s, it=400):
    sig = sig0.copy()
    mu_max = C / (2 * eps0) * (1 / h) * np.exp(-1.0)   # sup |H P| symbol
    step = 0.99 / mu_max**2
    r = d - Hz(z)
    for _ in range(it):
        grad = PT(Hstar(Hz(P(sig)) - r))
        sig = sig - step * grad
        sig = np.sign(sig) * np.maximum(np.abs(sig) - step * lam_s, 0.0)
    return sig

# alternating minimisation (closed-form ridge updates in Fourier space)
z_est = z_naive.copy()
sig_est = np.zeros(N)
lam_z = 2.0                                        # topography smoothness weight
lam_s = 1e26                                       # charge prior weight
q = C * kap * np.exp(-kap * h) / (2 * eps0 * E_bg) # symbol of H.P/E_bg
q[0] = 0.0
for it in range(30):
    # z-update (ridge in Fourier)
    denom = np.abs(H)**2 + lam_z * k**4
    denom[0] = 1.0
    z_est = np.real(ifft(np.conj(H) * fft(d - Hz(P(sig_est) / E_bg)) / denom))
    z_est = z_est - z_est.mean()
    # sigma-update (ridge in Fourier): sigma_hat = -q r / (q^2 + lam_s)
    r = d - Hz(z_est)
    sig_est = np.real(ifft(-q * fft(r) / (q**2 + lam_s)))
    sig_est = sig_est - sig_est.mean()

r_joint = corr(z_true, z_est)
r_sigma = corr(sigma_true, sig_est)
rmse_naive = np.sqrt(np.mean((z_naive - z_true)**2))
rmse_joint = np.sqrt(np.mean((z_est - z_true)**2))

print(f"charge-equivalent height peak: {np.max(np.abs(z_charge_eq))*1e9:.0f} nm")
print(f"r(naive z)   = {r_naive:.3f},  RMSE = {rmse_naive*1e9:.1f} nm")
print(f"r(joint z)   = {r_joint:.3f},  RMSE = {rmse_joint*1e9:.1f} nm")
print(f"r(joint sigma) = {r_sigma:.3f}")

# ---- figure ----
xum = x * 1e6
fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
axes[0].plot(xum, z_true * 1e9, 'b-', lw=2, label='$z_s$ (true)')
axes[0].plot(xum, z_charge_eq * 1e9, 'r-', lw=2, label=r'$\phi_{\rm charge}/E_{\rm bg}$')
axes[0].set_xlabel('$x$ [µm]'); axes[0].set_ylabel('[nm]')
axes[0].set_title('(a) True topography and charge-equivalent height')
axes[0].legend(fontsize=8); axes[0].grid(True, alpha=0.3)

axes[1].plot(xum, z_true * 1e9, 'b-', lw=2, label='true $z_s$')
axes[1].plot(xum, z_naive * 1e9, 'k--', lw=1.5, label='naive (charge ignored)')
axes[1].set_xlabel('$x$ [µm]'); axes[1].set_ylabel('[nm]')
axes[1].set_title(f'(b) Single-channel inversion: $r={r_naive:.2f}$')
axes[1].legend(fontsize=8); axes[1].grid(True, alpha=0.3)

axes[2].plot(xum, z_true * 1e9, 'b-', lw=2, label='true $z_s$')
axes[2].plot(xum, z_est * 1e9, 'g--', lw=1.5, label='joint MAP $\\hat z_s$')
ax2 = axes[2].twinx()
ax2.plot(xum, sigma_true * 1e8, 'r-', lw=1.5, label='true $\\sigma$')
ax2.plot(xum, sig_est * 1e8, 'm:', lw=1.5, label='joint MAP $\\hat\\sigma$')
axes[2].set_xlabel('$x$ [µm]'); axes[2].set_ylabel('$z$ [nm]')
ax2.set_ylabel('$\\sigma$ [$10^{-8}$ C/m$^2$]')
axes[2].set_title(f'(c) Joint inversion: $r_z={r_joint:.2f}$, $r_\\sigma={r_sigma:.2f}$')
axes[2].legend(loc='upper left', fontsize=7)
ax2.legend(loc='upper right', fontsize=7)
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('manuscript/figV13_joint_tc.pdf', dpi=150, bbox_inches='tight')
plt.close()
print('Saved: manuscript/figV13_joint_tc.pdf')
