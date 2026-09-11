#!/usr/bin/env python3
"""QESPM synthetic benchmark dataset.

Generates reference surfaces, simulated ion-signal maps (linear ITF) and
reference Tikhonov reconstructions exactly as used in the manuscript
(Sec. "Validation on synthetic surface textures"):

  L = 64 um field of view, N = 128 pixels, h = 40 um,
  C = e E_bg / (2 m omega_x) = 1.92e3 m/s  (E_bg = 1e4 V/m, omega_x = 2pi*1 MHz),
  H_x(k; h) = -C k^2 exp(-|k| h).

Signals are secular-frequency shifts Delta_omega_x in rad/s.
FFT convention: numpy (unnormalised) FFT of the centred arrays.

Outputs:
  benchmark/qespm_benchmark_data.npz
  benchmark/metadata.json
"""
import json
import numpy as np
from numpy.fft import fft, fftfreq, ifft

rng = np.random.default_rng(20260911)

# ---------- nominal QESPM parameters (manuscript) ----------
L = 64e-6            # field of view [m]
N = 128              # pixels
DX = L / N           # 0.5 um
h = 40e-6            # ion height [m]
C_ITF = 1.92e3       # m/s
NOISE_HZ = 5.0       # conservative technical floor [Hz]

x = (np.arange(N) - N // 2) * DX           # centred coordinate
kx = 2 * np.pi * fftfreq(N, d=DX)          # angular wavenumbers


def itf(k):
    return -C_ITF * k ** 2 * np.exp(-np.abs(k) * h)


# ---------- reference surfaces ----------
def grating(lam, A):
    return A * np.sin(2 * np.pi * x / lam)


# resolvable coarse grating (lambda = 100 um > lambda_min = 1.51h = 60.4 um)
z_grating = grating(100e-6, 100e-9)
# fine grating below the MTF cutoff (lambda = 4 um, tests regularisation)
z_fine = grating(4e-6, 50e-9)
# Gaussian plateau (resolvable)
z_plateau = 200e-9 * np.exp(-0.5 * (x / 30e-6) ** 2)
# broadband power-law roughness, normalised to Ra = 100 nm
lc, beta = 5e-6, 1.5
S = np.where(kx != 0, (1 + (kx * lc) ** 2) ** (-beta), 0.0)
z_rough = np.fft.ifft(np.sqrt(S) * np.exp(2j * np.pi * rng.random(N))).real
z_rough -= z_rough.mean()
z_rough *= 100e-9 / (np.mean(np.abs(z_rough)) + 1e-30)

z_ref = z_grating + z_plateau + z_rough + z_fine


# ---------- forward model and reference reconstruction ----------
def forward(z, sigma_noise=0.0):
    dw = ifft(itf(kx) * fft(z)).real
    if sigma_noise > 0:
        dw += rng.normal(0.0, sigma_noise, N)
    return dw


def tikhonov(dw, lam):
    H = itf(kx)
    zf = np.conj(H) / (np.abs(H) ** 2 + lam) * fft(dw)
    return ifft(zf).real


dw_clean = forward(z_ref)
dw_noisy = forward(z_ref, sigma_noise=2 * np.pi * NOISE_HZ)

LAM_CLEAN = 1.0     # Tikhonov parameter (documented; small since noise-free)
LAM_NOISY = 5e5     # Tikhonov parameter for the 5 Hz-floor data
z_rec_clean = tikhonov(dw_clean, LAM_CLEAN)
z_rec_noisy = tikhonov(dw_noisy, LAM_NOISY)

# ---------- output ----------
np.savez(
    "benchmark/qespm_benchmark_data.npz",
    x=x, kx=kx, h=h, L=L, N=N, dx=DX, C_ITF=C_ITF, noise_hz=NOISE_HZ,
    z_ref=z_ref,
    z_grating=z_grating, z_fine=z_fine, z_plateau=z_plateau, z_rough=z_rough,
    dw_clean=dw_clean, dw_noisy=dw_noisy,
    z_rec_clean=z_rec_clean, z_rec_noisy=z_rec_noisy,
)

metadata = {
    "description": (
        "QESPM synthetic benchmark: reference surfaces, simulated ion-signal "
        "maps and reference reconstructions for the linear instrument transfer "
        "function H_x(k;h) = -C k^2 exp(-|k| h), C = 1.92e3 m/s."
    ),
    "units": {
        "x": "m (centred, dx = 0.5 um)", "kx": "rad/m (numpy fftfreq order)",
        "z_*": "m", "dw_*": "rad/s (secular frequency shift Delta omega_x)",
    },
    "parameters": {
        "h": h, "L": L, "N": N, "C_ITF": C_ITF,
        "noise_hz": NOISE_HZ,
        "E_bg": 1e4, "omega_x_rad_s": 2 * np.pi * 1e6,
    },
    "surfaces": {
        "z_grating": "sinusoid, lambda = 100 um, A = 100 nm (resolvable)",
        "z_plateau": "Gaussian plateau, sigma = 30 um, A = 200 nm (resolvable)",
        "z_rough": "power-law roughness, lc = 5 um, beta = 1.5, Ra = 100 nm",
        "z_fine": "sinusoid, lambda = 4 um, A = 50 nm (below MTF cutoff 1.51h = 60.4 um)",
        "z_ref": "sum of the four components",
    },
    "reference_reconstructions": {
        "z_rec_clean": "Tikhonov, lambda = 1.0 (noise-free data)",
        "z_rec_noisy": "Tikhonov, lambda = 5e5 (data with 5 Hz white-noise floor)",
    },
    "conventions": (
        "numpy unnormalised FFT; dw = ifft(H(kx) fft(z)).real; "
        "z = ifft(conj(H)/(|H|^2 + lambda) fft(dw)).real; "
        "phase conventions follow the centred-coordinate FFT used in the manuscript."
    ),
}

with open("benchmark/metadata.json", "w") as f:
    json.dump(metadata, f, indent=2)

print("saved benchmark/qespm_benchmark_data.npz and benchmark/metadata.json")
print(f"  z_ref: mean={z_ref.mean():.3e} m, Ra={np.abs(z_ref - z_ref.mean()).mean()*1e9:.1f} nm")
print(f"  dw_clean: max|dw|/2pi = {np.abs(dw_clean).max()/(2*np.pi):.1f} Hz")
print(f"  dw_noisy: max|dw|/2pi = {np.abs(dw_noisy).max()/(2*np.pi):.1f} Hz")
print(f"  corr(z_rec_clean, z_ref) = {np.corrcoef(z_rec_clean, z_ref)[0, 1]:.4f}")
print(f"  corr(z_rec_noisy, z_ref) = {np.corrcoef(z_rec_noisy, z_ref)[0, 1]:.4f}")
