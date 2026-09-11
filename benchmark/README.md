# QESPM Synthetic Benchmark Dataset

Reference data for reproducing the QESPM forward/inverse chain of the
manuscript *"Quantum Limits of Surface Metrology: Spatial Information
Extraction by a Single Trapped Ion"*.

## Contents

- `qespm_benchmark_data.npz` — reference surfaces, simulated ion-signal
  maps (secular-frequency shifts) and reference Tikhonov reconstructions.
- `metadata.json` — parameter values, units, FFT conventions and
  descriptions of every array.
- `generate_benchmark_data.py` — deterministic generator (NumPy only).

## Model

Linear instrument transfer function

    H_x(k; h) = -C k^2 exp(-|k| h),   C = e E_bg / (2 m omega_x) = 1.92e3 m/s,

with E_bg = 1e4 V/m, omega_x = 2π·1 MHz, field of view L = 64 µm,
N = 128 pixels, ion height h = 40 µm. Forward map (numpy FFT convention):

    Δω̃_x(k) = H_x(k; h) z̃(k),    Δω_x(x) = ifft(H_x z̃).

Reference reconstructions use Tikhonov regularisation
ẑ = ifft( conj(H)/(|H|² + λ) · fft(d) ), with λ = 1.0 (noise-free)
and λ = 5e5 (5 Hz white-noise floor data).

## Reproduce

    python generate_benchmark_data.py   # regenerates the .npz + metadata.json

    import numpy as np
    d = np.load("qespm_benchmark_data.npz")
    # keys: x, kx, h, L, N, dx, C_ITF, noise_hz, z_ref,
    #       z_grating, z_fine, z_plateau, z_rough,
    #       dw_clean, dw_noisy, z_rec_clean, z_rec_noisy

## Licence

MIT — same as the repository.
