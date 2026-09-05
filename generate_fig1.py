#!/usr/bin/env python3
"""Generate Figure 1 only — direct numerical integration of ion trajectory."""
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy.fft import fft, fftfreq

E_CHARGE = 1.602176634e-19
M_CA40   = 40.0 * 1.660539e-27


class TrapConfig:
    def __init__(self, freq_sec=1.0e6, ion_mass=M_CA40, ion_height=40e-6,
                 E_bg=1.0e4):
        self.freq_sec = freq_sec
        self.omega_sec = 2.0 * np.pi * freq_sec
        self.ion_mass = ion_mass
        self.ion_height = ion_height
        self.E_bg = E_bg
        self.C = E_CHARGE * E_bg / (2.0 * ion_mass * self.omega_sec)


def surface_force(x, A, k_s, h, E_bg):
    return -E_CHARGE * E_bg * A * k_s * np.cos(k_s * x) * np.exp(-k_s * h)


def secular_ode(t, state, trap, A, k_s, x_scan):
    x, v = state
    a_trap = -trap.omega_sec**2 * (x - x_scan)
    a_surf = surface_force(x, A, k_s, trap.ion_height, trap.E_bg) / trap.ion_mass
    return [v, a_trap + a_surf]


def integrate_and_measure(trap, A, k_s, x_scan, t_max=2e-3, n_pts=50000,
                          x0_offset=50e-9):
    t_eval = np.linspace(0, t_max, n_pts)
    sol = solve_ivp(secular_ode, (0, t_max), [x_scan + x0_offset, 0.0],
                    args=(trap, A, k_s, x_scan), t_eval=t_eval,
                    method='RK45', rtol=1e-4, atol=1e-9)
    t = sol.t
    x = sol.y[0]
    n_half = len(t) // 2
    x_ss = x[n_half:]
    t_ss = t[n_half:]
    x_centred = x_ss - np.mean(x_ss)
    signs = np.sign(x_centred)
    zc = []
    for i in range(1, len(signs)):
        if signs[i-1] < 0 and signs[i] >= 0:
            frac = -x_centred[i-1] / (x_centred[i] - x_centred[i-1] + 1e-30)
            zc.append(t_ss[i-1] + frac * (t_ss[i] - t_ss[i-1]))
    if len(zc) < 3:
        return np.nan, np.nan, t, x
    periods = np.diff(zc)
    T_mean = np.mean(periods)
    T_std = np.std(periods) / np.sqrt(len(periods))
    return 1.0 / T_mean, (1.0 / T_mean) * T_std / T_mean, t, x


def analytical_dw(trap, A, k_s, x_pos):
    return -trap.C * A * k_s**2 * np.sin(k_s * x_pos) * np.exp(
        -k_s * trap.ion_height)


if __name__ == "__main__":
    plt.rcParams.update({
        'figure.dpi': 150, 'savefig.dpi': 150, 'savefig.bbox': 'tight',
        'font.family': 'serif', 'font.size': 11, 'mathtext.fontset': 'stix'
    })

    trap = TrapConfig()
    A_test, lam_test = 200e-9, 40e-6
    k_test = 2.0 * np.pi / lam_test
    x_demo = lam_test / 4

    print('Integrating with surface (5 ms, 200k pts, rtol=1e-4)...')
    f_demo, df_demo, t_demo, x_demo_arr = integrate_and_measure(
        trap, A_test, k_test, x_demo, t_max=5e-3, n_pts=200000,
        x0_offset=100e-9)
    print(f'  f_meas = {f_demo*1e-6:.6f} MHz, df = {df_demo:.1f} Hz')

    print('Integrating without surface...')
    f_ref, _, t_ref, x_ref = integrate_and_measure(
        trap, 0.0, k_test, x_demo, t_max=5e-3, n_pts=200000,
        x0_offset=100e-9)
    print(f'  f_ref  = {f_ref*1e-6:.6f} MHz')

    dw_demo = 2 * np.pi * (f_demo - trap.freq_sec)
    dw_ana_demo = analytical_dw(trap, A_test, k_test, x_demo)
    print(f'  dw/2pi (sim) = {dw_demo/(2*np.pi):.1f} Hz')
    print(f'  dw/2pi (ITF) = {dw_ana_demo/(2*np.pi):.1f} Hz')
    print(f'  Agreement: {abs(dw_demo-dw_ana_demo)/abs(dw_ana_demo)*100:.1f}%')

    fig1, ((ax_a, ax_b), (ax_c, ax_d)) = plt.subplots(2, 2, figsize=(16, 12))
    t_us = t_demo * 1e6

    # (a) Full 5 ms trajectory
    ax_a.plot(t_us, x_demo_arr * 1e9, 'b-', lw=0.3, alpha=0.8)
    ax_a.set_xlabel('Time [µs]', fontsize=12)
    ax_a.set_ylabel('$x$ position [nm]', fontsize=12)
    ax_a.set_title(
        r'(a) Ion trajectory at $x_{\rm scan} = \lambda/4$ (max shift)',
        fontsize=12, fontweight='bold')
    ax_a.grid(True, alpha=0.3)

    # (b) Zoom: 15 periods
    pts_per_period = int(1.0 / (trap.freq_sec * (t_us[1] - t_us[0]) * 1e-6))
    n_periods_zoom, start_idx = 15, 8000
    n_zoom = n_periods_zoom * pts_per_period
    t_zoom = t_us[start_idx:start_idx + n_zoom]
    ax_b.plot(t_zoom, x_demo_arr[start_idx:start_idx + n_zoom] * 1e9,
              '#2166AC', lw=2.0, label='With surface ($A=200$ nm)')
    ax_b.plot(t_zoom, x_ref[start_idx:start_idx + n_zoom] * 1e9,
              '#D6604D', lw=2.0, ls='--', label='Without surface ($A=0$)')
    ax_b.set_xlabel('Time [µs]', fontsize=12)
    ax_b.set_ylabel('$x$ position [nm]', fontsize=12)
    ax_b.set_title('(b) Zoom: 15 secular periods', fontsize=12,
                   fontweight='bold')
    ax_b.legend(fontsize=9, loc='upper right', framealpha=0.9)
    ax_b.grid(True, alpha=0.3)

    # (c) Full PSD
    n_fft = len(t_demo) // 2
    x_fft = x_demo_arr[-n_fft:]
    X = np.abs(fft(x_fft * np.hanning(n_fft)))**2
    f_fft = fftfreq(n_fft, t_demo[1] - t_demo[0])
    pos = f_fft > 0
    f_ita = trap.freq_sec + dw_ana_demo / (2 * np.pi)
    ax_c.semilogy(f_fft[pos] * 1e-6, X[pos] / X[pos].max(), 'b-', lw=1)
    ax_c.axvline(f_ita * 1e-6, color='#D6604D', ls='--', lw=2.5,
                 label=f'ITF: {f_ita*1e-6:.6f} MHz')
    ax_c.axvline(trap.freq_sec * 1e-6, color='gray', ls=':', lw=2.0,
                 label=f'Unperturbed: {trap.freq_sec*1e-6:.4f} MHz')
    ax_c.set_xlim(0.8, 1.2)
    ax_c.set_xlabel('Frequency [MHz]', fontsize=12)
    ax_c.set_ylabel('Normalised PSD', fontsize=12)
    ax_c.set_title('(c) Power spectral density', fontsize=12, fontweight='bold')
    ax_c.legend(fontsize=9, loc='upper right')
    ax_c.grid(True, alpha=0.3)

    # (d) Zoom on secular peak
    mask_d = (f_fft > 0.994e6) & (f_fft < 1.002e6)
    ax_d.plot(f_fft[mask_d] * 1e-6, X[mask_d] / X[mask_d].max(), 'b-', lw=1.8)
    ax_d.axvline(f_ita * 1e-6, color='#D6604D', ls='--', lw=2.5,
                 label='ITF prediction')
    ax_d.axvline(trap.freq_sec * 1e-6, color='gray', ls=':', lw=2.0,
                 label='Unperturbed')
    ax_d.annotate(
        f"ITF: $|\\Delta f| = {abs(dw_ana_demo/(2*np.pi)):.0f}$ Hz\n"
        f"$(\\Delta\\omega/\\omega = "
        f"{abs(dw_ana_demo)/trap.omega_sec*100:.3f}\\%)$",
        xy=(f_ita * 1e-6, 0.90), fontsize=11, color='#B2182B',
        fontweight='bold',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.90,
                  edgecolor='#D6604D', lw=1.5))
    sim_note = (f"Sim. zero-crossing: $\\Delta f = "
                f"{dw_demo/(2*np.pi):.0f}$ Hz "
                f"(FFT res. {1/(t_demo[-1]-t_demo[0]):.0f} Hz)")
    ax_d.annotate(sim_note, xy=(0.02, 0.10), xycoords='axes fraction',
                  fontsize=8, color='gray', fontstyle='italic')
    ax_d.set_xlabel('Frequency [MHz]', fontsize=12)
    ax_d.set_ylabel('Normalised PSD', fontsize=12)
    ax_d.set_title('(d) Secular peak — ITF prediction vs. unperturbed',
                   fontsize=12, fontweight='bold')
    ax_d.legend(fontsize=8, loc='upper left', framealpha=0.85)
    ax_d.grid(True, alpha=0.3)

    fig1.suptitle(
        'Direct Numerical Integration: Ion Dynamics Near a Structured Surface\n'
        r'$^{40}$Ca$^+$, $f_{\rm sec}=1$ MHz, $h=40$ $\mu$m, '
        r'$A=200$ nm, $\lambda=40$ $\mu$m, '
        r'$x_{\rm scan}=\lambda/4$, $5$ ms integration',
        fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig1.savefig('figS1_trajectory.pdf')
    print('Saved figS1_trajectory.pdf')
