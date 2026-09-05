#!/usr/bin/env python3
"""
paul_trap_simulation.py — Ion Dynamics + Surface Sensing
=========================================================
Direct numerical integration of the secular (pseudopotential) equation
of motion for a trapped ion near a sinusoidal surface.

Validates the analytical ITF H(k) = -(e*E_bg/2m*w)*k^2*exp(-k*h) by:
  1. Integrating the equation of motion at several scan positions
  2. Extracting the secular frequency via zero-crossing period measurement
  3. Comparing simulated frequency shift with analytical prediction

Physical model:
  d2x/dt2 + w_sec^2 * x = -(e/m) * dPhi_surface/dx
  Phi_surface(x) = E_bg * A * sin(k_s*x) * exp(-k_s*h)

Author: QESPM theory paper
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
import matplotlib.pyplot as plt

E_CHARGE = 1.602176634e-19
M_CA40   = 40.0 * 1.660539e-27


class TrapConfig:
    def __init__(self, freq_sec=1.0e6, ion_mass=M_CA40,
                 ion_height=40e-6, E_bg=1.0e4):
        self.freq_sec = freq_sec
        self.omega_sec = 2.0 * np.pi * freq_sec
        self.ion_mass = ion_mass
        self.ion_height = ion_height
        self.E_bg = E_bg
        # ITF prefactor
        self.C = E_CHARGE * E_bg / (2.0 * ion_mass * self.omega_sec)


def surface_force(x, A, k_s, h, E_bg):
    """F_x = -e * dPhi/dx for Phi = E_bg*A*sin(k_s*x)*exp(-k_s*h)."""
    return -E_CHARGE * E_bg * A * k_s * np.cos(k_s * x) * np.exp(-k_s * h)


def secular_ode(t, state, trap, A, k_s, x_scan=0.0):
    """Secular EOM: dx/dt=v, dv/dt = -w^2*(x-x_scan) + F_surface/m."""
    x, v = state
    a_trap = -trap.omega_sec**2 * (x - x_scan)
    a_surf = surface_force(x, A, k_s, trap.ion_height, trap.E_bg) / trap.ion_mass
    return [v, a_trap + a_surf]


def integrate_and_measure(trap, A, k_s, x_scan, t_max=2e-3, n_pts=50000,
                          x0_offset=50e-9):
    """
    Integrate secular EOM and extract frequency via zero-crossing period.

    The ion is trapped with equilibrium at x_scan (trap centre position).
    A small initial offset x0_offset starts the secular oscillation.
    Uses the last 50% of the trajectory for steady-state measurement.
    """
    t_eval = np.linspace(0, t_max, n_pts)
    sol = solve_ivp(secular_ode, (0, t_max), [x_scan + x0_offset, 0.0],
                    args=(trap, A, k_s, x_scan), t_eval=t_eval,
                    method="RK45", rtol=1e-8, atol=1e-11)

    t = sol.t
    x = sol.y[0]

    n_half = len(t) // 2
    x_ss = x[n_half:]
    t_ss = t[n_half:]

    # Detect zero crossings around the mean position (not x=0)
    x_centred = x_ss - np.mean(x_ss)
    signs = np.sign(x_centred)
    zero_crossings = []
    for i in range(1, len(signs)):
        if signs[i-1] < 0 and signs[i] >= 0:
            frac = -x_centred[i-1] / (x_centred[i] - x_centred[i-1] + 1e-30)
            t_cross = t_ss[i-1] + frac * (t_ss[i] - t_ss[i-1])
            zero_crossings.append(t_cross)

    if len(zero_crossings) < 3:
        return np.nan, np.nan, t, x

    periods = np.diff(zero_crossings)
    T_mean = np.mean(periods)
    T_std = np.std(periods) / np.sqrt(len(periods))

    f_meas = 1.0 / T_mean
    df = f_meas * T_std / T_mean

    return f_meas, df, t, x


def analytical_dw(trap, A, k_s, x_pos):
    """Analytical frequency shift from ITF."""
    return -trap.C * A * k_s**2 * np.sin(k_s * x_pos) * np.exp(-k_s * trap.ion_height)


# =====================================================================
# MAIN
# =====================================================================
if __name__ == "__main__":
    matplotlib.use("Agg")
    plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 150,
                         "savefig.bbox": "tight", "font.family": "serif",
                         "font.size": 11, "mathtext.fontset": "stix"})

    trap = TrapConfig()
    A_test = 200e-9
    lam_test = 40e-6
    k_test = 2.0 * np.pi / lam_test

    print("=" * 60)
    print("  Ion Dynamics + Surface Sensing — Direct Integration")
    print(f"  f_sec = {trap.freq_sec/1e6:.1f} MHz  h = {trap.ion_height*1e6:.0f} um")
    print(f"  A = {A_test*1e9:.0f} nm  lambda = {lam_test*1e6:.0f} um")
    print(f"  ITF prefactor C = {trap.C:.2f}")
    print("=" * 60)

    # ---- Test at 5 positions ----
    x_test_positions = np.array([0, lam_test/8, lam_test/4,
                                  3*lam_test/8, lam_test/2])
    f_sim = np.zeros(len(x_test_positions))
    df_sim = np.zeros(len(x_test_positions))
    dw_ana = np.zeros(len(x_test_positions))

    for i, xp in enumerate(x_test_positions):
        f, df, t, x = integrate_and_measure(trap, A_test, k_test, xp,
                                             t_max=3e-3, n_pts=120000)
        f_sim[i] = f
        df_sim[i] = df
        dw_ana[i] = analytical_dw(trap, A_test, k_test, xp)
        dw_s = 2*np.pi*(f - trap.freq_sec) if np.isfinite(f) else np.nan
        print(f"  x = {xp*1e6:6.1f} um:  f = {f*1e-6:.6f} MHz,  "
              f"dw/2pi(sim) = {dw_s/(2*np.pi):7.1f} Hz,  "
              f"dw/2pi(ana) = {dw_ana[i]/(2*np.pi):7.1f} Hz")

    dw_sim_arr = 2.0 * np.pi * (f_sim - trap.freq_sec)
    valid = np.isfinite(dw_sim_arr)
    if np.sum(valid) >= 3:
        corr = np.corrcoef(dw_sim_arr[valid], dw_ana[valid])[0, 1]
        print(f"  Correlation: r = {corr:.4f}")

    # ================================================================
    # Figure S1: Single trajectory at x = lambda/4 (max shift)
    # ================================================================
    x_demo = lam_test / 4
    f_demo, df_demo, t_demo, x_demo_arr = integrate_and_measure(
        trap, A_test, k_test, x_demo, t_max=5e-3, n_pts=200000,
        x0_offset=100e-9)

    # Also integrate without surface (A=0) for reference
    f_ref, _, t_ref, x_ref = integrate_and_measure(
        trap, 0.0, k_test, x_demo, t_max=5e-3, n_pts=200000,
        x0_offset=100e-9)

    dw_demo = 2*np.pi*(f_demo - trap.freq_sec)
    dw_ana_demo = analytical_dw(trap, A_test, k_test, x_demo)
    print(f"\n  Demo at x=lambda/4:")
    print(f"    dw/2pi (sim, zero-crossing) = {dw_demo/(2*np.pi):.1f} Hz")
    print(f"    dw/2pi (analytical ITF)     = {dw_ana_demo/(2*np.pi):.1f} Hz")
    print(f"    Agreement: {abs(dw_demo-dw_ana_demo)/abs(dw_ana_demo)*100:.1f}%")

    fig1, ((ax_a, ax_b), (ax_c, ax_d)) = plt.subplots(2, 2, figsize=(16, 12))

    # (a) Full trajectory with surface — show first 5 ms of 10 ms
    t_us = t_demo * 1e6
    mask_a = t_us <= 5000
    ax_a.plot(t_us[mask_a], x_demo_arr[mask_a] * 1e9, "b-", lw=0.3, alpha=0.8)
    ax_a.set_xlabel("Time [µs]", fontsize=12)
    ax_a.set_ylabel("$x$ position [nm]", fontsize=12)
    ax_a.set_title("(a) Ion trajectory at $x_{\\rm scan} = \\lambda/4$ (max shift)",
                   fontsize=12, fontweight="bold")
    ax_a.grid(True, alpha=0.3)

    # (b) Zoom: 15 periods with clear visual separation
    n_periods_zoom = 15
    pts_per_period = int(1.0 / (trap.freq_sec * (t_us[1]-t_us[0]) * 1e-6))
    n_zoom = n_periods_zoom * pts_per_period
    start_idx = 8000
    t_zoom = t_us[start_idx:start_idx+n_zoom]
    x_zoom_surf = x_demo_arr[start_idx:start_idx+n_zoom] * 1e9
    x_zoom_ref  = x_ref[start_idx:start_idx+n_zoom] * 1e9
    ax_b.plot(t_zoom, x_zoom_surf, "#2166AC", lw=2.0, label="With surface ($A=200$ nm)")
    ax_b.plot(t_zoom, x_zoom_ref,  "#D6604D", lw=2.0, ls="--",
              label="Without surface ($A=0$)")
    ax_b.set_xlabel("Time [µs]", fontsize=12)
    ax_b.set_ylabel("$x$ position [nm]", fontsize=12)
    ax_b.set_title("(b) Zoom: 15 secular periods",
                   fontsize=12, fontweight="bold")
    ax_b.legend(fontsize=9, loc="upper right", framealpha=0.9)
    ax_b.grid(True, alpha=0.3)

    # (c) FFT comparison — full PSD
    from numpy.fft import fft, fftfreq
    n_fft = len(t_demo) // 2
    x_fft = x_demo_arr[-n_fft:]
    X = np.abs(fft(x_fft * np.hanning(n_fft)))**2
    f_fft = fftfreq(n_fft, t_demo[1]-t_demo[0])
    pos = f_fft > 0
    f_ita = trap.freq_sec + dw_ana_demo/(2*np.pi)  # ITF-predicted frequency
    ax_c.semilogy(f_fft[pos]*1e-6, X[pos]/X[pos].max(), "b-", lw=1)
    ax_c.axvline(f_ita*1e-6, color="#D6604D", ls="--", lw=2.5,
                 label=f"ITF prediction: {f_ita*1e-6:.6f} MHz")
    ax_c.axvline(trap.freq_sec*1e-6, color="gray", ls=":", lw=2.0,
                 label=f"Unperturbed: {trap.freq_sec*1e-6:.4f} MHz")
    ax_c.set_xlim(0.8, 1.2)
    ax_c.set_xlabel("Frequency [MHz]", fontsize=12)
    ax_c.set_ylabel("Normalised PSD", fontsize=12)
    ax_c.set_title("(c) Power spectral density", fontsize=12, fontweight="bold")
    ax_c.legend(fontsize=9, loc="upper right")
    ax_c.grid(True, alpha=0.3)

    # (d) Zoom on secular peak — show both ITF prediction and measured shift
    mask_d = (f_fft > 0.99e6) & (f_fft < 1.01e6)
    ax_d.plot(f_fft[mask_d]*1e-6, X[mask_d]/X[mask_d].max(), "b-", lw=1.8)
    # ITF-predicted shifted peak
    ax_d.axvline(f_ita*1e-6, color="#D6604D", ls="--", lw=2.5,
                 label=f"ITF: $f_{{\\rm sec}}+\\Delta f$ = {f_ita*1e-6:.6f} MHz")
    # Unperturbed reference
    ax_d.axvline(trap.freq_sec*1e-6, color="gray", ls=":", lw=2.0,
                 label=f"$f_{{\\rm sec}}$ = {trap.freq_sec*1e-6:.4f} MHz")
    # Annotation with ITF prediction
    ax_d.annotate(
        f"ITF prediction:\n"
        f"$|\\Delta f| = {abs(dw_ana_demo/(2*np.pi)):.0f}$ Hz\n"
        f"$(\\Delta\\omega/\\omega = {abs(dw_ana_demo)/trap.omega_sec*100:.3f}\\%)$",
        xy=(f_ita*1e-6, 0.90), fontsize=11, color="#B2182B", fontweight="bold",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.90,
                  edgecolor="#D6604D", lw=1.5)
    )
    # Small note about simulation resolution
    ax_d.annotate(
        f"(Sim. $\\Delta f$ = {dw_demo/(2*np.pi):.0f} Hz;\n"
        f"FFT res. = {1/(t_demo[-1]-t_demo[0])*2:.0f} Hz)",
        xy=(0.02, 0.10), xycoords="axes fraction",
        fontsize=8, color="gray", fontstyle="italic"
    )
    ax_d.set_xlabel("Frequency [MHz]", fontsize=12)
    ax_d.set_ylabel("Normalised PSD", fontsize=12)
    ax_d.set_title("(d) Secular peak — ITF prediction vs.~unperturbed",
                   fontsize=12, fontweight="bold")
    ax_d.legend(fontsize=8, loc="upper left", framealpha=0.85)
    ax_d.grid(True, alpha=0.3)

    fig1.suptitle(
        "Direct Numerical Integration: Ion Dynamics Near a Structured Surface\n"
        "$^{40}$Ca$^+$, $f_{\\rm sec}=1$ MHz, $h=40$ µm, "
        "$A=200$ nm, $\\lambda=40$ µm, $x_{\\rm scan}=\\lambda/4$, "
        "integration time $=5$ ms",
        fontsize=14, fontweight="bold",
    )
    plt.tight_layout()
    fig1.savefig("figS1_trajectory.pdf")
    print("  Saved figS1_trajectory.pdf")

    # ================================================================
    # Figure S2: Scan validation — simulation points + analytical curve
    # ================================================================
    fig2, (ax2a, ax2b) = plt.subplots(1, 2, figsize=(14, 6))

    # Analytical scan (dense)
    x_dense = np.linspace(0, lam_test, 200)
    dw_dense = analytical_dw(trap, A_test, k_test, x_dense)

    # Simulation points
    ax2a.plot(x_dense * 1e6, dw_dense / (2*np.pi), "r-", lw=2,
              label="Analytical ITF")
    ax2a.errorbar(x_test_positions * 1e6, dw_sim_arr / (2*np.pi),
                   yerr=df_sim / (2*np.pi), fmt="bo", ms=10, capsize=6,
                   label="Simulation (zero-crossing)")
    ax2a.set_xlabel("Scan position $x$ [µm]", fontsize=12)
    ax2a.set_ylabel("$\\Delta\\omega_x / 2\\pi$ [Hz]", fontsize=12)
    ax2a.set_title("(a) Frequency shift vs.~scan position", fontsize=12,
                   fontweight="bold")
    ax2a.legend(fontsize=10)
    ax2a.grid(True, alpha=0.3)

    # Scatter: simulation vs analytical
    ax2b.errorbar(dw_ana[valid] / (2*np.pi), dw_sim_arr[valid] / (2*np.pi),
                   xerr=0, yerr=df_sim[valid] / (2*np.pi),
                   fmt="ko", ms=8, capsize=4)
    lim = np.max(np.abs(dw_dense)) / (2*np.pi) * 1.3
    ax2b.plot([-lim, lim], [-lim, lim], "r--", lw=1.5)
    ax2b.set_xlabel("$\\Delta\\omega_x^{\\rm ITF} / 2\\pi$ [Hz]", fontsize=12)
    ax2b.set_ylabel("$\\Delta\\omega_x^{\\rm sim} / 2\\pi$ [Hz]", fontsize=12)
    if np.sum(valid) >= 3:
        ax2b.set_title(f"(b) Simulation vs.~analytical  "
                       f"($r = {corr:.4f}$)", fontsize=12, fontweight="bold")
    ax2b.set_aspect("equal")
    ax2b.grid(True, alpha=0.3)

    fig2.suptitle(
        "Validation of the Analytical ITF by Direct Numerical Integration\n"
        f"$h={trap.ion_height*1e6:.0f}$ µm, $A={A_test*1e9:.0f}$ nm, "
        f"$\\lambda={lam_test*1e6:.0f}$ µm",
        fontsize=13, fontweight="bold",
    )
    plt.tight_layout()
    fig2.savefig("figS2_scan_validation.pdf")
    print("  Saved figS2_scan_validation.pdf")

    # ================================================================
    # Figure S3: Frequency shift vs surface amplitude (linearity)
    # ================================================================
    A_vals = np.array([10, 25, 50, 100, 200, 400, 800]) * 1e-9
    x_lin = lam_test / 4
    dw_lin_sim = np.zeros(len(A_vals))
    df_lin = np.zeros(len(A_vals))

    for i, Ai in enumerate(A_vals):
        f_i, df_i, _, _ = integrate_and_measure(
            trap, Ai, k_test, x_lin, t_max=3e-3, n_pts=120000)
        dw_lin_sim[i] = 2*np.pi*(f_i - trap.freq_sec)
        df_lin[i] = df_i

    dw_lin_ana = analytical_dw(trap, A_vals, k_test, x_lin)

    fig3, ax3 = plt.subplots(1, 1, figsize=(8, 6))
    ax3.loglog(A_vals * 1e9, np.abs(dw_lin_sim) / (2*np.pi),
               "bo-", ms=8, lw=2, label="Simulation (zero-crossing)")
    ax3.loglog(A_vals * 1e9, np.abs(dw_lin_ana) / (2*np.pi),
               "r--", lw=2, label="$|\\Delta\\omega_x| \\propto A$ (linear ITF)")
    ax3.set_xlabel("Surface amplitude $A$ [nm]", fontsize=12)
    ax3.set_ylabel("$|\\Delta\\omega_x| / 2\\pi$ [Hz]", fontsize=12)
    ax3.set_title("Linearity of Surface Perturbation\n"
                  f"$h={trap.ion_height*1e6:.0f}$ µm, "
                  f"$x=\\lambda/4$, $\\lambda={lam_test*1e6:.0f}$ µm",
                  fontsize=12, fontweight="bold")
    ax3.legend(fontsize=10)
    ax3.grid(True, alpha=0.3, which="both")
    # Validity boundary
    ax3.axvline(trap.ion_height * 0.1 * 1e9, color="gray", ls=":", alpha=0.5, lw=1.5)
    ax3.text(trap.ion_height * 0.1 * 1e9 * 1.15, 0.5,
             "$A = h/10$\n(perturbation\nvalidity limit)",
             fontsize=9, color="gray")
    fig3.tight_layout()
    fig3.savefig("figS3_linearity.pdf")
    print("  Saved figS3_linearity.pdf")

    # ================================================================
    # Figure S4: Comprehensive signal map (h, lambda) from analytical ITF
    # ================================================================
    h_vals = np.logspace(-5.5, -3.5, 40)  # 3 µm → 300 µm
    lam_vals = np.logspace(-5.5, -3.5, 40)
    HH, LL = np.meshgrid(h_vals, lam_vals)
    KK = 2.0 * np.pi / LL
    A_map = 100e-9
    dw_map = trap.C * A_map * KK**2 * np.exp(-KK * HH)
    dw_hz = dw_map / (2.0 * np.pi)

    fig4, ax4 = plt.subplots(1, 1, figsize=(9, 7))
    cs = ax4.contourf(LL * 1e6, HH * 1e6, np.log10(dw_hz + 1e-30),
                       levels=20, cmap="viridis")
    # SNR=1 contour (noise floor 5 Hz)
    ct = ax4.contour(LL * 1e6, HH * 1e6, dw_hz,
                      levels=[5, 50, 500, 5000, 50000],
                      colors=["white", "yellow", "orange", "red", "magenta"],
                      linewidths=[0.8, 1.0, 1.3, 1.6, 2.0])
    ax4.clabel(ct, inline=True, fontsize=7, fmt="%.0f Hz")
    ax4.set_xlabel("Surface wavelength $\\lambda$ [µm]", fontsize=12)
    ax4.set_ylabel("Ion height $h$ [µm]", fontsize=12)
    ax4.set_title(
        "Signal amplitude $|\\Delta\\omega_x|/2\\pi$ [Hz]\n"
        "$^{40}$Ca$^+$, $A=100$ nm, $f_{\\rm sec}=1$ MHz, $E_{\\rm bg}=10^4$ V/m",
        fontsize=12, fontweight="bold",
    )
    cbar = plt.colorbar(cs, ax=ax4, label="$\\log_{10}(|\\Delta\\omega_x|/2\\pi$ [Hz])")
    # Mark nominal operating point
    ax4.plot(lam_test * 1e6, trap.ion_height * 1e6, "r*", ms=18, mec="white", mew=2)
    ax4.text(lam_test * 1e6 * 1.3, trap.ion_height * 1e6 * 0.8,
             "Nominal\npoint", fontsize=9, color="white", fontweight="bold")
    # Optimal k = 2/h line
    h_line = np.logspace(-5.5, -3.5, 100)
    lam_opt = np.pi * h_line
    ax4.plot(lam_opt * 1e6, h_line * 1e6, "w--", lw=1, alpha=0.4)
    ax4.text(lam_opt[60] * 1e6, h_line[60] * 1e6 * 0.7,
             "$\\lambda_{\\rm opt}=\\pi h$", fontsize=8, color="white",
             rotation=-35, alpha=0.7)
    ax4.set_xscale("log")
    ax4.set_yscale("log")
    fig4.tight_layout()
    fig4.savefig("figS4_signal_map.pdf")
    print("  Saved figS4_signal_map.pdf")

    # ================================================================
    # Summary
    # ================================================================
    print(f"\n{'='*60}")
    print(f"  Simulation Complete")
    if np.sum(valid) >= 3:
        print(f"  ITF validation correlation: r = {corr:.4f}")
        if corr > 0.99:
            print(f"  ✓ Excellent agreement — perturbation theory validated.")
    print(f"\n  Figures: figS1_trajectory.pdf, figS2_scan_validation.pdf,")
    print(f"           figS3_linearity.pdf, figS4_signal_map.pdf")
