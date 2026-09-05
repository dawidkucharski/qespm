#!/usr/bin/env python3
"""
stability_and_roughness.py — Stability Diagram & Roughness Response
====================================================================
Two key additions for the QESPM theory paper:

  1. MATHIEU STABILITY DIAGRAM:
     How close can a rough surface approach before the ion is lost?
     Computes the effective Mathieu parameters (a_eff, q_eff) including
     the surface perturbation, and maps the stability boundary in
     (ion height h, roughness amplitude A) space.

  2. RESPONSE vs SURFACE ROUGHNESS PARAMETERS (Ra, Rz):
     How does the RMS frequency shift scale with ISO 25178 roughness
     parameters?  Generates random rough surfaces with controlled Ra,
     propagates through the ITF, and extracts the Δω_RMS vs Ra scaling.

Author: QESPM project
"""

import numpy as np
from numpy.fft import fft, ifft, fftfreq
from scipy.special import mathieu_a, mathieu_b
from scipy.optimize import brentq
import matplotlib
import matplotlib.pyplot as plt

# =====================================================================
# Physical constants
# =====================================================================
E_CHARGE = 1.602176634e-19
M_CA40   = 40.0 * 1.660539e-27
EPS0     = 8.8541878128e-12

# =====================================================================
# 1. MATHIEU STABILITY ANALYSIS
# =====================================================================

def mathieu_stability_boundary_a0(q):
    """Characteristic value a_0(q) — lower stability boundary."""
    # For q < 0.9, approximate: a_0 ≈ -q²/2 (pseudopotential)
    # For exact values we use scipy's mathieu_a
    if isinstance(q, np.ndarray):
        result = np.zeros_like(q)
        for i, qi in enumerate(q):
            if qi < 0.9:
                result[i] = mathieu_a(0, qi)  # characteristic value for order 0
            else:
                result[i] = mathieu_a(0, qi)
        return result
    else:
        return mathieu_a(0, q) if q < 0.9 else mathieu_a(0, q)


def mathieu_stability_boundary_b1(q):
    """Characteristic value b_1(q) — upper stability boundary."""
    if isinstance(q, np.ndarray):
        result = np.zeros_like(q)
        for i, qi in enumerate(q):
            result[i] = mathieu_b(1, qi)  # characteristic value for order 1
        return result
    else:
        return mathieu_b(1, q)


def surface_perturbation_delta_a(A, k, h, omega_rf, R_rf):
    """
    Effective shift in Mathieu 'a' parameter from surface roughness.

    The surface potential Φ = E_bg·A·sin(kx)·exp(-kh) adds a force gradient
    ∂F/∂x = e·E_bg·A·k²·sin(kx)·exp(-kh).

    In dimensionless Mathieu form, this adds to 'a':
    δa = (8e/mΩ²R²) · (∂Φ_perturbation/∂x linear term)

    For the worst case (sin(kx) = 1 at surface crest):
    δa_max = (8eE_bg A k² exp(-kh)) / (m Ω² R²)

    Parameters
    ----------
    A : surface amplitude [m]
    k : surface wavenumber [rad/m]
    h : ion height [m]
    omega_rf : rf angular frequency [rad/s]
    R_rf : rf electrode distance [m]
    """
    E_bg = 1e4  # V/m
    return (8.0 * E_CHARGE * E_bg * A * k**2 * np.exp(-k * h)) / (
        M_CA40 * omega_rf**2 * R_rf**2
    )


def compute_stability_boundary(
    omega_rf=2*np.pi*20e6, R_rf=500e-6,
    A_surface=100e-9, lam_surface=40e-6,
):
    """
    Compute the stability boundary in (h, A) space.

    Returns h_stable, A_max_at_h for plotting.
    """
    k_s = 2 * np.pi / lam_surface
    E_bg = 1e4

    # Nominal operating point (a_0, q_0) — well inside stability region
    q_0 = 4.0 * E_CHARGE * 200.0 / (M_CA40 * omega_rf**2 * R_rf**2)
    U_dc = 5.0
    d_dc = 2e-3
    kappa = 0.3
    a_0 = 8.0 * E_CHARGE * kappa * U_dc / (M_CA40 * omega_rf**2 * d_dc**2)

    # Stability boundaries at q_0
    a_lower = mathieu_a(0, q_0)  # a_0(q)
    a_upper = mathieu_b(1, q_0)  # b_1(q)

    # Surface perturbation shifts a: a_eff = a_0 + δa
    # Stability requires: a_lower < a_0 + δa_max < a_upper
    # (assuming δa is the worst-case positive shift)

    # For a sinusoidal surface, δa oscillates with sin(kx).
    # Maximum positive shift occurs at surface crests.
    # The condition δa_max < a_upper - a_0 gives the stability limit.

    h_vals = np.logspace(-5.5, -3, 200)  # 3 µm → 1 mm
    A_vals = np.logspace(-10, -5, 200)   # 0.1 nm → 10 µm

    # For each (h, A), compute δa and check stability
    stability_map = np.zeros((len(A_vals), len(h_vals)))
    for i, A_i in enumerate(A_vals):
        for j, h_j in enumerate(h_vals):
            da = surface_perturbation_delta_a(A_i, k_s, h_j, omega_rf, R_rf)
            # Check both upper and lower bounds
            a_eff_max = a_0 + da  # worst case (sin(kx)=1 at crest)
            a_eff_min = a_0 - da  # worst case (sin(kx)=-1 at trough)
            if a_lower < a_eff_min and a_eff_max < a_upper:
                stability_map[i, j] = 1

    return h_vals, A_vals, stability_map, a_0, q_0, a_lower, a_upper


# =====================================================================
# 2. RESPONSE vs SURFACE ROUGHNESS (Ra, Rz)
# =====================================================================

def generate_rough_surface_1d(n_points, dx, Ra_target, correlation_length,
                               rng=None):
    """
    Generate a 1D random rough surface with specified Ra.

    Uses a power-law spectrum S(k) ∝ (1 + (k·ℓ_c)²)^(-β) filtered to
    achieve the target Ra, where ℓ_c is the correlation length.

    Parameters
    ----------
    n_points : int
    dx : float, pixel size [m]
    Ra_target : float, target arithmetic mean roughness [m]
    correlation_length : float, correlation length [m]
    rng : np.random.Generator

    Returns
    -------
    z : (n_points,) ndarray, surface height [m]
    Ra_actual : float, achieved Ra [m]
    """
    if rng is None:
        rng = np.random.default_rng()

    # Generate white noise in Fourier space
    noise = rng.normal(0, 1, n_points)
    noise_hat = fft(noise)

    # Power-law filter
    k = 2 * np.pi * fftfreq(n_points, dx)
    k[0] = 1e-30  # avoid division by zero
    lc = correlation_length
    # Filter: S(k) ∝ 1 / (1 + (k·lc)^2)^(beta)
    beta = 1.5  # spectral exponent (1 < beta < 2 for realistic surfaces)
    H_k = 1.0 / (1.0 + (np.abs(k) * lc)**2)**(beta / 2.0)
    H_k[0] = 0.0  # DC = 0

    # Apply filter
    z_hat = noise_hat * H_k
    z = np.real(ifft(z_hat))

    # Normalise to target Ra
    Ra_current = np.mean(np.abs(z - np.mean(z)))
    if Ra_current > 1e-30:
        z = z * (Ra_target / Ra_current)

    Ra_actual = np.mean(np.abs(z - np.mean(z)))
    return z, Ra_actual


def compute_iso_parameters(z, dx):
    """Compute ISO 25178-2 height parameters for a 1D profile."""
    z_flat = z - np.mean(z)
    Ra = np.mean(np.abs(z_flat))
    Rq = np.std(z_flat)
    Rz_iso = np.max(z_flat) - np.min(z_flat)
    Rsk = np.mean(z_flat**3) / Rq**3 if Rq > 0 else 0
    Rku = np.mean(z_flat**4) / Rq**4 if Rq > 0 else 0
    return {'Ra': Ra, 'Rq': Rq, 'Rz': Rz_iso, 'Rsk': Rsk, 'Rku': Rku}


def propagate_through_itf(z_surface, dx, h, E_bg, omega_sec, ion_mass):
    """
    Compute Δω_x from surface profile z_surface using the ITF.

    Δω̃_x(k) = -(e·E_bg/2mω) · k² · e^{-k·h} · z̃(k)
    """
    n = len(z_surface)
    C = E_CHARGE * E_bg / (2.0 * ion_mass * omega_sec)
    k = 2 * np.pi * fftfreq(n, dx)
    k[0] = 1e-30

    z_hat = fft(z_surface)
    dw_hat = -C * k**2 * np.exp(-np.abs(k) * h) * z_hat
    dw_hat[0] = 0.0
    dw = np.real(ifft(dw_hat))
    return dw


# =====================================================================
# FIGURE GENERATION
# =====================================================================
if __name__ == "__main__":
    matplotlib.use("Agg")
    plt.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 150, "savefig.bbox": "tight",
        "font.family": "serif", "font.size": 11,
        "mathtext.fontset": "stix",
    })

    omega_rf = 2 * np.pi * 20e6
    R_rf = 500e-6
    omega_sec = 2 * np.pi * 1e6
    E_bg = 1e4

    print("=" * 60)
    print("  Stability Diagram & Roughness Response Analysis")
    print("=" * 60)

    # ================================================================
    # Figure R1: Mathieu Stability Diagram with Surface Perturbation
    # ================================================================
    print("\nComputing Mathieu stability diagram...")

    # Compute the full (a,q) stability region
    q_vals = np.linspace(0, 0.92, 200)
    a_lower_vals = np.array([mathieu_a(0, qi) for qi in q_vals])
    a_upper_vals = np.array([mathieu_b(1, qi) for qi in q_vals])

    # Nominal operating points for different trap configurations
    configs = [
        {"U_dc": 5.0, "V_rf": 200.0, "label": "Nominal", "color": "blue", "marker": "o"},
        {"U_dc": 2.0, "V_rf": 140.0, "label": "Low dc",  "color": "green", "marker": "s"},
        {"U_dc": 10.0, "V_rf": 260.0, "label": "High dc", "color": "orange", "marker": "D"},
    ]

    d_dc = 2e-3
    kappa = 0.3

    fig1, (ax1a, ax1b) = plt.subplots(1, 2, figsize=(16, 7))

    # (a) Full (a,q) stability diagram
    ax1a.fill_between(q_vals, a_lower_vals, a_upper_vals,
                       alpha=0.2, color="green", label="Stable region")
    ax1a.plot(q_vals, a_lower_vals, "k-", lw=1.5)
    ax1a.plot(q_vals, a_upper_vals, "k-", lw=1.5)

    for cfg in configs:
        U_dc = cfg["U_dc"]
        V_rf = cfg["V_rf"]
        q0 = 4.0 * E_CHARGE * V_rf / (M_CA40 * omega_rf**2 * R_rf**2)
        a0 = 8.0 * E_CHARGE * kappa * U_dc / (M_CA40 * omega_rf**2 * d_dc**2)
        ax1a.plot(q0, a0, cfg["marker"], color=cfg["color"], ms=12,
                   mec="white", mew=1.2, label=cfg["label"])
        # Arrow showing perturbation direction
        ax1a.annotate("", xy=(q0, a0 + 0.02), xytext=(q0, a0),
                       arrowprops=dict(arrowstyle="->", color=cfg["color"], lw=2))

    ax1a.set_xlabel("Mathieu parameter $q$", fontsize=12)
    ax1a.set_ylabel("Mathieu parameter $a$", fontsize=12)
    ax1a.set_title("(a) Mathieu stability diagram\n"
                   "Green: stable region. Arrows: surface perturbation direction",
                   fontsize=12, fontweight="bold")
    ax1a.legend(fontsize=9, loc="lower right", bbox_to_anchor=(0.99, 0.01))
    ax1a.grid(True, alpha=0.3)
    ax1a.set_xlim(0, 0.95)
    ax1a.set_ylim(-0.45, 0.25)

    # (b) Stability in (h, A) space
    lam_s = 40e-6
    k_s = 2 * np.pi / lam_s

    h_stab, A_stab, stab_map, a0_nom, q0_nom, a_low, a_up = \
        compute_stability_boundary(omega_rf, R_rf, 100e-9, lam_s)

    HH, AA = np.meshgrid(h_stab, A_stab)
    ax1b.contourf(HH * 1e6, AA * 1e9, stab_map,
                   levels=[0.5, 1.5], colors=["lightcoral", "lightgreen"],
                   alpha=0.4)
    # Stability boundary contour — manual label to avoid overlap
    ax1b.contour(HH * 1e6, AA * 1e9, stab_map,
                  levels=[0.5], colors=["red"], linewidths=2)
    ax1b.text(50, 200, "Stability limit", fontsize=9, color="red",
              fontweight="bold", rotation=-30, alpha=0.9)

    # Region labels — placed away from contour
    ax1b.text(15, 2000, "UNSTABLE", fontsize=11, color="darkred",
              fontweight="bold", ha="center", va="center", alpha=0.7)
    ax1b.text(120, 5, "STABLE", fontsize=11, color="darkgreen",
              fontweight="bold", ha="center", va="center", alpha=0.7)

    # Operating points at different heights for A=100nm
    h_ops = np.array([5, 10, 20, 40, 80, 160])
    for h_op in h_ops:
        da = surface_perturbation_delta_a(100e-9, k_s, h_op*1e-6, omega_rf, R_rf)
        a_eff = a0_nom + da
        stable = a_low < a0_nom - da and a0_nom + da < a_up
        marker = "o" if stable else "x"
        color = "green" if stable else "red"
        ax1b.plot(h_op, 100, marker, color=color, ms=10, mec="white", mew=1)

    ax1b.set_xscale("log")
    ax1b.set_yscale("log")
    ax1b.set_xlabel("Ion height $h$ [µm]", fontsize=12)
    ax1b.set_ylabel("Surface amplitude $A$ [nm]", fontsize=12)
    ax1b.set_title("(b) Stability in $(h, A)$ space\n"
                   f"$\\lambda = {lam_s*1e6:.0f}$ µm, "
                   f"$a_0={a0_nom:.4f}$, $q_0={q0_nom:.4f}$",
                   fontsize=12, fontweight="bold")
    ax1b.grid(True, alpha=0.3, which="both")

    fig1.suptitle("Mathieu Stability Under Surface Perturbation\n"
                  "$^{40}$Ca$^+$, $f_{\\rm rf}=20$ MHz, $R_{\\rm rf}=500$ µm",
                  fontsize=14, fontweight="bold")
    plt.tight_layout()
    fig1.savefig("figR1_stability.pdf")
    print("  Saved figR1_stability.pdf")

    # ================================================================
    # Figure R2: Response vs Surface Roughness Parameters
    # ================================================================
    print("Computing roughness response...")

    n_pts = 512
    dx = 0.5e-6
    lc = 5e-6  # correlation length
    ion_heights = np.array([5, 10, 20, 40, 80]) * 1e-6
    Ra_values = np.logspace(-10, -6, 15)  # 0.1 nm → 1 µm
    n_realisations = 10
    rng = np.random.default_rng(42)

    dw_rms = np.zeros((len(ion_heights), len(Ra_values)))
    dw_rms_std = np.zeros((len(ion_heights), len(Ra_values)))

    for j, h_j in enumerate(ion_heights):
        print(f"  h = {h_j*1e6:.0f} µm...")
        for i, Ra_i in enumerate(Ra_values):
            vals = np.zeros(n_realisations)
            for r in range(n_realisations):
                z_surf, _ = generate_rough_surface_1d(
                    n_pts, dx, Ra_i, lc, rng)
                dw = propagate_through_itf(
                    z_surf, dx, h_j, E_bg, omega_sec, M_CA40)
                vals[r] = np.std(dw) / (2 * np.pi)  # RMS in Hz
            dw_rms[j, i] = np.mean(vals)
            dw_rms_std[j, i] = np.std(vals)

    # Analytical prediction: for Ra ∝ A, σ_Δω ≈ |H(k_peak)| · σ_z
    # where k_peak ≈ 1/lc (dominant spatial frequency of roughness)
    k_peak = 2 * np.pi / lc  # correlation length → peak wavenumber
    C = E_CHARGE * E_bg / (2.0 * M_CA40 * omega_sec)
    dw_rms_analytic = np.zeros((len(ion_heights), len(Ra_values)))
    for j, h_j in enumerate(ion_heights):
        H_peak = C * k_peak**2 * np.exp(-k_peak * h_j)
        # Ra ≈ 0.8 × Rq for sinusoidal; σ_z ≈ Rq
        dw_rms_analytic[j, :] = H_peak * Ra_values * 0.8 / (2 * np.pi)

    fig2, (ax2a, ax2b) = plt.subplots(1, 2, figsize=(16, 7))
    colours = plt.cm.viridis(np.linspace(0.15, 0.9, len(ion_heights)))

    # (a) RMS frequency shift vs Ra
    for j, (h_j, c) in enumerate(zip(ion_heights, colours)):
        ax2a.errorbar(Ra_values * 1e9, dw_rms[j],
                       yerr=dw_rms_std[j], fmt="o-", color=c, ms=6,
                       capsize=3, lw=1.5,
                       label=f"$h = {h_j*1e6:.0f}$ µm")

    ax2a.set_xscale("log")
    ax2a.set_yscale("log")
    ax2a.set_xlabel("Arithmetic mean roughness $R_a$ [nm]", fontsize=12)
    ax2a.set_ylabel("RMS frequency shift $\\sigma_{\\Delta\\omega}/2\\pi$ [Hz]",
                    fontsize=12)
    ax2a.set_title("(a) Frequency shift response vs.~surface roughness\n"
                   f"Correlation length $\\ell_c = {lc*1e6:.0f}$ µm, "
                   f"{n_realisations} realisations each",
                   fontsize=12, fontweight="bold")
    ax2a.legend(fontsize=9, loc="lower right")
    ax2a.grid(True, alpha=0.3, which="both")

    # (b) Response normalised by Ra (sensitivity)
    ax2b_titles = []
    for j, (h_j, c) in enumerate(zip(ion_heights, colours)):
        sensitivity = dw_rms[j] / (Ra_values * 1e9)  # Hz/nm
        ax2b.semilogx(Ra_values * 1e9, sensitivity, "o-", color=c, ms=6, lw=1.5,
                       label=f"$h = {h_j*1e6:.0f}$ µm")

    ax2b.set_xlabel("Arithmetic mean roughness $R_a$ [nm]", fontsize=12)
    ax2b.set_ylabel("Sensitivity $\\sigma_{\\Delta\\omega}/R_a$ [Hz/nm]",
                    fontsize=12)
    ax2b.set_title("(b) Normalised sensitivity vs.~roughness\n"
                   "Constant sensitivity = linear regime",
                   fontsize=12, fontweight="bold")
    ax2b.legend(fontsize=9)
    ax2b.grid(True, alpha=0.3)

    fig2.suptitle("QESPM Response to Broadband Surface Roughness\n"
                  "$^{40}$Ca$^+$, $f_{\\rm sec}=1$ MHz, "
                  "$E_{\\rm bg}=10^4$ V/m",
                  fontsize=14, fontweight="bold")
    plt.tight_layout()
    fig2.savefig("figR2_roughness_response.pdf")
    print("  Saved figR2_roughness_response.pdf")

    # ================================================================
    # Figure R3: Example rough surfaces and their Δω response
    # ================================================================
    Ra_demo = np.array([1, 10, 100]) * 1e-9
    h_demo = 40e-6
    n_demo = 256

    fig3, axes3 = plt.subplots(len(Ra_demo), 2, figsize=(14, 10))

    for i, Ra_i in enumerate(Ra_demo):
        z_s, Ra_act = generate_rough_surface_1d(n_demo, dx, Ra_i, lc, rng)
        dw = propagate_through_itf(z_s, dx, h_demo, E_bg, omega_sec, M_CA40)
        x_um = np.arange(n_demo) * dx * 1e6

        # Surface profile
        ax_z = axes3[i, 0]
        ax_z.plot(x_um, z_s * 1e9, "k-", lw=0.8)
        ax_z.fill_between(x_um, 0, z_s * 1e9, alpha=0.3)
        ax_z.set_ylabel("$z$ [nm]", fontsize=10)
        ax_z.set_title(f"$R_a = {Ra_act*1e9:.1f}$ nm, "
                       f"$R_z = {(np.max(z_s)-np.min(z_s))*1e9:.0f}$ nm",
                       fontsize=10)
        ax_z.grid(True, alpha=0.3)

        # Frequency shift
        ax_dw = axes3[i, 1]
        ax_dw.plot(x_um, dw / (2*np.pi), "r-", lw=0.8)
        ax_dw.set_ylabel("$\\Delta\\omega_x/2\\pi$ [Hz]", fontsize=10)
        ax_dw.set_title(f"$\\sigma_{{\\Delta\\omega}}/2\\pi = "
                        f"{np.std(dw)/(2*np.pi):.1f}$ Hz",
                        fontsize=10)
        ax_dw.grid(True, alpha=0.3)

    axes3[-1, 0].set_xlabel("$x$ [µm]", fontsize=11)
    axes3[-1, 1].set_xlabel("$x$ [µm]", fontsize=11)

    fig3.suptitle(f"Surface Profiles and QESPM Response\n"
                  f"$h = {h_demo*1e6:.0f}$ µm, "
                  f"$\\ell_c = {lc*1e6:.0f}$ µm",
                  fontsize=14, fontweight="bold")
    plt.tight_layout()
    fig3.savefig("figR3_example_profiles.pdf")
    print("  Saved figR3_example_profiles.pdf")

    # ================================================================
    # Numerical summary
    # ================================================================
    print(f"\n{'='*60}")
    print(f"  Stability & Roughness Analysis Complete")
    print(f"{'='*60}")

    # Stability summary
    h_test_vals = np.array([3, 5, 7, 10, 15, 20, 30, 50, 80, 150]) * 1e-6
    print(f"\n  Stability check for A=100nm, lambda=40um:")
    print(f"  {'h [µm]':>8s}  {'δa':>10s}  {'a_eff':>10s}  {'Stable?':>8s}")
    for h_t in h_test_vals:
        da = surface_perturbation_delta_a(100e-9, k_s, h_t, omega_rf, R_rf)
        a_eff = a0_nom + da
        stable = "YES" if (a_low < a0_nom - abs(da) and a0_nom + abs(da) < a_up) else "NO"
        print(f"  {h_t*1e6:8.1f}  {da:10.2e}  {a_eff:10.4f}  {stable:>8s}")

    # Roughness response summary
    print(f"\n  Roughness response: RMS Δω [Hz] vs Ra [nm] at h=40um:")
    idx_40 = np.argmin(np.abs(ion_heights - 40e-6))
    for i, Ra_i in enumerate(Ra_values):
        print(f"    Ra = {Ra_i*1e9:8.2f} nm  →  "
              f"σ_Δω/2π = {dw_rms[idx_40, i]:8.2f} Hz")

    print(f"\n  Figures: figR1_stability.pdf, figR2_roughness_response.pdf,")
    print(f"           figR3_example_profiles.pdf")
