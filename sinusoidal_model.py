"""
sinusoidal_model.py — Simplified Analytical Model
===================================================
Single trapped Ca⁺ ion as a probe of sinusoidal surface roughness.

Geometry:
  Surface:  z(x) = A sin(kx)   [conducting, grounded]
  Ion:      trapped at height h above the mean surface plane
  Trap:     linear Paul trap, secular frequency ω_sec/2π = 1 MHz

Physics chain:
  z(x) → φ_s(x) → Φ(x,h) → ∂²Φ/∂x² → Δω(x) → measured signal

All quantities in SI units unless otherwise noted.

Author: Research simulation
"""

import numpy as np
from numpy.fft import fftfreq
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from dataclasses import dataclass

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
E_CHARGE = 1.602176634e-19     # C
M_CA40   = 40.0 * 1.660539e-27  # kg  (⁴⁰Ca⁺)
EPS0     = 8.8541878128e-12     # F/m

# ---------------------------------------------------------------------------
# Experimental parameters
# ---------------------------------------------------------------------------
@dataclass
class TrapParams:
    """Nominal experimental parameters for a ⁴⁰Ca⁺ surface trap."""
    ion_mass: float = M_CA40            # kg
    freq_sec: float = 1.0e6             # secular frequency, Hz (ω/2π)
    E_bg: float = 1.0e4                 # background field at surface, V/m
    delta_omega_min: float = 2.0 * np.pi * 5.0  # freq. resolution, rad/s (5 Hz conservative)

    @property
    def omega_sec(self) -> float:
        return 2.0 * np.pi * self.freq_sec

    @property
    def prefactor(self) -> float:
        """e / (2 m ω_sec) — converts ∂²Φ/∂x² to Δω [rad/s per V/m²]"""
        return E_CHARGE / (2.0 * self.ion_mass * self.omega_sec)


# ---------------------------------------------------------------------------
# Analytical model: single sinusoidal mode
# ---------------------------------------------------------------------------

def surface_potential_amplitude(A: float, E_bg: float) -> float:
    """Surface potential amplitude: φ_s = E_bg · A  [V]"""
    return E_bg * A


def potential_at_ion(
    x: np.ndarray,
    A: float,
    k: float,
    h: float,
    E_bg: float,
) -> np.ndarray:
    """
    Electrostatic potential at ion height h above a sinusoidal surface.

    Φ(x, h) = E_bg · A · sin(kx) · exp(-k h)

    This is the exact solution of Laplace's equation ∇²Φ = 0 in the
    half-space z > 0 with Dirichlet boundary condition
    Φ(x, 0) = E_bg · A · sin(kx).

    Parameters
    ----------
    x : ndarray, positions along surface [m]
    A : float, surface roughness amplitude [m]
    k : float, surface wavenumber [rad/m]  (k = 2π/λ)
    h : float, ion height above mean surface [m]
    E_bg : float, background electric field [V/m]

    Returns
    -------
    Phi : ndarray, potential at (x, h) [V]
    """
    return E_bg * A * np.sin(k * x) * np.exp(-k * h)


def secular_frequency_shift(
    x: np.ndarray,
    A: float,
    k: float,
    h: float,
    params: TrapParams,
) -> np.ndarray:
    """
    Secular frequency shift due to surface potential.

    Δω_x(x) = (e / 2m ω_sec) · ∂²Φ/∂x²
            = -(e / 2m ω_sec) · E_bg · A · k² · sin(kx) · exp(-k h)

    Parameters
    ----------
    x, A, k, h : as above
    params : TrapParams

    Returns
    -------
    delta_omega : ndarray, frequency shift [rad/s]
    """
    prefactor = params.prefactor * params.E_bg * A * k**2 * np.exp(-k * h)
    return -prefactor * np.sin(k * x)


def frequency_shift_amplitude(
    A: float,
    k: float,
    h: float,
    params: TrapParams,
) -> float:
    """
    Amplitude of secular frequency shift oscillation.

    |Δω_x| = (e / 2m ω_sec) · E_bg · A · k² · exp(-k h)

    Returns scalar [rad/s].
    """
    return params.prefactor * params.E_bg * A * k**2 * np.exp(-k * h)


def minimum_detectable_amplitude(
    k: np.ndarray,
    h: float,
    params: TrapParams,
) -> np.ndarray:
    """
    Minimum surface roughness amplitude detectable at SNR = 1.

    From |Δω_x| = δω_min:
        A_min = (2m ω_sec / e) · δω_min / (E_bg · k² · exp(-k h))

    Returns A_min [m].  Returns inf where k ≈ 0.
    """
    k = np.asarray(k, dtype=float)
    with np.errstate(divide="ignore", over="ignore"):
        result = (params.delta_omega_min / params.prefactor) / (
            params.E_bg * k**2 * np.exp(-k * h)
        )
    result = np.where(np.isfinite(result), result, np.inf)
    return result


def optimal_wavenumber(h: float) -> float:
    """
    Wavenumber giving maximum frequency shift for a given height.

    Maximises k² · exp(-k h):
        d/dk [k² exp(-k h)] = 0  →  k_opt = 2/h

    Returns k_opt [rad/m].
    """
    return 2.0 / h


def response_function(k: float, h: float) -> float:
    """
    Spatial frequency response: R(k; h) = k² · exp(-k h).

    This is the transfer function from surface potential amplitude to
    frequency shift amplitude (up to the prefactor).

    Returns dimensionless response (× k² has units of 1/m²).
    """
    return k**2 * np.exp(-k * h)


def snr(A: float, k: float, h: float, params: TrapParams) -> float:
    """Signal-to-noise ratio for a single sinusoidal mode."""
    return frequency_shift_amplitude(A, k, h, params) / params.delta_omega_min


# ---------------------------------------------------------------------------
# Reconstruction analysis
# ---------------------------------------------------------------------------

def reconstruction_band(
    h: float,
    A: float,
    params: TrapParams,
    snr_threshold: float = 1.0,
    k_min: float = 1.0,
    k_max: float = 1e8,
    n_k: int = 10000,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """
    Determine the band of reconstructable wavenumbers.

    A spatial frequency k is reconstructable if SNR(k) > snr_threshold.

    Parameters
    ----------
    h : float, ion height [m]
    A : float, surface roughness amplitude [m]
    params : TrapParams
    snr_threshold : float, minimum SNR for detection
    k_min, k_max : float, search range [rad/m]
    n_k : int, number of k samples

    Returns
    -------
    k_vals : ndarray, wavenumbers
    snr_vals : ndarray, SNR at each k
    k_low : float, lower bound of reconstructable band (or NaN)
    k_high : float, upper bound of reconstructable band (or NaN)
    """
    k_vals = np.logspace(np.log10(k_min), np.log10(k_max), n_k)
    snr_vals = snr(A, k_vals, h, params)

    detectable = snr_vals >= snr_threshold
    if not np.any(detectable):
        return k_vals, snr_vals, np.nan, np.nan

    idx = np.where(detectable)[0]
    k_low = k_vals[idx[0]]
    k_high = k_vals[idx[-1]]
    return k_vals, snr_vals, k_low, k_high


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def set_pub_style(ax, xlabel="", ylabel="", title="", grid=True):
    """Apply consistent publication-ready styling."""
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=13, fontweight="bold")
    if grid:
        ax.grid(True, alpha=0.3, lw=0.5)
    ax.tick_params(labelsize=10)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)


# ===================================================================
# Main: generate all figures
# ===================================================================
if __name__ == "__main__":
    matplotlib.use("Agg")
    plt.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 150,
        "savefig.bbox": "tight",
        "font.family": "serif",
        "font.size": 11,
        "mathtext.fontset": "stix",
    })

    params = TrapParams()
    A_nom = 100e-9      # 100 nm roughness amplitude
    k_nom = 2e5          # k = 2π/λ with λ ≈ 31 µm
    h_nom = 40e-6        # 40 µm ion height

    print("=" * 65)
    print("  Simplified Analytical Model — Sinusoidal Surface")
    print("=" * 65)
    print(f"  Ion:  ⁴⁰Ca⁺, m = 40 amu")
    print(f"  Secular frequency:  ω_sec/2π = {params.freq_sec/1e6:.1f} MHz")
    print(f"  Freq. resolution:   δω = {params.delta_omega_min:.1f} rad/s")
    print(f"                       = {params.delta_omega_min/(2*np.pi):.2f} Hz")
    print(f"  Background field:   E_bg = {params.E_bg:.0e} V/m")
    print(f"  Prefactor e/(2mω):  {params.prefactor:.4f} (V·s²/m²·rad)⁻¹? "
          f"= {params.prefactor:.4f} m²/(V·s²)")
    print()

    # --- Nominal case ---
    S_nom = frequency_shift_amplitude(A_nom, k_nom, h_nom, params)
    k_opt = optimal_wavenumber(h_nom)
    S_opt = frequency_shift_amplitude(A_nom, k_opt, h_nom, params)
    A_min_nom = minimum_detectable_amplitude(k_nom, h_nom, params)
    SNR_nom = snr(A_nom, k_nom, h_nom, params)

    print(f"--- Nominal case ---")
    print(f"  A = {A_nom*1e9:.0f} nm,  k = {k_nom:.1e} rad/m")
    print(f"  λ = 2π/k = {2*np.pi/k_nom*1e6:.1f} µm")
    print(f"  h = {h_nom*1e6:.0f} µm")
    print(f"  k·h = {k_nom*h_nom:.2f}")
    print(f"  exp(-k·h) = {np.exp(-k_nom*h_nom):.2e}")
    print(f"  |Δω| = {S_nom:.2f} rad/s  = {S_nom/(2*np.pi):.2f} Hz")
    print(f"  SNR  = {SNR_nom:.1f}")
    print(f"  A_min (SNR=1) = {A_min_nom*1e9:.2f} nm")
    print(f"  Optimal k = 2/h = {k_opt:.1e} rad/m  (λ = {2*np.pi/k_opt*1e6:.1f} µm)")
    print(f"  |Δω| at k_opt = {S_opt:.2f} rad/s  = {S_opt/(2*np.pi):.2f} Hz")
    print(f"  SNR at k_opt = {snr(A_nom, k_opt, h_nom, params):.1f}")
    print()

    # ================================================================
    # Figure 1: Single-mode spatial profiles
    # ================================================================
    x_span = 4 * 2 * np.pi / k_nom  # 4 periods
    x = np.linspace(0, x_span, 1000)

    z_surface = A_nom * np.sin(k_nom * x)
    Phi_ion = potential_at_ion(x, A_nom, k_nom, h_nom, params.E_bg)
    dw = secular_frequency_shift(x, A_nom, k_nom, h_nom, params)

    fig1, axes1 = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

    # Panel (a): surface topography
    ax = axes1[0]
    ax.plot(x * 1e6, z_surface * 1e9, "k-", lw=1.5)
    ax.axhline(0, color="gray", ls="--", lw=0.5)
    set_pub_style(ax, "", "$z(x)$ [nm]",
                  f"(a) Surface topography  $z(x) = A \\sin(kx)$")
    ax.text(0.02, 0.9, f"$A = {A_nom*1e9:.0f}$ nm,  "
            f"$\\lambda = {2*np.pi/k_nom*1e6:.0f}$ µm",
            transform=ax.transAxes, fontsize=10, va="top",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.7))

    # Panel (b): potential at ion height
    ax = axes1[1]
    ax.plot(x * 1e6, Phi_ion * 1e6, "b-", lw=1.5)
    ax.axhline(0, color="gray", ls="--", lw=0.5)
    set_pub_style(ax, "", "$\\Phi(x,h)$ [$\\mu$V]",
                  f"(b) Electrostatic potential at ion height  $h = {h_nom*1e6:.0f}$ µm")
    ax.text(0.02, 0.9, f"$\\Phi = E_{{\\rm bg}} A \\, e^{{-kh}} \\sin(kx)$\n"
            f"$e^{{-kh}} = {np.exp(-k_nom*h_nom):.1e}$",
            transform=ax.transAxes, fontsize=10, va="top",
            bbox=dict(boxstyle="round", facecolor="lightblue", alpha=0.7))

    # Panel (c): secular frequency shift
    ax = axes1[2]
    ax.plot(x * 1e6, dw / (2 * np.pi), "r-", lw=1.5)
    ax.axhline(0, color="gray", ls="--", lw=0.5)
    # Show noise floor
    noise_hz = params.delta_omega_min / (2 * np.pi)
    ax.axhline(+noise_hz, color="gray", ls=":", lw=0.8, alpha=0.5)
    ax.axhline(-noise_hz, color="gray", ls=":", lw=0.8, alpha=0.5,
               label=f"noise floor ±{noise_hz:.1f} Hz")
    ax.legend(fontsize=9, loc="upper right")
    set_pub_style(ax, "$x$ [µm]", "$\\Delta\\omega_x / 2\\pi$ [Hz]",
                  "(c) Secular frequency shift  $\\Delta\\omega_x(x)$")
    ax.text(0.02, 0.9, f"$|\\Delta\\omega_x| = {S_nom/(2*np.pi):.1f}$ Hz\n"
            f"SNR = {SNR_nom:.1f}",
            transform=ax.transAxes, fontsize=10, va="top",
            bbox=dict(boxstyle="round", facecolor="lightcoral", alpha=0.7))

    fig1.suptitle(
        "Single-Ion Sinusoidal Surface Sensing — Spatial Profiles\n"
        f"$^{{40}}$Ca$^+$, $\\omega_{{\\rm sec}}/2\\pi = {params.freq_sec/1e6:.0f}$ MHz, "
        f"$h = {h_nom*1e6:.0f}$ µm, "
        f"$A = {A_nom*1e9:.0f}$ nm, "
        f"$k = {k_nom:.1e}$ rad/m",
        fontsize=14, fontweight="bold", y=1.01,
    )
    plt.tight_layout()
    fig1.savefig("fig1_spatial_profiles.pdf")
    print("Saved fig1_spatial_profiles.pdf")

    # ================================================================
    # Figure 2: Transfer function R(k; h) = k² exp(-kh)
    # ================================================================
    heights = np.array([1, 5, 10, 40, 100]) * 1e-6  # m
    k_vals = np.logspace(3, 8, 1000)  # 10³ to 10⁸ rad/m  (λ: 6 mm → 60 nm)
    colours = plt.cm.viridis(np.linspace(0.15, 0.9, len(heights)))

    fig2, (ax2a, ax2b) = plt.subplots(1, 2, figsize=(14, 6))

    for h_i, c in zip(heights, colours):
        R = response_function(k_vals, h_i)
        ax2a.loglog(k_vals, R, color=c, lw=1.8,
                     label=f"$h = {h_i*1e6:.0f}$ µm")

        # Mark optimum
        k_opt_i = optimal_wavenumber(h_i)
        R_opt = response_function(k_opt_i, h_i)
        ax2a.plot(k_opt_i, R_opt, "o", color=c, ms=7, mec="white", mew=0.8)

    ax2a.set_xlabel("Wavenumber $k$ [rad/m]", fontsize=12)
    ax2a.set_ylabel("Response $R(k;h) = k^2 e^{-kh}$ [m$^{-2}$]", fontsize=12)
    ax2a.set_title("(a) Spatial frequency transfer function", fontsize=13,
                   fontweight="bold")
    ax2a.legend(fontsize=9, loc="lower left")
    ax2a.grid(True, alpha=0.3, lw=0.5, which="both")
    # Secondary x-axis: wavelength
    ax2a_top = ax2a.twiny()
    ax2a_top.set_xscale("log")
    ax2a_top.set_xlim(ax2a.get_xlim())
    lambda_ticks = np.array([0.1, 1, 10, 100, 1000]) * 1e-6
    k_ticks = 2 * np.pi / lambda_ticks
    ax2a_top.set_xticks(k_ticks)
    ax2a_top.set_xticklabels([f"{lt*1e6:.0f}" for lt in lambda_ticks])
    ax2a_top.set_xlabel("Wavelength $\\lambda = 2\\pi/k$ [µm]", fontsize=11)

    # Panel (b): Optimal k vs h and peak response vs h
    h_range = np.logspace(-6, -3, 200)  # 1 µm → 1 mm
    k_opt_range = optimal_wavenumber(h_range)
    R_peak = response_function(k_opt_range, h_range)

    ax2b_twin = ax2b.twinx()
    ax2b.loglog(h_range * 1e6, k_opt_range, "b-", lw=2,
                label="$k_{\\rm opt} = 2/h$")
    ax2b_twin.loglog(h_range * 1e6, R_peak, "r-", lw=2,
                     label="$R_{\\rm max} = 4e^{-2}/h^2$")
    ax2b.set_xlabel("Ion height $h$ [µm]", fontsize=12)
    ax2b.set_ylabel("Optimal wavenumber $k_{\\rm opt}$ [rad/m]", fontsize=12,
                    color="b")
    ax2b_twin.set_ylabel("Peak response $R_{\\rm max}$ [m$^{-2}$]", fontsize=12,
                         color="r")
    ax2b.set_title("(b) Optimal wavenumber and peak response", fontsize=13,
                   fontweight="bold")
    ax2b.grid(True, alpha=0.3, lw=0.5, which="both")
    ax2b.tick_params(axis="y", colors="b")
    ax2b_twin.tick_params(axis="y", colors="r")
    # Secondary x-axis for k_opt → λ_opt
    ax2b_top = ax2b.twiny()
    ax2b_top.set_xscale("log")
    ax2b_top.set_xlim(ax2b.get_xlim())
    ax2b_top.set_xlabel("")
    # Add wavelength labels at a few key heights
    for h_ref in [1e-6, 10e-6, 100e-6]:
        k_ref = 2.0 / h_ref
        lam_ref = 2 * np.pi / k_ref
        ax2b.axvline(h_ref * 1e6, color="gray", ls=":", alpha=0.3)
        ax2b.annotate(f"$\\lambda_{{\\rm opt}}$={lam_ref*1e6:.0f} µm",
                      (h_ref * 1e6, k_ref),
                      textcoords="offset points", xytext=(5, 10),
                      fontsize=8, color="gray",
                      arrowprops=dict(arrowstyle="->", color="gray", alpha=0.5))

    fig2.suptitle(
        "Spatial Frequency Response of the Ion Probe",
        fontsize=14, fontweight="bold",
    )
    plt.tight_layout()
    fig2.savefig("fig2_transfer_function.pdf")
    print("Saved fig2_transfer_function.pdf")

    # ================================================================
    # Figure 3: Detectability map in (k, A) space for various h
    # ================================================================
    fig3, axes3 = plt.subplots(2, 3, figsize=(18, 11))
    h_list = [1e-6, 5e-6, 10e-6, 40e-6, 100e-6, 200e-6]
    k_map = np.logspace(3, 8, 300)
    A_map = np.logspace(-11, -5, 300)  # 10 pm → 10 µm
    KK, AA = np.meshgrid(k_map, A_map)

    for ax, h_i in zip(axes3.ravel(), h_list):
        SNR_map = snr(AA, KK, h_i, params)

        # Contour at SNR = 1, 10, 100
        levels = [1, 10, 100, 1000]
        cs = ax.contour(KK, AA, SNR_map, levels=levels,
                         colors=["white", "yellow", "orange", "red"],
                         linewidths=[0.8, 1.2, 1.6, 2.0])
        ax.clabel(cs, inline=True, fontsize=7, fmt="SNR=%.0f")

        # Shaded region: SNR > 1
        ax.contourf(KK, AA, SNR_map, levels=[1, 1e6],
                     colors=["lightgreen"], alpha=0.25)

        # Optimal k line
        k_opt_i = optimal_wavenumber(h_i)
        ax.axvline(k_opt_i, color="blue", ls="--", lw=1.2, alpha=0.7,
                    label=f"$k_{{\\rm opt}} = 2/h$")

        # Mark nominal case
        if abs(h_i - h_nom) < 1e-9:
            ax.plot(k_nom, A_nom, "r*", ms=12, mec="white", mew=1,
                    zorder=10, label="nominal")

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("$k$ [rad/m]", fontsize=10)
        ax.set_ylabel("$A$ [m]", fontsize=10)
        ax.set_title(f"$h = {h_i*1e6:.0f}$ µm", fontsize=11, fontweight="bold")
        ax.grid(True, alpha=0.3, lw=0.5, which="both")
        if h_i == 1e-6:
            ax.legend(fontsize=8, loc="lower right")

        # Add λ ticks on top
        ax_top = ax.twiny()
        ax_top.set_xscale("log")
        ax_top.set_xlim(ax.get_xlim())
        lam_ticks = [0.1e-6, 1e-6, 10e-6, 100e-6, 1e-3]
        ax_top.set_xticks([2 * np.pi / lt for lt in lam_ticks])
        ax_top.set_xticklabels([f"{lt*1e6:.0f}" for lt in lam_ticks])
        if h_i == h_list[0]:
            ax_top.set_xlabel("$\\lambda$ [µm]", fontsize=10)

    fig3.suptitle(
        "Detectability Map: SNR$(k, A; h)$\n"
        "Green shaded region = reconstructable ($\\rm SNR \\geq 1$)",
        fontsize=14, fontweight="bold",
    )
    plt.tight_layout()
    fig3.savefig("fig3_detectability_map.pdf")
    print("Saved fig3_detectability_map.pdf")

    # ================================================================
    # Figure 4: Minimum detectable amplitude and reconstructable band
    # ================================================================
    fig4, (ax4a, ax4b) = plt.subplots(1, 2, figsize=(14, 6))

    # Panel (a): A_min(k) for various h
    h_vals_Amin = [5e-6, 10e-6, 40e-6, 100e-6]
    k_vals_4 = np.logspace(3, 8, 500)
    colours_4 = plt.cm.plasma(np.linspace(0.15, 0.9, len(h_vals_Amin)))

    for h_i, c in zip(h_vals_Amin, colours_4):
        A_min = minimum_detectable_amplitude(k_vals_4, h_i, params)
        ax4a.loglog(k_vals_4, A_min * 1e9, color=c, lw=1.8,
                     label=f"$h = {h_i*1e6:.0f}$ µm")
        k_opt_i = optimal_wavenumber(h_i)
        A_min_opt = minimum_detectable_amplitude(k_opt_i, h_i, params)
        ax4a.plot(k_opt_i, A_min_opt * 1e9, "o", color=c, ms=8,
                   mec="white", mew=0.8)

    # Reference: typical AFM roughness
    ax4a.axhline(0.1, color="gray", ls=":", alpha=0.5)
    ax4a.text(1e3, 0.12, "AFM resolution (~1 Å)", fontsize=8, color="gray")
    # Reference: optical profilometry
    ax4a.axhline(1, color="gray", ls=":", alpha=0.5)
    ax4a.text(1e3, 1.1, "Optical profilometry (~1 nm)", fontsize=8, color="gray")

    ax4a.set_xlabel("Wavenumber $k$ [rad/m]", fontsize=12)
    ax4a.set_ylabel("Minimum detectable $A_{\\rm min}$ [nm]", fontsize=12)
    ax4a.set_title("(a) Detection threshold  $A_{\\rm min}(k; h)$", fontsize=13,
                   fontweight="bold")
    ax4a.legend(fontsize=9, loc="upper left")
    ax4a.grid(True, alpha=0.3, lw=0.5, which="both")
    ax4a.set_ylim(1e-5, 1e5)

    # Panel (b): Reconstructable band as function of h for fixed A
    A_values = [1e-9, 10e-9, 100e-9, 1e-6]  # 1 nm → 1 µm
    h_range_4b = np.logspace(-6, -3, 200)
    colours_4b = plt.cm.cividis(np.linspace(0.15, 0.9, len(A_values)))

    for A_i, c in zip(A_values, colours_4b):
        k_lo_vals = np.full_like(h_range_4b, np.nan)
        k_hi_vals = np.full_like(h_range_4b, np.nan)

        for j, h_j in enumerate(h_range_4b):
            _, _, k_lo, k_hi = reconstruction_band(
                h_j, A_i, params,
                k_min=2 * np.pi / 200e-6,  # λ_max = 200 µm (FOV-limited)
            )
            k_lo_vals[j] = k_lo
            k_hi_vals[j] = k_hi

        # Convert k to λ for plotting
        lam_lo = np.where(np.isfinite(k_hi_vals), 2 * np.pi / k_hi_vals, np.nan)
        lam_hi = np.where(np.isfinite(k_lo_vals), 2 * np.pi / k_lo_vals, np.nan)

        # Plot the band
        valid = np.isfinite(lam_lo) & np.isfinite(lam_hi)
        if np.any(valid):
            ax4b.fill_between(
                h_range_4b[valid] * 1e6,
                lam_lo[valid] * 1e6,
                lam_hi[valid] * 1e6,
                color=c, alpha=0.25,
            )
            # Upper and lower boundary lines
            ax4b.loglog(h_range_4b[valid] * 1e6, lam_lo[valid] * 1e6,
                        color=c, lw=1.5, ls="-")
            ax4b.loglog(h_range_4b[valid] * 1e6, lam_hi[valid] * 1e6,
                        color=c, lw=1.5, ls="-",
                        label=f"$A = {A_i*1e9:.0f}$ nm")

        # Mark the nominal point
        if abs(A_i - A_nom) < 1e-12:
            k_lo_nom, k_hi_nom = None, None
            for j, h_j in enumerate(h_range_4b):
                if abs(h_j - h_nom) < 5e-8:
                    _, _, k_lo_nom, k_hi_nom = reconstruction_band(
                        h_j, A_i, params,
                        k_min=2 * np.pi / 200e-6,
                    )
                    break
            if k_lo_nom is not None and np.isfinite(k_lo_nom):
                ax4b.plot(h_nom * 1e6, 2 * np.pi / k_hi_nom * 1e6,
                          "r*", ms=12, zorder=10)
                ax4b.plot(h_nom * 1e6, 2 * np.pi / k_lo_nom * 1e6,
                          "r*", ms=12, zorder=10)

    ax4b.set_xlabel("Ion height $h$ [µm]", fontsize=12)
    ax4b.set_ylabel("Reconstructable $\\lambda$ band [µm]", fontsize=12)
    ax4b.set_title("(b) Reconstructable wavelength band vs. height", fontsize=13,
                   fontweight="bold")
    ax4b.legend(fontsize=9, loc="lower right")
    ax4b.grid(True, alpha=0.3, lw=0.5, which="both")

    fig4.suptitle(
        "Reconstruction Feasibility: Detection Thresholds and Accessible "
        "Spatial Frequencies",
        fontsize=14, fontweight="bold",
    )
    plt.tight_layout()
    fig4.savefig("fig4_reconstruction_conditions.pdf")
    print("Saved fig4_reconstruction_conditions.pdf")

    # ================================================================
    # Figure 5: Signal vs. height — practical measurement regime
    # ================================================================
    fig5, ax5 = plt.subplots(1, 1, figsize=(10, 6))

    h_practical = np.logspace(-6, -3, 300)
    k_practical = [2 * np.pi / lam for lam in
                   [1e-6, 5e-6, 10e-6, 50e-6, 100e-6]]  # λ = 1–100 µm
    A_practical = 100e-9  # 100 nm
    colours_5 = plt.cm.turbo(np.linspace(0.1, 0.9, len(k_practical)))

    for k_i, c in zip(k_practical, colours_5):
        dw_hz = frequency_shift_amplitude(A_practical, k_i, h_practical, params)
        dw_hz /= (2 * np.pi)  # convert to Hz
        lam_i = 2 * np.pi / k_i
        ax5.loglog(h_practical * 1e6, dw_hz, color=c, lw=1.8,
                    label=f"$\\lambda = {lam_i*1e6:.0f}$ µm")

    # Noise floor
    ax5.axhline(params.delta_omega_min / (2 * np.pi), color="black",
                ls="--", lw=1.2, alpha=0.6)
    ax5.text(1.5, params.delta_omega_min / (2 * np.pi) * 1.5,
             "$\\delta\\omega$ noise floor (1 Hz)",
             fontsize=9, color="black")

    # Quantum projection noise limit (lower bound)
    # δω_QPN ≈ 1/(2η√N) — for η=0.1, N=1000, τ=1ms
    qpn_limit = 0.01  # Hz (aspirational)
    ax5.axhline(qpn_limit, color="red", ls=":", lw=1, alpha=0.5)
    ax5.text(80, qpn_limit * 1.5, "QPN limit (~0.01 Hz, aspirational)",
             fontsize=8, color="red", alpha=0.7)

    # Shading: "detectable" region
    ax5.fill_between(h_practical * 1e6,
                      params.delta_omega_min / (2 * np.pi), 1e10,
                      color="green", alpha=0.07)

    ax5.set_xlabel("Ion height $h$ [µm]", fontsize=12)
    ax5.set_ylabel("$|\\Delta\\omega_x| / 2\\pi$ [Hz]", fontsize=12)
    ax5.set_title(
        "Measurable Frequency Shift vs. Ion Height\n"
        f"$A = {A_practical*1e9:.0f}$ nm surface roughness",
        fontsize=13, fontweight="bold",
    )
    ax5.legend(fontsize=9, loc="lower left", ncol=2)
    ax5.grid(True, alpha=0.3, lw=0.5, which="both")
    ax5.set_ylim(1e-4, 1e8)

    fig5.tight_layout()
    fig5.savefig("fig5_signal_vs_height.pdf")
    print("Saved fig5_signal_vs_height.pdf")

    # ================================================================
    # Numerical summary
    # ================================================================
    print()
    print("=" * 65)
    print("  Key Numerical Results")
    print("=" * 65)
    print()
    print("  Transfer function peak:")
    for h_i in [1e-6, 5e-6, 10e-6, 40e-6, 100e-6]:
        k_opt_i = optimal_wavenumber(h_i)
        R_opt = response_function(k_opt_i, h_i)
        lam_opt = 2 * np.pi / k_opt_i
        A_min = minimum_detectable_amplitude(k_opt_i, h_i, params)
        S_opt = frequency_shift_amplitude(A_nom, k_opt_i, h_i, params)
        print(f"    h = {h_i*1e6:5.0f} µm  →  "
              f"λ_opt = {lam_opt*1e6:6.1f} µm,  "
              f"R_max = {R_opt:.2e} m⁻²,  "
              f"A_min = {A_min*1e9:.2e} nm,  "
              f"|Δω|(A={A_nom*1e9:.0f}nm) = {S_opt/(2*np.pi):.1f} Hz")

    print()
    print("  Reconstructable bands for A = 100 nm:")
    for h_i in [5e-6, 10e-6, 40e-6, 100e-6]:
        _, _, k_lo, k_hi = reconstruction_band(
            h_i, A_nom, params, k_min=2 * np.pi / 200e-6,
        )
        if np.isfinite(k_lo):
            lam_lo = 2 * np.pi / k_hi
            lam_hi = 2 * np.pi / k_lo
            print(f"    h = {h_i*1e6:5.0f} µm  →  "
                  f"λ ∈ [{lam_lo*1e6:.1f}, {lam_hi*1e6:.1f}] µm")
        else:
            print(f"    h = {h_i*1e6:5.0f} µm  →  "
                  f"NO reconstructable band at A = 100 nm")

    print()
    print("✓ All figures generated as PDF vector graphics.")
    print("  fig1_spatial_profiles.pdf       — spatial profiles z(x), Φ(x), Δω(x)")
    print("  fig2_transfer_function.pdf       — R(k;h) = k² e^{-kh}")
    print("  fig3_detectability_map.pdf       — SNR(k, A; h) contour maps")
    print("  fig4_reconstruction_conditions.pdf — A_min(k) and reconstructable band")
    print("  fig5_signal_vs_height.pdf        — |Δω| vs h for various λ")
