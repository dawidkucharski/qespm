"""
metrology_toolkit.py — QESPM Metrology Toolkit
================================================
Implements the metrology framework for Quantum Electrostatic Scanning
Probe Microscopy (QESPM):

  1. Instrument Transfer Function (ITF) computation and visualisation
  2. Spatial resolution analysis (Rayleigh, Sparrow, vertical)
  3. GUM-compliant uncertainty propagation
  4. ISO 25178 surface texture parameter computation from reconstructed z(x,y)
  5. Calibration factor extraction from sinusoidal standards

All computations in SI units.  Designed for publication-quality PDF figure
generation.

Author: QESPM metrology framework
"""

import numpy as np
from numpy.fft import fft2, ifft2, fftfreq
from dataclasses import dataclass, field
from typing import Optional, Dict, Tuple
import matplotlib
import matplotlib.pyplot as plt


# =====================================================================
# Physical constants
# =====================================================================
E_CHARGE = 1.602176634e-19       # C
M_CA40   = 40.0 * 1.660539e-27   # kg
EPS0     = 8.8541878128e-12      # F/m


# =====================================================================
# Instrument parameters
# =====================================================================
@dataclass
class QESPMParameters:
    """QESPM instrument parameters with nominal uncertainties."""

    # ---- Ion parameters ----
    ion_mass: float = M_CA40
    u_ion_mass: float = 1e-9 * M_CA40      # negligible (< 10⁻⁹ relative)

    # ---- Trap parameters ----
    freq_sec_x: float = 1.0e6              # Hz
    u_freq_sec_x: float = 10.0             # Hz (calibration uncertainty)
    E_bg: float = 1.0e4                    # V/m, background field at surface
    u_E_bg: float = 10.0                   # V/m (from trap calibration)

    # ---- Scan parameters ----
    ion_height: float = 40e-6              # m
    u_ion_height: float = 100e-9           # m (height calibration)
    dx_sample: float = 0.5e-6             # m, sampling interval
    fov: float = 128e-6                    # m, field of view

    # ---- Measurement precision ----
    delta_omega_noise: float = 2 * np.pi * 5.0   # rad/s (5 Hz conservative freq. resolution)
    u_delta_omega_typeA: float = 2 * np.pi * 0.3  # rad/s (repeatability for 100 ms, unchanged)

    # ---- Surface charge (dominant systematic) ----
    sigma_residual: float = 1.13e-8          # C/m², residual charge after UV cleaning

    # ---- Temperature ----
    temperature: float = 4.0               # K
    u_temperature: float = 0.1             # K

    @property
    def omega_sec_x(self) -> float:
        return 2.0 * np.pi * self.freq_sec_x

    @property
    def prefactor(self) -> float:
        """e / (2 m ω_x) — converts ∂²Φ/∂x² to Δω [m²/(V·s)]"""
        return E_CHARGE / (2.0 * self.ion_mass * self.omega_sec_x)

    @property
    def nx(self) -> int:
        """Number of x-pixels for the nominal FOV."""
        return int(self.fov / self.dx_sample)

    @property
    def u_surface_charge_equiv_height(self) -> float:
        """
        Equivalent height uncertainty from residual surface charge.
        Uses the Fourier-space formula: φ_charge(k) = σ/(2ε₀k),
        z_eq = φ_charge/E_bg.  Evaluated at k = 2×10⁵ rad/m for the
        single-point budget; the full k-dependence is used in figures.
        """
        k_eval = 2e5  # evaluation wavenumber from Table 3
        return self.sigma_residual / (2 * EPS0 * k_eval * self.E_bg)


# =====================================================================
# 1. Instrument Transfer Function
# =====================================================================

def itf_x(
    kx: np.ndarray,
    ky: np.ndarray,
    params: QESPMParameters,
) -> np.ndarray:
    """
    Instrument transfer function for the x-channel.

    H_x(k_x, k_y; h) = -(e E_bg / 2m ω_x) · k_x² · exp(-k h)

    Parameters
    ----------
    kx, ky : ndarray
        Wavevectors [rad/m].  Broadcastable.
    params : QESPMParameters

    Returns
    -------
    H : ndarray, complex (purely real, negative)
    """
    k = np.sqrt(kx**2 + ky**2)
    return -(params.prefactor * params.E_bg) * kx**2 * np.exp(-k * params.ion_height)


def itf_magnitude(kx, ky, params):
    """|H_x(k_x, k_y)|."""
    return np.abs(itf_x(kx, ky, params))


def itf_normalised(kx, ky, params):
    """Normalised ITF (MTF analogue): |H| / max(|H|)."""
    H_mag = itf_magnitude(kx, ky, params)
    return H_mag / np.max(H_mag)


def itf_optimal_wavenumber(params: QESPMParameters) -> float:
    """Wavenumber giving maximum |H| for a given height: k_opt = √2 / h."""
    return np.sqrt(2.0) / params.ion_height


def itf_cutoff_wavelengths(params: QESPMParameters) -> Tuple[float, float]:
    """
    50% MTF cut-off wavelengths.

    Returns
    -------
    lambda_min : float
        Short-wavelength cut-off (high-k side, 50% of peak).
    lambda_max : float
        Long-wavelength limit from field of view.
    """
    # High-k cut-off: solve k² exp(-k h) = 0.5 · max
    # Approximate: k_high ≈ 3.4/h
    k_high = 3.4 / params.ion_height
    lambda_min = 2 * np.pi / k_high

    # Low-k limit: field of view
    lambda_max = params.fov

    return lambda_min, lambda_max


# =====================================================================
# 2. Spatial Resolution
# =====================================================================

def rayleigh_resolution(params: QESPMParameters) -> float:
    """Rayleigh lateral resolution: δx ≈ h/2."""
    return params.ion_height / 2.0


def sparrow_resolution(params: QESPMParameters) -> float:
    """Sparrow lateral resolution: δx ≈ h/√2."""
    return params.ion_height / np.sqrt(2.0)


def vertical_resolution(
    k: float,
    params: QESPMParameters,
) -> float:
    """
    Noise-equivalent vertical resolution at spatial frequency k.

    δz(k) = (2m ω_x / e E_bg) · δω / (k² e^{-k h})
    """
    return (params.delta_omega_noise / params.prefactor) / (
        params.E_bg * k**2 * np.exp(-k * params.ion_height)
    )


def vertical_resolution_optimal(params: QESPMParameters) -> float:
    """Vertical resolution at the optimal spatial frequency k = 2/h."""
    k_opt = 2.0 / params.ion_height
    return vertical_resolution(k_opt, params)


def resolution_area_product(params: QESPMParameters) -> float:
    """Q = δz_min / δx_Rayleigh."""
    return vertical_resolution_optimal(params) / rayleigh_resolution(params)


# =====================================================================
# 3. GUM-Compliant Uncertainty Propagation
# =====================================================================

@dataclass
class UncertaintyBudget:
    """Container for an uncertainty budget evaluation."""

    sources: Dict[str, float] = field(default_factory=dict)
    sensitivities: Dict[str, float] = field(default_factory=dict)
    contributions: Dict[str, float] = field(default_factory=dict)
    u_combined: float = 0.0
    u_expanded: float = 0.0  # k = 2

    def report(self, z_nominal: float = 100e-9, h: float = 40e-6) -> str:
        """Generate a formatted uncertainty budget report."""
        lines = []
        lines.append("=" * 70)
        lines.append("  QESPM Uncertainty Budget (GUM-compliant)")
        lines.append("=" * 70)
        lines.append(f"  Nominal surface height: z = {z_nominal*1e9:.1f} nm")
        lines.append(f"  Ion height: h = {h*1e6:.0f} µm")
        lines.append("")
        lines.append(f"  {'Source':<35s} {'u(x_i)':>12s} {'|∂z/∂x_i|':>12s} "
                     f"{'u_i(z) [nm]':>12s}  {'% var':>6s}")
        lines.append("  " + "-" * 80)

        total_var = self.u_combined**2
        for name in self.sources:
            u_i = self.contributions.get(name, 0)
            pct = u_i**2 / total_var * 100 if total_var > 0 else 0
            lines.append(
                f"  {name:<35s} {self.sources[name]:>12.3e} "
                f"{self.sensitivities[name]:>12.3e} {u_i*1e9:>12.3f}  "
                f"{pct:>5.1f}%"
            )

        lines.append("  " + "-" * 80)
        lines.append(f"  Combined standard uncertainty u_c(z):  "
                     f"{self.u_combined*1e9:.3f} nm")
        lines.append(f"  Expanded uncertainty U (k=2):          "
                     f"{self.u_expanded*1e9:.3f} nm")
        lines.append("=" * 70)
        return "\n".join(lines)


def evaluate_uncertainty(
    z: float,
    k: float,
    params: QESPMParameters,
) -> UncertaintyBudget:
    """
    Evaluate the GUM uncertainty budget for a surface height measurement
    at a single spatial frequency k.

    Parameters
    ----------
    z : float
        Nominal surface height amplitude [m].
    k : float
        Surface spatial frequency [rad/m].
    params : QESPMParameters

    Returns
    -------
    UncertaintyBudget
    """
    budget = UncertaintyBudget()

    h = params.ion_height
    E = params.E_bg
    w = params.omega_sec_x
    m = params.ion_mass
    e = E_CHARGE
    dw = params.delta_omega_noise

    # ---- Source values (standard uncertainties) ----
    u_sources = {
        "Ion height h": params.u_ion_height,
        "Background field E_bg": params.u_E_bg,
        "Frequency meas. (Type A)": params.u_delta_omega_typeA,
        "Secular frequency ω_x": params.u_freq_sec_x * 2 * np.pi,
        "Ion mass m": params.u_ion_mass,
    }
    # Surface charge: k-dependent Fourier formula σ/(2ε₀kE_bg)
    u_sources["Surface charge σ (equiv. height)"] = (
        params.sigma_residual / (2 * EPS0 * k * params.E_bg))

    # ---- Absolute sensitivity coefficients |∂z/∂x_i| ----
    # Derived from z = (2m ω_x / e E_bg) · Δω · e^{+kh} / k²
    ekh = np.exp(-k * h)

    sensitivities = {
        # Charge: z_eq = σ/(2ε₀kE_bg), ∂z/∂σ = 1/(2ε₀kE_bg)
        "Surface charge σ (equiv. height)": 1.0,  # u_source is pre-converted to equiv. height
        "Ion height h": abs(z * k),                # ∂z/∂h = z·k
        "Background field E_bg": abs(z / E),       # ∂z/∂E = -z/E
        "Frequency meas. (Type A)": (2 * m * w / (e * E * k**2 * ekh)),
        "Secular frequency ω_x": abs(z / w),       # ∂z/∂ω ≈ z/ω
        "Ion mass m": abs(z / m),                  # ∂z/∂m ≈ z/m
    }

    # ---- Uncertainty contributions u_i(z) = |∂z/∂x_i| · u(x_i) ----
    contributions = {}
    for name in u_sources:
        contributions[name] = min(
            sensitivities[name] * u_sources[name],
            1e-3,  # cap at 1 mm to prevent overflow
        )

    # ---- Combined uncertainty ----
    u_c = np.sqrt(sum(min(c, 1e-3)**2 for c in contributions.values()))

    budget.sources = u_sources
    budget.sensitivities = sensitivities
    budget.contributions = contributions
    budget.u_combined = u_c
    budget.u_expanded = 2.0 * u_c

    return budget


# =====================================================================
# 4. ISO 25178 Surface Texture Parameters
# =====================================================================

def iso25178_height_parameters(
    z: np.ndarray,
    dx: float,
    dy: Optional[float] = None,
) -> Dict[str, float]:
    """
    Compute ISO 25178-2 areal height parameters from a reconstructed
    surface height map z(x,y).

    Parameters
    ----------
    z : (ny, nx) ndarray
        Surface height map [m].
    dx, dy : float
        Pixel spacing [m].

    Returns
    -------
    params : dict
        Dictionary of ISO 25178 height parameters.
    """
    if dy is None:
        dy = dx

    z_flat = z.ravel()
    z_mean = np.mean(z_flat)
    z_centred = z - z_mean

    # -- Height parameters (ISO 25178-2, §4) --
    Sa = np.mean(np.abs(z_centred))           # Arithmetic mean height
    Sq = np.std(z_flat, ddof=0)               # Root mean square height
    Ssk = np.mean(z_centred**3) / Sq**3 if Sq > 0 else 0.0  # Skewness
    Sku = np.mean(z_centred**4) / Sq**4 if Sq > 0 else 0.0  # Kurtosis
    Sp = np.max(z_flat) - z_mean              # Maximum peak height
    Sv = abs(np.min(z_flat) - z_mean)         # Maximum pit height
    Sz = Sp + Sv                              # Maximum height

    # -- Spatial parameters (ISO 25178-2, §5) --
    # Auto-correlation length Sal: distance for ACF to decay to 0.2
    # Simplified: compute radial ACF and find 0.2 crossing
    from scipy.signal import correlate2d
    acf = correlate2d(z_centred, z_centred, mode="same")
    acf /= np.max(acf)
    centre_y, centre_x = acf.shape[0] // 2, acf.shape[1] // 2

    # Radial profile
    r_max = min(centre_x, centre_y)
    radial_acf = np.zeros(r_max)
    for r in range(r_max):
        y_idx, x_idx = np.ogrid[-centre_y:acf.shape[0]-centre_y,
                                 -centre_x:acf.shape[1]-centre_x]
        mask = (np.sqrt(x_idx**2 + y_idx**2) >= r) & \
               (np.sqrt(x_idx**2 + y_idx**2) < r + 1)
        if np.any(mask):
            radial_acf[r] = np.mean(acf[mask])

    # Sal: distance where ACF drops below 0.2
    below_02 = np.where(radial_acf < 0.2)[0]
    Sal = below_02[0] * dx if len(below_02) > 0 else r_max * dx

    # -- Hybrid parameters (ISO 25178-2, §6) --
    # RMS gradient Sdq
    dzdx = np.gradient(z, dx, axis=1)
    dzdy = np.gradient(z, dy, axis=0)
    Sdq = np.sqrt(np.mean(dzdx**2 + dzdy**2))

    # Developed interfacial area ratio Sdr
    A_sampling = dx * dy
    A_surface = np.sum(np.sqrt(1 + dzdx**2 + dzdy**2)) * A_sampling
    A_nominal = z.size * A_sampling
    Sdr = (A_surface / A_nominal - 1.0) * 100  # percent

    return {
        "Sa": Sa, "Sq": Sq, "Ssk": Ssk, "Sku": Sku,
        "Sp": Sp, "Sv": Sv, "Sz": Sz,
        "Sal": Sal, "Sdq": Sdq, "Sdr": Sdr,
    }


# =====================================================================
# 5. Calibration Factor Extraction
# =====================================================================

def calibrate_from_grating(
    k_ref: float,
    A_ref: float,
    delta_omega_measured: float,
    params: QESPMParameters,
) -> Tuple[float, float]:
    """
    Extract calibration factor from a sinusoidal grating standard.

    Given a grating with known wavenumber k_ref and amplitude A_ref,
    and the measured frequency shift amplitude, determine:
      C_cal = e E_bg / (2 m ω_x)   (the calibrated prefactor)

    Parameters
    ----------
    k_ref : float
        Grating wavenumber [rad/m].
    A_ref : float
        Certified grating amplitude [m].
    delta_omega_measured : float
        Measured |Δω_x| amplitude [rad/s].
    params : QESPMParameters

    Returns
    -------
    C_cal : float
        Calibrated prefactor [m²/(V·s²)].
    u_C_cal : float
        Standard uncertainty in C_cal.
    """
    h = params.ion_height
    ekh = np.exp(-k_ref * h)

    # Expected: |Δω| = C · E_bg · A_ref · k_ref² · e^{-k_ref h}
    # But actually C already includes E_bg: C = e·E_bg/(2mω)
    # So: |Δω| = C_cal · A_ref · k_ref² · e^{-k_ref h}
    C_cal = delta_omega_measured / (A_ref * k_ref**2 * ekh)

    # Uncertainty propagation (simplified)
    u_A = 0.5e-9  # 0.5 nm standard uncertainty in A_ref (typical PTB/NIST)
    u_dw = params.u_delta_omega_typeA
    u_h = params.u_ion_height

    u_C = C_cal * np.sqrt(
        (u_A / A_ref)**2 +
        (u_dw / delta_omega_measured)**2 +
        (k_ref * u_h)**2   # from e^{-k h} sensitivity
    )

    return C_cal, u_C


# =====================================================================
# 6. Main: generate metrology figures
# =====================================================================
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

    params = QESPMParameters()

    # ==================================================================
    # Figure M1: ITF magnitude and phase
    # ==================================================================
    k_vals = np.logspace(3, 8, 500)
    heights = [5e-6, 10e-6, 40e-6, 100e-6]
    colours = plt.cm.viridis(np.linspace(0.15, 0.9, len(heights)))

    fig1, (ax1a, ax1b) = plt.subplots(1, 2, figsize=(14, 6))

    for h_i, c in zip(heights, colours):
        p = QESPMParameters(ion_height=h_i)
        H_mag = itf_magnitude(k_vals, np.zeros_like(k_vals), p)
        H_mag_norm = H_mag / np.max(H_mag) if np.max(H_mag) > 0 else H_mag

        ax1a.loglog(k_vals, H_mag, color=c, lw=1.8,
                     label=f"$h = {h_i*1e6:.0f}$ µm")
        # Mark k_opt
        k_opt_i = itf_optimal_wavenumber(p)
        ax1a.axvline(k_opt_i, color=c, ls=":", lw=1, alpha=0.5)

        ax1b.semilogx(k_vals, H_mag_norm, color=c, lw=1.8,
                       label=f"$h = {h_i*1e6:.0f}$ µm")

    # 50% MTF line
    ax1b.axhline(0.5, color="gray", ls="--", lw=1, alpha=0.5)
    ax1b.text(2e3, 0.52, "50% MTF", fontsize=9, color="gray")

    ax1a.set_xlabel("Wavenumber $k$ [rad/m]", fontsize=12)
    ax1a.set_ylabel("$|H_x(k)|$  [rad·s$^{-1}$·m$^{-1}$]", fontsize=12)
    ax1a.set_title("(a) ITF magnitude  $|H_x(k;h)|$", fontsize=13, fontweight="bold")
    ax1a.legend(fontsize=9)
    ax1a.grid(True, alpha=0.3, lw=0.5, which="both")

    ax1b.set_xlabel("Wavenumber $k$ [rad/m]", fontsize=12)
    ax1b.set_ylabel("Normalised ITF (MTF)", fontsize=12)
    ax1b.set_title("(b) Normalised ITF  $|H_x|/\\max|H_x|$", fontsize=13,
                   fontweight="bold")
    ax1b.legend(fontsize=9, loc="upper right")
    ax1b.grid(True, alpha=0.3, lw=0.5)
    ax1b.set_ylim(0, 1.05)

    # Secondary x-axis: wavelength
    for ax in [ax1a, ax1b]:
        ax_top = ax.twiny()
        ax_top.set_xscale("log")
        ax_top.set_xlim(ax.get_xlim())
        lam_ticks = np.array([0.1, 1, 10, 100, 1000]) * 1e-6
        ax_top.set_xticks([2*np.pi/lt for lt in lam_ticks])
        ax_top.set_xticklabels([f"{lt*1e6:.0f}" for lt in lam_ticks])
        ax_top.set_xlabel("$\\lambda$ [µm]", fontsize=10)

    fig1.suptitle("QESPM Instrument Transfer Function", fontsize=14,
                  fontweight="bold")
    plt.tight_layout()
    fig1.savefig("figM1_itf.pdf")
    print("Saved figM1_itf.pdf")

    # ==================================================================
    # Figure M2: Spatial resolution vs ion height
    # ==================================================================
    h_range = np.logspace(-6, -3, 200)
    rayleigh = np.array([rayleigh_resolution(QESPMParameters(ion_height=h))
                          for h in h_range])
    sparrow = np.array([sparrow_resolution(QESPMParameters(ion_height=h))
                         for h in h_range])
    vert_res = np.array([vertical_resolution_optimal(QESPMParameters(ion_height=h))
                          for h in h_range])

    fig2, (ax2a, ax2b) = plt.subplots(1, 2, figsize=(14, 6))

    ax2a.loglog(h_range * 1e6, rayleigh * 1e6, "b-", lw=2,
                label="Rayleigh: $\\delta x \\approx h/2$")
    ax2a.loglog(h_range * 1e6, sparrow * 1e6, "b--", lw=1.5,
                label="Sparrow: $\\delta x \\approx h/\\sqrt{2}$")
    ax2a.set_xlabel("Ion height $h$ [µm]", fontsize=12)
    ax2a.set_ylabel("Lateral resolution [µm]", fontsize=12)
    ax2a.set_title("(a) Lateral resolution", fontsize=13, fontweight="bold")
    ax2a.legend(fontsize=9)
    ax2a.grid(True, alpha=0.3, lw=0.5, which="both")

    # Reference: AFM, optical
    ax2a.axhline(0.001, color="gray", ls=":", alpha=0.5)
    ax2a.text(2, 0.0015, "AFM (~1 nm)", fontsize=8, color="gray")
    ax2a.axhline(0.5, color="gray", ls=":", alpha=0.5)
    ax2a.text(2, 0.6, "Optical (~0.5 µm)", fontsize=8, color="gray")

    ax2b.loglog(h_range * 1e6, vert_res * 1e9, "r-", lw=2,
                label="$\\delta z_{\\min}$ at $k_{\\rm opt}$")
    ax2b.set_xlabel("Ion height $h$ [µm]", fontsize=12)
    ax2b.set_ylabel("Vertical resolution [nm]", fontsize=12)
    ax2b.set_title("(b) Vertical resolution", fontsize=13, fontweight="bold")
    ax2b.grid(True, alpha=0.3, lw=0.5, which="both")

    # Reference
    ax2b.axhline(0.1, color="gray", ls=":", alpha=0.5)
    ax2b.text(2, 0.12, "AFM (~0.1 nm)", fontsize=8, color="gray")
    ax2b.axhline(1.0, color="gray", ls=":", alpha=0.5)
    ax2b.text(2, 1.1, "Optical (~1 nm)", fontsize=8, color="gray")

    # Power-law fit
    coeffs = np.polyfit(np.log(h_range), np.log(vert_res), 1)
    ax2b.text(0.5, 0.05, f"$\\delta z \\propto h^{{{coeffs[0]:.1f}}}$",
              transform=ax2b.transAxes, fontsize=10, color="red")

    fig2.suptitle("QESPM Spatial Resolution", fontsize=14, fontweight="bold")
    plt.tight_layout()
    fig2.savefig("figM2_resolution.pdf")
    print("Saved figM2_resolution.pdf")

    # ==================================================================
    # Figure M3: Uncertainty budget
    # ==================================================================
    # Evaluate at several (k, z) combinations
    k_test = np.logspace(3, 7, 6)
    z_test = [10e-9, 100e-9]

    fig3, axes3 = plt.subplots(1, 2, figsize=(14, 6))

    for ax, z_i in zip(axes3, z_test):
        u_combined = []
        u_charge = []
        u_height = []
        u_field = []
        u_typeA = []

        for k_j in k_test:
            budget = evaluate_uncertainty(z_i, k_j, params)
            u_combined.append(budget.u_combined)
            u_charge.append(budget.contributions.get(
                "Surface charge σ (equiv. height)", 0))
            u_height.append(budget.contributions.get("Ion height h", 0))
            u_field.append(budget.contributions.get("Background field E_bg", 0))
            u_typeA.append(budget.contributions.get("Frequency meas. (Type A)", 0))

        u_combined = np.array(u_combined)
        u_charge = np.array(u_charge)
        u_height = np.array(u_height)
        u_field = np.array(u_field)
        u_typeA = np.array(u_typeA)

        ax.loglog(k_test, u_combined * 1e9, "k-", lw=2.5, label="$u_c$ (combined)")
        ax.loglog(k_test, u_charge * 1e9, ls="--", label="Surface charge")
        ax.loglog(k_test, u_height * 1e9, ls="--", label="Ion height $h$")
        ax.loglog(k_test, u_field * 1e9, ls="--", label="Field $E_{\\rm bg}$")
        ax.loglog(k_test, u_typeA * 1e9, ls=":", label="Type A (freq. meas.)")

        ax.set_xlabel("Wavenumber $k$ [rad/m]", fontsize=12)
        ax.set_ylabel("Standard uncertainty $u(z)$ [nm]", fontsize=12)
        ax.set_title(f"$z = {z_i*1e9:.0f}$ nm", fontsize=13, fontweight="bold")
        ax.legend(fontsize=8, loc="upper left")
        ax.grid(True, alpha=0.3, lw=0.5, which="both")

    fig3.suptitle("QESPM Uncertainty Budget vs. Spatial Frequency",
                  fontsize=14, fontweight="bold")
    plt.tight_layout()
    fig3.savefig("figM3_uncertainty.pdf")
    print("Saved figM3_uncertainty.pdf")

    # ==================================================================
    # Figure M4: ISO 25178 parameter demonstration
    # ==================================================================
    # Generate a synthetic reconstructed surface
    nx, ny = 128, 128
    dx = 0.5e-6
    x = np.arange(nx) * dx
    y = np.arange(ny) * dx
    xx, yy = np.meshgrid(x, y)
    xc, yc = nx * dx / 2, ny * dx / 2

    # Test surface: plateau + roughness
    z_surf = np.zeros((ny, nx))
    # Plateau
    z_surf += 200e-9 * np.exp(-((xx - xc)**2 + (yy - yc)**2) / (2 * (15e-6)**2))
    # Sinusoidal grating
    z_surf += 50e-9 * np.sin(2 * np.pi * xx / 20e-6)
    # Random roughness
    rng = np.random.default_rng(42)
    from scipy.ndimage import gaussian_filter
    roughness = rng.normal(0, 5e-9, (ny, nx))
    roughness = gaussian_filter(roughness, sigma=1e-6 / dx)
    z_surf += roughness

    iso_params = iso25178_height_parameters(z_surf, dx)

    fig4, axes4 = plt.subplots(2, 3, figsize=(16, 10))

    # Height map
    im0 = axes4[0, 0].imshow(z_surf * 1e9, cmap="viridis", origin="lower")
    axes4[0, 0].set_title("Reconstructed $z(x,y)$ [nm]")
    plt.colorbar(im0, ax=axes4[0, 0])

    # Height histogram
    axes4[0, 1].hist(z_surf.ravel() * 1e9, bins=60, color="steelblue",
                      edgecolor="white", lw=0.3, density=True)
    axes4[0, 1].axvline(0, color="gray", ls="--")
    axes4[0, 1].set_xlabel("Height [nm]")
    axes4[0, 1].set_ylabel("Probability density")
    axes4[0, 1].set_title(f"Height distribution\n"
                          f"$S_a={iso_params['Sa']*1e9:.1f}$ nm, "
                          f"$S_q={iso_params['Sq']*1e9:.1f}$ nm")

    # Abbott-Firestone (bearing area) curve
    z_sorted = np.sort(z_surf.ravel())
    bearing = np.linspace(0, 100, len(z_sorted))
    axes4[0, 2].plot(bearing, (z_sorted - z_sorted.min()) * 1e9, "k-", lw=1.5)
    axes4[0, 2].set_xlabel("Bearing area ratio [%]")
    axes4[0, 2].set_ylabel("Height [nm]")
    axes4[0, 2].set_title(f"Abbott-Firestone curve\n"
                          f"$S_{{p}}={iso_params['Sp']*1e9:.1f}$, "
                          f"$S_{{v}}={iso_params['Sv']*1e9:.1f}$, "
                          f"$S_{{z}}={iso_params['Sz']*1e9:.1f}$ nm")
    axes4[0, 2].grid(True, alpha=0.3)

    # Surface gradient magnitude
    dzdx = np.gradient(z_surf, dx, axis=1)
    dzdy = np.gradient(z_surf, dx, axis=0)
    grad_mag = np.sqrt(dzdx**2 + dzdy**2)
    im3 = axes4[1, 0].imshow(grad_mag, cmap="inferno", origin="lower")
    axes4[1, 0].set_title(f"Surface gradient $|\\nabla z|$\n"
                          f"$S_{{dq}}={iso_params['Sdq']:.4f}$")
    plt.colorbar(im3, ax=axes4[1, 0])

    # Auto-correlation function
    from scipy.signal import correlate2d
    zc = z_surf - z_surf.mean()
    acf = correlate2d(zc, zc, mode="same")
    acf /= acf.max()
    im4 = axes4[1, 1].imshow(acf, cmap="RdBu_r", origin="lower",
                              vmin=-0.2, vmax=1.0)
    axes4[1, 1].set_title(f"Auto-correlation function\n"
                          f"$S_{{al}}={iso_params['Sal']*1e6:.1f}$ µm")
    plt.colorbar(im4, ax=axes4[1, 1])

    # Parameter summary table
    axes4[1, 2].axis("off")
    table_data = [
        ["Parameter", "Symbol", "Value", "Unit"],
        ["Arith. mean height", "$S_a$", f"{iso_params['Sa']*1e9:.1f}", "nm"],
        ["RMS height", "$S_q$", f"{iso_params['Sq']*1e9:.1f}", "nm"],
        ["Skewness", "$S_{sk}$", f"{iso_params['Ssk']:.3f}", "—"],
        ["Kurtosis", "$S_{ku}$", f"{iso_params['Sku']:.2f}", "—"],
        ["Max. peak height", "$S_p$", f"{iso_params['Sp']*1e9:.1f}", "nm"],
        ["Max. pit height", "$S_v$", f"{iso_params['Sv']*1e9:.1f}", "nm"],
        ["Max. height", "$S_z$", f"{iso_params['Sz']*1e9:.1f}", "nm"],
        ["Autocorr. length", "$S_{al}$", f"{iso_params['Sal']*1e6:.1f}", "µm"],
        ["RMS gradient", "$S_{dq}$", f"{iso_params['Sdq']:.4f}", "—"],
        ["Dev. area ratio", "$S_{dr}$", f"{iso_params['Sdr']:.2f}", "%"],
    ]
    table = axes4[1, 2].table(
        cellText=table_data[1:],
        colLabels=table_data[0],
        cellLoc="center",
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.0, 1.3)
    axes4[1, 2].set_title("ISO 25178-2 Parameters", fontsize=12, fontweight="bold",
                           pad=20)

    fig4.suptitle("QESPM — ISO 25178 Surface Texture Parameter Demonstration",
                  fontsize=14, fontweight="bold")
    plt.tight_layout()
    fig4.savefig("figM4_iso25178.pdf")
    print("Saved figM4_iso25178.pdf")

    # ==================================================================
    # Print uncertainty budget for a representative case
    # ==================================================================
    print()
    budget_report = evaluate_uncertainty(100e-9, 2e5, params)
    print(budget_report.report(100e-9, params.ion_height))
    print()
    print("✓ All metrology figures generated:")
    print("  figM1_itf.pdf           — ITF magnitude and normalised MTF")
    print("  figM2_resolution.pdf     — Lateral and vertical resolution vs. h")
    print("  figM3_uncertainty.pdf    — Uncertainty budget vs. spatial frequency")
    print("  figM4_iso25178.pdf       — ISO 25178 parameter demonstration")
