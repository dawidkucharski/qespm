"""
forward_model.py — Single trapped ion surface topography forward model
=====================================================================

Implements the forward map: surface topography z_s(x,y) + charge distribution
sigma(x,y) → ion motional observables (secular frequency shift, micromotion
amplitude, micromotion phase, heating rate).

Numerical approach:
  1. Compute surface potential Phi_s(x,y,z) from topography and charge via
     the Poisson kernel (half-space Dirichlet Green's function) using FFT
     convolution.
  2. Apply perturbation theory to compute motional observable shifts.
  3. Add noise models for realistic sensitivity estimation.

Author: Research proposal simulation — NOT for production use.
Licence: MIT
"""

import numpy as np
from numpy.fft import fft2, ifft2, fftfreq
from scipy.ndimage import sobel, laplace
from scipy.signal.windows import gaussian
from dataclasses import dataclass
from typing import Optional, Tuple


# ---------------------------------------------------------------------------
# Physical constants and ion parameters
# ---------------------------------------------------------------------------

@dataclass
class IonParameters:
    """Physical parameters for a trapped ion."""
    species: str = "Ca40"
    mass_amu: float = 40.0       # atomic mass units
    charge_e: float = 1.0        # in units of elementary charge
    freq_sec_x: float = 1e6      # secular frequency, Hz
    freq_sec_y: float = 1e6
    freq_sec_z: float = 1e6
    omega_rf: float = 40e6       # rf drive frequency, rad/s (2π × 20 MHz)

    @property
    def mass_kg(self) -> float:
        return self.mass_amu * 1.660539e-27

    @property
    def omega_sec_x(self) -> float:
        return 2.0 * np.pi * self.freq_sec_x

    @property
    def omega_sec_y(self) -> float:
        return 2.0 * np.pi * self.freq_sec_y

    @property
    def omega_sec_z(self) -> float:
        return 2.0 * np.pi * self.freq_sec_z


@dataclass
class ScanParameters:
    """Scan geometry parameters."""
    nx: int = 256                 # grid points in x
    ny: int = 256                 # grid points in y
    dx: float = 0.5e-6            # pixel size, m (0.5 µm)
    ion_height: float = 40e-6     # ion height above surface, m (40 µm — comparable to
                                  #   Maiwald et al., Nat. Phys. 2009)
    surface_potential_rms: float = 0.1  # RMS surface potential, V


# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

EPS0 = 8.8541878128e-12          # vacuum permittivity, F/m
E_CHARGE = 1.602176634e-19       # elementary charge, C
KB = 1.380649e-23                # Boltzmann constant, J/K
HBAR = 1.054571817e-34           # reduced Planck constant, J·s


# ---------------------------------------------------------------------------
# Forward model: topography/charge → surface potential
# ---------------------------------------------------------------------------

def poisson_kernel_halfspace(
    zs: np.ndarray,
    sigma: np.ndarray,
    height: float,
    dx: float,
    dy: Optional[float] = None,
) -> np.ndarray:
    """
    Compute the electrostatic potential at height `height` above a surface.

    Uses the Poisson kernel (half-space Dirichlet Green's function) via FFT
    convolution.  The surface boundary condition is:
        Phi(x, y, z=0) = phi_s(x, y)
    where phi_s(x, y) encodes topographic and charge contributions.

    For a conducting surface with topography zs(x,y) and surface charge sigma(x,y):
        phi_s(x, y) ≈ phi_topography + phi_charge
    where phi_topography arises from the conductor boundary deformation and
    phi_charge = sigma(x,y) / (2*eps0) * (Coulomb kernel).

    Parameters
    ----------
    zs : (ny, nx) ndarray
        Surface height profile [m].
    sigma : (ny, nx) ndarray
        Surface charge density [C/m²].
    height : float
        Ion height above the reference plane [m].
    dx : float
        Pixel size in x [m].
    dy : float, optional
        Pixel size in y [m]; defaults to dx.

    Returns
    -------
    phi : (ny, nx) ndarray
        Electrostatic potential at the ion height [V].
    """
    if dy is None:
        dy = dx

    ny, nx = zs.shape

    # --- Topographic contribution: phi_topography ≈ E_background * zs ---
    # For a grounded conductor with an applied uniform background field
    # E_bg (e.g. from trap electrodes), the surface perturbation potential
    # is proportional to the height deviation.
    # We assume an effective background field of ~1 V/m at the ion height.
    # --- Topographic contribution ---
    # For a grounded conducting surface with topography zs, the surface
    # potential arises from displacement of the equipotential surface in
    # the presence of a background field E_bg from the trap electrodes.
    #     φ_topo(x,y) ≈ E_bg · zs(x,y)
    # For a surface trap at h = 40 µm, the rf pseudopotential gradient is
    # of order 10⁴–10⁵ V/m at the ion position; at the surface the field
    # is stronger.  We use ~10⁴ V/m as a conservative estimate.
    E_bg = 1e4  # V/m
    phi_topo = E_bg * zs

    # --- Charge/patch-potential contribution ---
    # Surface charge σ [C/m²] produces surface potential φ ≈ σ·d_eff/ε₀.
    # d_eff ~ 1 nm is the effective capacitive screening length.
    # For σ = 1 mC/m², φ ≈ 1e-3 × 1e-9 / 8.85e-12 ≈ 0.11 V.
    d_eff = 1e-9  # m
    phi_charge = sigma * d_eff / EPS0

    # --- Combined surface potential ---
    phi_surface = phi_topo + phi_charge

    # --- Propagate to ion height via Poisson kernel (FFT convolution) ---
    # The Poisson kernel in Fourier space: G(k; h) = exp(-k * h)
    kx = 2.0 * np.pi * fftfreq(nx, dx)
    ky = 2.0 * np.pi * fftfreq(ny, dy)
    kx_grid, ky_grid = np.meshgrid(kx, ky)
    k = np.sqrt(kx_grid**2 + ky_grid**2)

    # Propagator
    G_k = np.exp(-k * height)
    G_k[0, 0] = 0.0  # DC component: zero net potential at infinity

    phi_hat = fft2(phi_surface)
    phi_hat *= G_k
    phi_at_ion = np.real(ifft2(phi_hat))

    return phi_at_ion


# ---------------------------------------------------------------------------
# Forward model: surface potential → motional observables
# ---------------------------------------------------------------------------

def compute_observables(
    phi: np.ndarray,
    ion: IonParameters,
    dx: float,
    dy: Optional[float] = None,
    temperature: float = 4.0,
    surface_resistivity: float = 2.44e-8,  # gold, Ω·m
) -> dict:
    """
    Convert surface potential at ion height to motional observables.

    Parameters
    ----------
    phi : (ny, nx) ndarray
        Electrostatic potential at ion position [V].
    ion : IonParameters
        Ion physical parameters.
    dx, dy : float
        Grid spacing [m].
    temperature : float
        Surface temperature [K].
    surface_resistivity : float
        Electrical resistivity of surface material [Ω·m].

    Returns
    -------
    observables : dict
        Dictionary containing:
        - 'delta_omega_x', 'delta_omega_y': secular frequency shifts [rad/s]
        - 'beta_x', 'beta_y': micromotion modulation indices
        - 'phi_mu': micromotion phase [rad]
        - 'heating_rate': motional heating rate [quanta/s]
    """
    if dy is None:
        dy = dx

    m = ion.mass_kg
    e = E_CHARGE
    height = None  # would need to be stored; we use phi directly

    # --- Secular frequency shifts (second derivatives) ---
    # Eq. (12) of proposal: Δω_u = (e / 2m ω_u) * ∂²Φ_s/∂u²
    # Use FFT for exact second derivatives (avoids FD numerical precision
    # issues for the small potentials at ion height).
    ny, nx = phi.shape
    kx = 2.0 * np.pi * fftfreq(nx, dx)
    ky = 2.0 * np.pi * fftfreq(ny, dy)
    kx_grid, ky_grid = np.meshgrid(kx, ky)

    phi_hat = fft2(phi)
    # ∂²Φ/∂x² → -k_x²·Φ̃  in Fourier space
    d2phi_dx2 = np.real(ifft2(-kx_grid**2 * phi_hat))
    d2phi_dy2 = np.real(ifft2(-ky_grid**2 * phi_hat))

    delta_omega_x = (e / (2.0 * m * ion.omega_sec_x)) * d2phi_dx2
    delta_omega_y = (e / (2.0 * m * ion.omega_sec_y)) * d2phi_dy2

    # --- Micromotion amplitude (first derivatives) ---
    # β_u ∝ |∂Φ/∂u|.  Compute via FFT for consistency.
    # ∂Φ/∂x → i·k_x·Φ̃
    dphi_dx = np.real(ifft2(1j * kx_grid * phi_hat))
    dphi_dy = np.real(ifft2(1j * ky_grid * phi_hat))

    # Assume secular amplitude ~ 100 nm (Lamb-Dicke regime)
    secular_amplitude = 100e-9  # m
    # Mathieu q-parameter (typical)
    q_x = 0.3
    q_y = 0.3

    beta_x = (q_x * e / (2.0 * m * ion.omega_sec_x**2 * secular_amplitude)) * np.abs(dphi_dx)
    beta_y = (q_y * e / (2.0 * m * ion.omega_sec_y**2 * secular_amplitude)) * np.abs(dphi_dy)

    # --- Micromotion phase ---
    phi_mu = np.arctan2(dphi_dy, dphi_dx)

    # --- Heating rate ---
    # Johnson noise contribution: S_E = k_B T ρ / (π ε₀ ω h³)
    # h is height — would vary with position for a real topographical scan
    h_ref = 10e-6  # reference ion height (should match scan parameter)
    S_E = (KB * temperature * surface_resistivity) / (
        np.pi * EPS0 * ion.omega_sec_x * h_ref**3
    )
    heating_rate = (e**2 / (4.0 * m * HBAR * ion.omega_sec_x)) * S_E
    heating_rate_map = np.full_like(phi, heating_rate)

    return {
        'delta_omega_x': delta_omega_x,
        'delta_omega_y': delta_omega_y,
        'beta_x': beta_x,
        'beta_y': beta_y,
        'phi_mu': phi_mu,
        'heating_rate': heating_rate_map,
        'dphi_dx': dphi_dx,
        'dphi_dy': dphi_dy,
        'd2phi_dx2': d2phi_dx2,
        'd2phi_dy2': d2phi_dy2,
    }


# ---------------------------------------------------------------------------
# Test surface generation
# ---------------------------------------------------------------------------

def generate_test_surface(
    nx: int,
    ny: int,
    dx: float,
    dy: Optional[float] = None,
    mode: str = "mixed",
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate synthetic calibration surface with controlled test patterns.

    Modes:
      'topo_only'  — pure topography, zero charge (validates topo → Δω)
      'charge_only' — pure charge patches, zero topography
      'mixed'       — both contributions (demonstrates the degeneracy)

    Topography features:
      - Central plateau (200 nm height, 30 µm Gaussian FWHM)
      - Linear ramp across the field of view (100 nm over 100 µm)
      - Fine sinusoidal grating (5 nm amplitude, 4 µm period)

    Charge features:
      - Localised Gaussian patches, ~8–12 µm FWHM
      - Surface potentials of ~50–250 mV per patch

    Returns
    -------
    zs : (ny, nx) ndarray, surface height [m]
    sigma : (ny, nx) ndarray, surface charge density [C/m²]
    """
    if dy is None:
        dy = dx

    rng = np.random.default_rng(42)

    x = np.arange(nx) * dx
    y = np.arange(ny) * dy
    xx, yy = np.meshgrid(x, y)

    xc = nx * dx / 2.0
    yc = ny * dy / 2.0
    fov_x = nx * dx

    # ==================================================================
    # Topography
    # ==================================================================
    if mode in ("topo_only", "mixed"):
        zs = np.zeros((ny, nx))

        # Central plateau: 200 nm height, Gaussian with σ = 7.5 µm
        plateau_height = 200e-9
        plateau_sigma = 7.5e-6
        plateau = plateau_height * np.exp(
            -((xx - xc) ** 2 + (yy - yc) ** 2) / (2 * plateau_sigma**2)
        )
        zs += plateau

        # Linear ramp: 100 nm height change across full FOV
        ramp_height = 100e-9
        ramp = ramp_height * (xx - xc) / (fov_x / 2.0)
        zs += ramp

        # Fine grating: 5 nm amplitude, 4 µm period (sub-resolution test)
        grating_period = 4e-6
        grating_amplitude = 5e-9
        zs += grating_amplitude * np.sin(2.0 * np.pi * (xx + yy) / grating_period)
    else:
        zs = np.zeros((ny, nx))

    # ==================================================================
    # Surface charge patches
    # ==================================================================
    sigma = np.zeros((ny, nx))

    if mode in ("charge_only", "mixed"):
        # For d_eff = 1 nm: σ = 1 mC/m² → φ_surface ≈ 0.11 V
        patch_specs = [
            # (cx, cy, sigma_FWHM [m], charge_density [C/m²])
            (xc - 25e-6, yc, 8e-6, 2.0e-3),          # ~0.22 V
            (xc + 30e-6, yc - 15e-6, 10e-6, -1.5e-3), # ~−0.17 V
            (xc + 15e-6, yc + 20e-6, 6e-6, 1.0e-3),   # ~0.11 V
            (xc - 10e-6, yc - 25e-6, 12e-6, -2.5e-3), # ~−0.28 V
        ]

        for cx, cy, fwhm, charge_dens in patch_specs:
            patch_sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
            gauss = np.exp(
                -((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * patch_sigma**2)
            )
            sigma += charge_dens * gauss

    return zs, sigma


# ---------------------------------------------------------------------------
# Noise models
# ---------------------------------------------------------------------------

def add_measurement_noise(
    observables: dict,
    ion: IonParameters,
    integration_time: float = 1.0,
    rng: Optional[np.random.Generator] = None,
) -> dict:
    """
    Add realistic measurement noise to observables.

    Parameters
    ----------
    observables : dict
        Noise-free observables from compute_observables().
    ion : IonParameters
    integration_time : float
        Per-pixel integration time [s].
    rng : np.random.Generator, optional

    Returns
    -------
    noisy : dict, same keys as observables with added noise.
    """
    if rng is None:
        rng = np.random.default_rng()

    ny, nx = observables['delta_omega_x'].shape
    noisy = {}

    # Frequency shift noise: δω ~ 2π × 1 Hz (typical for 1 ms integration)
    # Scales as 1/√τ
    delta_omega_noise_1ms = 2.0 * np.pi * 1.0  # rad/s for 1 ms
    delta_omega_noise = delta_omega_noise_1ms / np.sqrt(integration_time / 1e-3)
    noisy['delta_omega_x'] = observables['delta_omega_x'] + rng.normal(
        0, delta_omega_noise, (ny, nx)
    )
    noisy['delta_omega_y'] = observables['delta_omega_y'] + rng.normal(
        0, delta_omega_noise, (ny, nx)
    )

    # Micromotion modulation index noise: δβ ~ 10⁻⁴
    beta_noise = 1e-4 / np.sqrt(integration_time / 1e-3)
    noisy['beta_x'] = np.abs(observables['beta_x'] + rng.normal(0, beta_noise, (ny, nx)))
    noisy['beta_y'] = np.abs(observables['beta_y'] + rng.normal(0, beta_noise, (ny, nx)))

    # Phase noise: δφ ~ 0.01 rad
    phi_noise = 0.01 / np.sqrt(integration_time / 1e-3)
    noisy['phi_mu'] = observables['phi_mu'] + rng.normal(0, phi_noise, (ny, nx))

    # Heating rate noise: dominated by Poisson statistics of motional quanta
    heating_rate = observables['heating_rate']
    heating_noise = np.sqrt(np.abs(heating_rate) * integration_time) / integration_time
    noisy['heating_rate'] = heating_rate + rng.normal(0, heating_noise, (ny, nx))

    # Keep noise-free derivatives for debugging
    noisy['dphi_dx'] = observables['dphi_dx']
    noisy['dphi_dy'] = observables['dphi_dy']
    noisy['d2phi_dx2'] = observables['d2phi_dx2']
    noisy['d2phi_dy2'] = observables['d2phi_dy2']

    return noisy


# ---------------------------------------------------------------------------
# Main simulation
# ---------------------------------------------------------------------------

def run_forward_simulation(
    ion: Optional[IonParameters] = None,
    scan: Optional[ScanParameters] = None,
    mode: str = "mixed",
) -> dict:
    """
    Run the complete forward simulation: topography → observables.

    Parameters
    ----------
    mode : str
        'topo_only', 'charge_only', or 'mixed'.

    Returns
    -------
    dict with keys: zs, sigma, phi, observables, observables_noisy,
                    ion, scan, mode
    """
    if ion is None:
        ion = IonParameters()
    if scan is None:
        scan = ScanParameters()

    print(f"=== Forward Model Simulation [{mode}] ===")
    print(f"Ion:  {ion.species}, m = {ion.mass_amu} amu")
    print(f"Secular frequencies:  ω_x/2π = {ion.freq_sec_x/1e6:.1f} MHz")
    print(f"                       ω_y/2π = {ion.freq_sec_y/1e6:.1f} MHz")
    print(f"Grid: {scan.nx}×{scan.ny}, dx = {scan.dx*1e6:.2f} µm")
    print(f"Ion height: {scan.ion_height*1e6:.1f} µm")
    print(f"FOV: {scan.nx*scan.dx*1e6:.0f} × {scan.ny*scan.dx*1e6:.0f} µm")
    print()

    # 1. Generate test surface
    print("Generating test surface...")
    zs, sigma = generate_test_surface(scan.nx, scan.ny, scan.dx, mode=mode)
    print(f"  Topography range: [{zs.min()*1e9:.1f}, {zs.max()*1e9:.1f}] nm")
    print(f"  Charge range:     [{sigma.min()*1e6:.2f}, {sigma.max()*1e6:.2f}] µC/m²")

    # 2. Compute surface potential at ion height
    print("Computing surface potential at ion height...")
    phi = poisson_kernel_halfspace(zs, sigma, scan.ion_height, scan.dx)
    print(f"  Φ at ion height:  [{phi.min()*1e6:.2f}, {phi.max()*1e6:.2f}] µV")

    # 3. Compute motional observables
    print("Computing motional observables...")
    obs = compute_observables(phi, ion, scan.dx)

    # 4. Add measurement noise
    print("Adding measurement noise...")
    obs_noisy = add_measurement_noise(obs, ion, integration_time=0.1)

    # Report
    omega_x = ion.omega_sec_x
    dw_max = max(abs(obs['delta_omega_x'].min()), abs(obs['delta_omega_x'].max()))
    print()
    print("--- Observable ranges (noise-free) ---")
    print(f"  Δω_x:  [{obs['delta_omega_x'].min():.2f}, "
          f"{obs['delta_omega_x'].max():.2f}] rad/s")
    print(f"         = [{obs['delta_omega_x'].min()/(2*np.pi):.1f}, "
          f"{obs['delta_omega_x'].max()/(2*np.pi):.1f}] Hz")
    print(f"  |Δω_x| / ω_x = {dw_max/omega_x*100:.2f}%"
          f"  {'✓ (perturbation valid)' if dw_max/omega_x < 0.1 else '⚠ (non-perturbative)'}")
    print(f"  β_x:   [{obs['beta_x'].min():.2e}, {obs['beta_x'].max():.2e}]")
    print(f"  ḣ:     {obs['heating_rate'][0,0]:.1f} quanta/s (uniform)")

    return {
        'zs': zs,
        'sigma': sigma,
        'phi': phi,
        'observables': obs,
        'observables_noisy': obs_noisy,
        'ion': ion,
        'scan': scan,
        'mode': mode,
    }


# ---------------------------------------------------------------------------
# Simple Tikhonov inversion (demonstration)
# ---------------------------------------------------------------------------

def tikhonov_inversion(
    delta_omega: np.ndarray,
    dx: float,
    omega_sec: float,
    ion_mass: float,
    lambda_reg: Optional[float] = None,
    dy: Optional[float] = None,
) -> np.ndarray:
    """
    Tikhonov-regularised inversion of Δω_x → Φ_at_ion in Fourier space.

    Forward model (from potential at ion height to x-frequency shift):
        Δω̃_x(k) = -(e / 2mω) · k_x² · Φ̃_at_ion(k)

    Note: This reconstructs the potential AT THE ION HEIGHT, not at the
    surface.  The e^{-k·h} propagation from surface to ion height is a
    separate, subsequent step if surface potential/topography is desired.

    Tikhonov-regularised inverse:
        Φ̃_λ(k) = H*(k) · Δω̃_x(k) / (|H(k)|² + λ²)

    where H(k) = -(e/2mω) · k_x².

    λ is chosen automatically via a heuristic based on max(|H|) if not
    provided.

    Parameters
    ----------
    delta_omega : (ny, nx) ndarray
        Measured secular frequency shift [rad/s].
    dx : float
        Pixel size [m].
    omega_sec : float
        Secular frequency [rad/s].
    ion_mass : float
        Ion mass [kg].
    lambda_reg : float, optional
        Tikhonov regularisation parameter [same units as H].
        If None, uses λ = 1e-2 × max(|H|) as a reasonable default.
    dy : float, optional

    Returns
    -------
    phi_recon : (ny, nx) ndarray
        Reconstructed electrostatic potential at the ion height [V].
    """
    if dy is None:
        dy = dx

    ny, nx = delta_omega.shape
    e = E_CHARGE

    # Fourier coordinates
    kx = 2.0 * np.pi * fftfreq(nx, dx)
    ky = 2.0 * np.pi * fftfreq(ny, dy)
    kx_grid, ky_grid = np.meshgrid(kx, ky)

    # Forward operator (no e^{-k·h} — operates on Φ at ion height)
    # H(k) = -(e / 2mω) * k_x²
    H = -(e / (2.0 * ion_mass * omega_sec)) * kx_grid**2
    H[0, 0] = 0.0  # DC mode not recoverable

    # Fourier transform of data
    dw_hat = fft2(delta_omega)

    # Choose λ if not provided
    if lambda_reg is None:
        # Heuristic: λ should be ~noise_level in Δω space, mapped to H space.
        # For a sensible default, use λ = 10⁻⁴ × max(|H|), which empirically
        # gives good results for our test cases (balances signal recovery
        # against noise amplification from the 1/k_x² divergence).
        H_mag = np.abs(H)
        nonzero = H_mag[H_mag > 0]
        if len(nonzero) > 0:
            lambda_reg = 1e-4 * np.max(nonzero)
        else:
            lambda_reg = 1e-10

    # Tikhonov filter: H* / (|H|² + λ²)
    H_sq = np.abs(H) ** 2
    denom = H_sq + lambda_reg**2
    inv_filter = np.zeros_like(H, dtype=complex)
    mask = denom > 1e-30
    inv_filter[mask] = np.conj(H[mask]) / denom[mask]

    # Apply filter
    phi_hat = dw_hat * inv_filter
    phi_hat[0, 0] = 0.0
    phi_recon = np.real(ifft2(phi_hat))

    return phi_recon


# ---------------------------------------------------------------------------
# Run if executed as script
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")  # non-interactive backend
    import matplotlib.pyplot as plt

    # ==================================================================
    # Run simulations for all three modes
    # ==================================================================
    all_results = {}
    for mode in ["topo_only", "charge_only", "mixed"]:
        all_results[mode] = run_forward_simulation(mode=mode)
        print()

    # ==================================================================
    # Multi-panel figure: forward model overview
    # ==================================================================
    fig, axes = plt.subplots(3, 5, figsize=(22, 13))

    mode_labels = {
        "topo_only": "Topography only",
        "charge_only": "Charge only",
        "mixed": "Mixed (topo + charge)",
    }

    for row, mode in enumerate(["topo_only", "charge_only", "mixed"]):
        r = all_results[mode]
        zs = r['zs']
        sigma = r['sigma']
        phi = r['phi']
        obs = r['observables']
        obs_n = r['observables_noisy']

        # Column 0: topography
        im00 = axes[row, 0].imshow(zs * 1e9, cmap='viridis', origin='lower')
        axes[row, 0].set_title(f'{mode_labels[mode]}\n$z_s$ [nm]', fontsize=9)
        plt.colorbar(im00, ax=axes[row, 0])

        # Column 1: charge
        im01 = axes[row, 1].imshow(sigma * 1e6, cmap='RdBu_r', origin='lower')
        axes[row, 1].set_title('$\\sigma$ [$\\mu$C/m²]', fontsize=9)
        plt.colorbar(im01, ax=axes[row, 1])

        # Column 2: potential at ion height
        im02 = axes[row, 2].imshow(phi * 1e6, cmap='RdBu_r', origin='lower')
        axes[row, 2].set_title('$\\Phi_s$ at ion [$\\mu$V]', fontsize=9)
        plt.colorbar(im02, ax=axes[row, 2])

        # Column 3: Δω_x (noisy)
        im03 = axes[row, 3].imshow(
            obs_n['delta_omega_x'] / (2 * np.pi), cmap='coolwarm', origin='lower'
        )
        axes[row, 3].set_title('$\\Delta\\omega_x/2\\pi$ [Hz]', fontsize=9)
        plt.colorbar(im03, ax=axes[row, 3])

        # Column 4: micromotion phase
        im04 = axes[row, 4].imshow(
            obs['phi_mu'], cmap='twilight', origin='lower',
            vmin=-np.pi, vmax=np.pi,
        )
        axes[row, 4].set_title('$\\phi_\\mu$ [rad]', fontsize=9)
        plt.colorbar(im04, ax=axes[row, 4])

    fig.suptitle(
        "Single-Ion Surface Topography Forward Model\n"
        "$^{40}$Ca$^+$, $h = 40$ µm, $\\Delta x = 0.5$ µm, FOV $128 \\times 128$ µm²",
        fontsize=14, fontweight='bold',
    )
    plt.tight_layout()
    plt.savefig("forward_model_output.pdf", dpi=150, bbox_inches='tight')
    print("Figure saved to forward_model_output.pdf")

    # ==================================================================
    # Manuscript Figure: forward-model chain for the topography-only case
    # ==================================================================
    r = all_results["topo_only"]
    zs = r['zs']
    phi = r['phi']
    obs = r['observables']
    dw_hz = obs['delta_omega_x'] / (2 * np.pi)

    fig1b, axes1b = plt.subplots(1, 3, figsize=(18, 5.5))

    im_a = axes1b[0].imshow(zs * 1e9, cmap='terrain', origin='lower',
                            extent=[0, 128, 0, 128])
    axes1b[0].set_title('(a) Surface height $z_s(x,y)$ [nm]', fontsize=12)
    plt.colorbar(im_a, ax=axes1b[0], fraction=0.046)
    axes1b[0].set_xlabel('$x$ [µm]', fontsize=11)
    axes1b[0].set_ylabel('$y$ [µm]', fontsize=11)

    vmax_p = np.percentile(np.abs(phi), 99)
    im_b = axes1b[1].imshow(phi * 1e6, cmap='RdBu_r', origin='lower',
                            vmin=-vmax_p * 1e6, vmax=vmax_p * 1e6,
                            extent=[0, 128, 0, 128])
    axes1b[1].set_title('(b) Potential at ion height $\\Phi(x,y,h)$ [µV]',
                        fontsize=12)
    plt.colorbar(im_b, ax=axes1b[1], fraction=0.046)
    axes1b[1].set_xlabel('$x$ [µm]', fontsize=11)
    axes1b[1].set_ylabel('$y$ [µm]', fontsize=11)

    vmax_d = np.percentile(np.abs(dw_hz), 99)
    im_c = axes1b[2].imshow(dw_hz, cmap='coolwarm', origin='lower',
                            vmin=-vmax_d, vmax=vmax_d,
                            extent=[0, 128, 0, 128])
    axes1b[2].set_title('(c) Secular shift $\\Delta\\omega_x/2\\pi$ [Hz]',
                        fontsize=12)
    plt.colorbar(im_c, ax=axes1b[2], fraction=0.046)
    axes1b[2].set_xlabel('$x$ [µm]', fontsize=11)
    axes1b[2].set_ylabel('$y$ [µm]', fontsize=11)

    fig1b.suptitle(
        "Forward-Model Chain (topography only)\n"
        "$^{40}$Ca$^+$, $f_{\\rm sec}=1$ MHz, $h=40$ µm, "
        "$\\Delta x = 0.5$ µm, FOV $128 \\times 128$ µm$^2$",
        fontsize=14, fontweight='bold',
    )
    plt.tight_layout()
    fig1b.savefig("fig1b_spatial.pdf", dpi=150, bbox_inches='tight')
    print("Figure saved to fig1b_spatial.pdf")

    # ==================================================================
    # Tikhonov inversion comparison (topo_only vs mixed)
    # ==================================================================
    fig2, axes2 = plt.subplots(2, 4, figsize=(20, 10))

    for row, mode in enumerate(["topo_only", "mixed"]):
        r = all_results[mode]
        phi_true = r['phi']
        obs_n = r['observables_noisy']

        # Invert
        phi_recon = tikhonov_inversion(
            obs_n['delta_omega_x'],
            dx=r['scan'].dx,
            omega_sec=r['ion'].omega_sec_x,
            ion_mass=r['ion'].mass_kg,
        )

        # Metrics
        phi_t_zm = phi_true - phi_true.mean()
        phi_r_zm = phi_recon - phi_recon.mean()
        corr = np.corrcoef(phi_t_zm.ravel(), phi_r_zm.ravel())[0, 1]
        mse = np.mean((phi_true - phi_recon) ** 2)
        signal_range = phi_true.max() - phi_true.min()

        # Plot
        im0 = axes2[row, 0].imshow(phi_true * 1e6, cmap='RdBu_r', origin='lower')
        axes2[row, 0].set_title(f'{mode_labels[mode]}: True $\\Phi_s$ [$\\mu$V]', fontsize=9)
        plt.colorbar(im0, ax=axes2[row, 0])

        im1 = axes2[row, 1].imshow(phi_recon * 1e6, cmap='RdBu_r', origin='lower')
        axes2[row, 1].set_title(
            f'Tikhonov recon. [$\\mu$V]\ncorr = {corr:.3f}, '
            f'RMSE/Δ = {np.sqrt(mse)/(signal_range+1e-30)*100:.0f}%',
            fontsize=9,
        )
        plt.colorbar(im1, ax=axes2[row, 1])

        im2 = axes2[row, 2].imshow((phi_true - phi_recon) * 1e6, cmap='RdBu_r', origin='lower')
        axes2[row, 2].set_title('Residual [$\\mu$V]', fontsize=9)
        plt.colorbar(im2, ax=axes2[row, 2])

        # Cross-section through centre
        mid_y = phi_true.shape[0] // 2
        axes2[row, 3].plot(
            np.arange(phi_true.shape[1]) * r['scan'].dx * 1e6,
            phi_true[mid_y, :] * 1e6, 'k-', label='True', lw=1.5,
        )
        axes2[row, 3].plot(
            np.arange(phi_true.shape[1]) * r['scan'].dx * 1e6,
            phi_recon[mid_y, :] * 1e6, 'r--', label='Recon.', lw=1.5,
        )
        axes2[row, 3].set_xlabel('x [µm]', fontsize=9)
        axes2[row, 3].set_ylabel('$\\Phi_s$ [$\\mu$V]', fontsize=9)
        axes2[row, 3].set_title(f'Cross-section (y = centre)', fontsize=9)
        axes2[row, 3].legend(fontsize=8)
        axes2[row, 3].grid(True, alpha=0.3)

    fig2.suptitle(
        "Tikhonov-Regularised Inversion: $\\Delta\\omega_x \\to \\Phi_s$\n"
        "Topography-only vs. Mixed (topo + charge) — the degeneracy problem",
        fontsize=13, fontweight='bold',
    )
    plt.tight_layout()
    plt.savefig("tikhonov_inversion_demo.pdf", dpi=150, bbox_inches='tight')
    print("Figure saved to tikhonov_inversion_demo.pdf")

    # ==================================================================
    # Manuscript Figure: 4-panel inversion for the topography-only case
    # ==================================================================
    r_top = all_results["topo_only"]
    phi_true = r_top['phi']
    phi_recon = tikhonov_inversion(
        r_top['observables_noisy']['delta_omega_x'],
        dx=r_top['scan'].dx,
        omega_sec=r_top['ion'].omega_sec_x,
        ion_mass=r_top['ion'].mass_kg,
    )
    corr_inv = np.corrcoef((phi_true - phi_true.mean()).ravel(),
                           (phi_recon - phi_recon.mean()).ravel())[0, 1]
    resid = phi_true - phi_recon
    vmax_pt = np.percentile(np.abs(phi_true), 99)

    fig_inv, axes_inv = plt.subplots(1, 4, figsize=(20, 5))

    im0 = axes_inv[0].imshow(phi_true * 1e6, cmap='RdBu_r', origin='lower',
                             vmin=-vmax_pt * 1e6, vmax=vmax_pt * 1e6,
                             extent=[0, 128, 0, 128])
    axes_inv[0].set_title('(a) True potential $\\Phi(x,y,h)$ [µV]', fontsize=11)
    plt.colorbar(im0, ax=axes_inv[0], fraction=0.046)
    axes_inv[0].set_xlabel('$x$ [µm]', fontsize=10)
    axes_inv[0].set_ylabel('$y$ [µm]', fontsize=10)

    im1 = axes_inv[1].imshow(phi_recon * 1e6, cmap='RdBu_r', origin='lower',
                             vmin=-vmax_pt * 1e6, vmax=vmax_pt * 1e6,
                             extent=[0, 128, 0, 128])
    axes_inv[1].set_title(f'(b) Tikhonov reconstruction\n$r = {corr_inv:.3f}$',
                          fontsize=11)
    plt.colorbar(im1, ax=axes_inv[1], fraction=0.046)
    axes_inv[1].set_xlabel('$x$ [µm]', fontsize=10)
    axes_inv[1].set_ylabel('$y$ [µm]', fontsize=10)

    im2 = axes_inv[2].imshow(resid * 1e6, cmap='RdBu_r', origin='lower',
                             extent=[0, 128, 0, 128])
    axes_inv[2].set_title('(c) Residual $\\Phi - \\Phi_{\\rm rec}$ [µV]',
                          fontsize=11)
    plt.colorbar(im2, ax=axes_inv[2], fraction=0.046)
    axes_inv[2].set_xlabel('$x$ [µm]', fontsize=10)
    axes_inv[2].set_ylabel('$y$ [µm]', fontsize=10)

    mid_y = phi_true.shape[0] // 2
    x_prof = np.arange(phi_true.shape[1]) * r_top['scan'].dx * 1e6
    axes_inv[3].plot(x_prof, phi_true[mid_y, :] * 1e6, 'k-', lw=1.5,
                     label='True')
    axes_inv[3].plot(x_prof, phi_recon[mid_y, :] * 1e6, 'r--', lw=1.5,
                     label='Recon.')
    axes_inv[3].set_xlabel('$x$ [µm]', fontsize=10)
    axes_inv[3].set_ylabel('$\\Phi$ [µV]', fontsize=10)
    axes_inv[3].set_title('(d) Cross-section through centre', fontsize=11)
    axes_inv[3].legend(fontsize=9)
    axes_inv[3].grid(True, alpha=0.3)

    fig_inv.suptitle(
        "Tikhonov-Regularised Inversion: $\\Delta\\omega_x \\to \\Phi$ (topography only)\n"
        "$^{40}$Ca$^+$, $f_{\\rm sec}=1$ MHz, $h=40$ µm, "
        "$\\Delta x = 0.5$ µm",
        fontsize=13, fontweight='bold',
    )
    plt.tight_layout()
    fig_inv.savefig("fig2_inversion.pdf", dpi=150, bbox_inches='tight')
    print("Figure saved to fig2_inversion.pdf")

    # ==================================================================
    # Resolution analysis: vary ion height
    # ==================================================================
    print("\n--- Resolution vs. Ion Height ---")
    heights = np.array([5, 10, 20, 40, 80]) * 1e-6  # metres
    correlations = []
    rmses = []

    for h in heights:
        scan_h = ScanParameters(ion_height=h)
        r_h = run_forward_simulation(mode="topo_only", scan=scan_h)

        phi_recon = tikhonov_inversion(
            r_h['observables_noisy']['delta_omega_x'],
            dx=r_h['scan'].dx,
            omega_sec=r_h['ion'].omega_sec_x,
            ion_mass=r_h['ion'].mass_kg,
        )

        phi_t = r_h['phi']
        corr = np.corrcoef(
            (phi_t - phi_t.mean()).ravel(),
            (phi_recon - phi_recon.mean()).ravel(),
        )[0, 1]
        rmse = np.sqrt(np.mean((phi_t - phi_recon) ** 2))
        correlations.append(corr)
        rmses.append(rmse)

        signal = phi_t.max() - phi_t.min()
        print(f"  h = {h*1e6:4.0f} µm:  corr = {corr:.4f},  "
              f"RMSE/Δ = {rmse/(signal+1e-30)*100:5.1f}%,  "
              f"Φ_peak = {signal*1e6:.2f} µV")

    # Resolution plot
    fig3, (ax3a, ax3b) = plt.subplots(1, 2, figsize=(12, 5))

    ax3a.semilogx(heights * 1e6, correlations, 'ko-', ms=8)
    ax3a.set_xlabel('Ion height $h$ [µm]')
    ax3a.set_ylabel('Pearson correlation (true vs. recon. $\\Phi$ at ion)')
    ax3a.set_title('Reconstruction quality vs. ion height')
    ax3a.axhline(y=0.9, color='gray', ls='--', alpha=0.5, label='0.9 threshold')
    ax3a.legend()
    ax3a.grid(True, alpha=0.3)
    ax3a.set_ylim(-0.1, 1.1)

    ax3b.loglog(heights * 1e6, rmses, 'rs-', ms=8, label='RMSE [V]')
    ax3b.set_xlabel('Ion height $h$ [µm]')
    ax3b.set_ylabel('RMSE [V]')
    ax3b.set_title('Reconstruction error vs. ion height')
    # Fit power law
    coeffs = np.polyfit(np.log(heights), np.log(rmses), 1)
    ax3b.loglog(
        heights * 1e6,
        np.exp(coeffs[1]) * heights ** coeffs[0],
        'r--', alpha=0.5,
        label=f'RMSE ∝ h^{coeffs[0]:.1f}',
    )
    ax3b.legend()
    ax3b.grid(True, alpha=0.3)

    fig3.suptitle(
        "Spatial Resolution Limits: Ion Height Dependence\n"
        "Forward model $e^{-kh}$ low-pass filtering dominates",
        fontsize=13, fontweight='bold',
    )
    plt.tight_layout()
    plt.savefig("resolution_analysis.pdf", dpi=150, bbox_inches='tight')
    print("Figure saved to resolution_analysis.pdf")

    print("\n✓ All simulations complete.")
