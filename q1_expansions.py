"""
q1_expansions.py — Q1 Journal Expansion Computations
=====================================================
Implements key theoretical additions for Q1-level manuscript:

  1. Quantum Fisher Information for surface height estimation
     (coherent, squeezed, and NOON states)
  2. Singular Value Decomposition of the discretised forward operator
     (effective rank, information content, resolution matrix)
  3. Micromotion ITF — first-derivative complementary channel
  4. Two-height condition number vs. spatial frequency

Generates figures: figQ1_qfi.pdf, figQ1_svd.pdf, figQ1_micromotion.pdf

Author: QESPM project
"""

import numpy as np
from numpy.fft import fft2, ifft2, fftfreq
from dataclasses import dataclass
from typing import Optional, Tuple, Dict
import matplotlib
import matplotlib.pyplot as plt

# =====================================================================
# Constants and parameters
# =====================================================================
E_CHARGE = 1.602176634e-19
M_CA40   = 40.0 * 1.660539e-27
EPS0     = 8.8541878128e-12
HBAR     = 1.054571817e-34

@dataclass
class QESPMParams:
    ion_mass: float = M_CA40
    freq_sec: float = 1.0e6       # Hz
    E_bg: float = 1.0e4           # V/m
    ion_height: float = 40e-6     # m
    eta: float = 0.1              # Lamb-Dicke parameter
    omega_carrier: float = 2 * np.pi * 1e6  # carrier Rabi freq, rad/s

    @property
    def omega_sec(self):
        return 2 * np.pi * self.freq_sec

    @property
    def prefactor(self):
        return E_CHARGE / (2 * self.ion_mass * self.omega_sec)


# =====================================================================
# 1. Quantum Fisher Information Analysis
# =====================================================================

def qfi_coherent_state(
    A: float, k: float, h: float, tau: float, params: QESPMParams
) -> float:
    """
    Quantum Fisher Information for estimating surface amplitude A from
    a coherent state |α⟩ evolving for time τ.

    Measurement model: Δω = C·A·k²·e^{-kh}, where C = e·E_bg/(2mω).
    For a coherent state with amplitude α, the QFI for frequency estimation is:
        F_ω = 4|α|² τ²   (in units of rad⁻²·s² → (rad/s)⁻²)

    Then by error propagation: F_A = F_ω · (∂Δω/∂A)².

    Returns
    -------
    F_A : float
        QFI for surface amplitude A [(m)⁻²].
    delta_A_min : float
        Minimum detectable amplitude (Cramér--Rao bound) [m].
    """
    alpha = 10.0  # coherent state amplitude (typical for Doppler-cooled ion)
    C = params.prefactor * params.E_bg
    domega_dA = C * k**2 * np.exp(-k * h)

    # QFI for frequency: F_ω = 4|α|² τ² (coherent state)
    F_omega = 4.0 * alpha**2 * tau**2

    # Propagate to amplitude
    F_A = F_omega * domega_dA**2
    delta_A_min = 1.0 / np.sqrt(max(F_A, 1e-300))

    return F_A, delta_A_min


def qfi_squeezed_state(
    A: float, k: float, h: float, tau: float, params: QESPMParams, r: float = 1.5
) -> tuple:
    """
    QFI for squeezed vacuum with squeezing parameter r.
    F_ω^{(sq)} = e^{2r} · F_ω^{(coh)}  (for phase estimation in the squeezed quadrature).
    """
    F_A_coh, delta_A_coh = qfi_coherent_state(A, k, h, tau, params)
    F_A_sq = np.exp(2 * r) * F_A_coh
    delta_A_sq = 1.0 / np.sqrt(max(F_A_sq, 1e-300))
    return F_A_sq, delta_A_sq


def qfi_noon_state(
    A: float, k: float, h: float, tau: float, params: QESPMParams, N: int = 20
) -> tuple:
    """
    QFI for a NOON state |N,0⟩ + |0,N⟩ with N quanta.
    Heisenberg limit: F_ω^{(NOON)} = N² · F_ω^{(single quantum)}.

    For fair comparison with coherent state |α⟩ having ⟨n⟩ = |α|²:
    use N = ⟨n⟩_coh = |α|² so both states have the same average quanta.
    With |α|² = 100 and N = 100: enhancement = N²/(4|α|²) = N/4 = 25×.
    """
    C = params.prefactor * params.E_bg
    domega_dA = C * k**2 * np.exp(-k * h)
    F_omega_1 = 4.0 * 1.0 * tau**2  # single quantum Fisher information
    F_omega_N = N**2 * F_omega_1    # Heisenberg scaling
    F_A_noon = F_omega_N * domega_dA**2
    delta_A_noon = 1.0 / np.sqrt(max(F_A_noon, 1e-300))
    return F_A_noon, delta_A_noon


# =====================================================================
# 2. Singular Value Decomposition of Forward Operator
# =====================================================================

def build_forward_matrix_1d(
    nx: int, dx: float, h: float, params: QESPMParams
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build the 1D discretised forward operator H: z_s → Δω_x.

    In 1D, H is a circulant matrix (Toeplitz) whose action in Fourier
    space is H(k) = -C · k² · e^{-k|h|}.

    Returns
    -------
    H_matrix : (nx, nx) ndarray
        Forward operator matrix.
    k_vals : (nx,) ndarray
        Wavenumbers [rad/m].
    """
    C = params.prefactor * params.E_bg
    k_vals = 2 * np.pi * fftfreq(nx, dx)

    # Build circulant matrix: H[i,j] = h(x_i - x_j)
    H_k = -C * k_vals**2 * np.exp(-np.abs(k_vals) * h)
    H_k[0] = 0.0  # DC mode

    # First column of the circulant matrix (inverse FFT of H_k)
    h_col = np.real(ifft2(H_k.reshape(-1, 1))).ravel()

    # Build Toeplitz/circulant matrix
    from scipy.linalg import circulant
    H_matrix = circulant(h_col)

    return H_matrix, k_vals


def svd_analysis(
    H_matrix: np.ndarray, k_vals: np.ndarray, params: QESPMParams
) -> Dict:
    """
    Perform SVD analysis of the forward operator.

    Returns dictionary with singular values, effective rank, resolution
    matrix, and information content.
    """
    U, S, Vh = np.linalg.svd(H_matrix, full_matrices=False)

    # Effective rank: number of singular values above ε·S_max
    eps_rank = 1e-3
    effective_rank = np.sum(S > eps_rank * S[0])

    # Information content: I = -½ Σ log(1 - R_ii)
    # where R = V V^† (for no regularisation) or with Tikhonov
    lam = 1e-4 * S[0]
    # Tikhonov-filtered singular values
    S_filt = S**2 / (S**2 + lam**2)  # filter factors
    R_diag = np.sum((Vh.T)**2 * S_filt, axis=1)  # diagonal of resolution matrix
    R_diag = np.clip(R_diag, 0, 0.9999)
    info_content = -0.5 * np.sum(np.log(1.0 - R_diag))

    return {
        'S': S,
        'effective_rank': effective_rank,
        'info_content': info_content,
        'S_filt': S_filt,
        'U': U,
        'Vh': Vh,
    }


# =====================================================================
# 3. Micromotion ITF
# =====================================================================

def micromotion_itf_1d(
    k: np.ndarray, h: float, params: QESPMParams
) -> np.ndarray:
    """
    ITF for micromotion amplitude: β_x ∝ |ik_x e^{-kh} φ̃_s|.

    The micromotion amplitude β_x(x) ∝ |∂Φ/∂x|.
    In Fourier space: ∂Φ/∂x → ik_x · Φ̃.
    So the effective transfer function for β_x is:
        H_β(k) ∝ |k_x| e^{-k h}.

    Returns |H_β(k)| (arbitrary normalisation).
    """
    return np.abs(k) * np.exp(-np.abs(k) * h)


def compare_itf_channels(
    k: np.ndarray, h: float, params: QESPMParams
) -> Dict:
    """Compare Δω and β ITF magnitudes."""
    C = params.prefactor * params.E_bg
    H_dw = np.abs(-C * k**2 * np.exp(-np.abs(k) * h))
    H_beta = micromotion_itf_1d(k, h, params)

    # Normalise to peak = 1 for comparison
    H_dw_norm = H_dw / np.max(H_dw[1:]) if np.max(H_dw[1:]) > 0 else H_dw
    H_beta_norm = H_beta / np.max(H_beta[1:]) if np.max(H_beta[1:]) > 0 else H_beta

    return {
        'H_dw': H_dw_norm,
        'H_beta': H_beta_norm,
        'k': k,
    }


# =====================================================================
# 4. Two-height condition number
# =====================================================================

def two_height_condition_number(
    k: np.ndarray, h1: float, h2: float
) -> np.ndarray:
    """
    Condition number of the two-height measurement matrix for
    topography--charge separation.

    The measurement matrix M(k) couples (z_s, φ_charge) to
    (Δω(h1), Δω(h2)).  The condition number κ(M) governs
    how well the two contributions can be separated.
    """
    # M = [[e^{-k h1}, e^{-k h1}],
    #      [e^{-k h2}, e^{-k h2}]]
    # Actually both columns are proportional — this matrix is rank-1!
    # The trick is that the prefactors differ: topographic contribution
    # ∝ E_bg, charge contribution ∝ 1.  So:
    # M ∝ [[E_bg e^{-k h1}, e^{-k h1}],
    #       [E_bg e^{-k h2}, e^{-k h2}]]
    # This is STILL rank-1 if E_bg is the same at both heights!
    #
    # The real separability comes from the fact that φ_charge(k) is the
    # surface potential, which is propagated by e^{-kz} just like z_s.
    # Separation requires a DIFFERENT h-dependence, e.g.:
    # - Topography → Δω via k² e^{-kh} · E_bg
    # - Charge → Δω via k² e^{-kh} · (1/ε₀) · σ · d_eff
    # If E_bg ≠ (d_eff/ε₀), the columns differ.
    # But in our model, both scale as e^{-kh} — separation is IMPOSSIBLE
    # from heights alone if the charge is unknown!
    #
    # The resolution: charge produces a SURFACE potential φ_charge = σ·d_eff/ε₀
    # which propagates identically to topography.  At a single height,
    # they are indistinguishable.  At two heights, we get:
    # Δω(h1) = C·k²·e^{-k h1}·(E_bg·z_s + φ_charge)
    # Δω(h2) = C·k²·e^{-k h2}·(E_bg·z_s + φ_charge)
    # Ratio: Δω(h1)/Δω(h2) = e^{-k(h1-h2)} — SAME for both components!
    #
    # This proves that TWO-HEIGHT SCANNING ALONE CANNOT SEPARATE
    # TOPOGRAPHY FROM CHARGE!  This is a key theoretical result.
    #
    # Additional information is needed: e.g. known charge distribution,
    # or in-situ charge elimination (UV, ion bombardment), or
    # independent measurement of φ_charge.

    # The condition number is formally infinite (rank-1 matrix).
    # Return a diagnostic.
    k = np.abs(k)
    ratio = np.exp(-k * np.abs(h1 - h2))
    return ratio, np.ones_like(k) * np.inf  # formally infinite condition number


# =====================================================================
# Main: generate Q1 expansion figures
# =====================================================================
if __name__ == "__main__":
    matplotlib.use("Agg")
    plt.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 150, "savefig.bbox": "tight",
        "font.family": "serif", "font.size": 11,
        "mathtext.fontset": "stix",
    })

    params = QESPMParams()

    # ================================================================
    # Figure Q1: Quantum Fisher Information
    # ================================================================
    k_vals = np.logspace(3, 6.5, 300)
    h_test = 40e-6
    A_test = 100e-9
    tau_test = 10e-3  # 10 ms

    delta_A_coh = np.zeros(len(k_vals))
    delta_A_sq = np.zeros(len(k_vals))
    delta_A_noon = np.zeros(len(k_vals))

    for i, k_i in enumerate(k_vals):
        _, dA = qfi_coherent_state(A_test, k_i, h_test, tau_test, params)
        delta_A_coh[i] = dA
        _, dA = qfi_squeezed_state(A_test, k_i, h_test, tau_test, params, r=1.5)
        delta_A_sq[i] = dA
        _, dA = qfi_noon_state(A_test, k_i, h_test, tau_test, params, N=100)
        delta_A_noon[i] = dA

    fig1, (ax1a, ax1b) = plt.subplots(1, 2, figsize=(14, 6))

    # Panel (a): minimum detectable amplitude vs k
    ax1a.loglog(k_vals, delta_A_coh * 1e9, "b-", lw=2, label="$|\\alpha=10\\rangle$ (classical)")
    ax1a.loglog(k_vals, delta_A_sq * 1e9, "r--", lw=2, label="Squeezed ($r=1.5$, $\\sim 13$ dB)")
    ax1a.loglog(k_vals, delta_A_noon * 1e9, "g:", lw=2, label="NOON state ($N=100$, Heisenberg)")
    ax1a.axhline(A_test * 1e9, color="gray", ls=":", alpha=0.5)
    ax1a.text(2e3, A_test * 1e9 * 1.3, "$A = 100$ nm", fontsize=9, color="gray")
    ax1a.set_xlabel("Wavenumber $k$ [rad/m]", fontsize=12)
    ax1a.set_ylabel("$\\delta A_{\\min}$ (CRB) [nm]", fontsize=12)
    ax1a.set_title("(a) Quantum-limited surface amplitude sensitivity\n"
                   f"$h = {h_test*1e6:.0f}$ µm, $\\tau = {tau_test*1e3:.0f}$ ms",
                   fontsize=12, fontweight="bold")
    ax1a.legend(fontsize=9, loc="upper left")
    ax1a.grid(True, alpha=0.3, lw=0.5, which="both")

    # Panel (b): squeezing enhancement factor
    enhancement_sq = delta_A_coh / delta_A_sq
    enhancement_noon = delta_A_coh / delta_A_noon
    ax1b.semilogx(k_vals, enhancement_sq, "r-", lw=2, label=f"Squeezing ($r=1.5$)")
    ax1b.semilogx(k_vals, enhancement_noon, "g--", lw=2, label=f"NOON ($N=100$)")
    ax1b.axhline(np.exp(1.5), color="r", ls=":", alpha=0.3)
    ax1b.axhline(25.0, color="g", ls=":", alpha=0.3)
    ax1b.set_xlabel("Wavenumber $k$ [rad/m]", fontsize=12)
    ax1b.set_ylabel("Sensitivity enhancement $\\delta A_{\\rm coh}/\\delta A_{\\rm q}$",
                    fontsize=12)
    ax1b.set_title("(b) Quantum enhancement factor", fontsize=12, fontweight="bold")
    ax1b.grid(True, alpha=0.3, lw=0.5)

    fig1.suptitle("Quantum-Enhanced Surface Topography Estimation\n"
                  "Quantum Fisher Information analysis for ${}^{40}$Ca$^+$ ion probe",
                  fontsize=14, fontweight="bold")
    plt.tight_layout()
    # Legend after tight_layout to prevent repositioning
    ax1b.legend(fontsize=9, loc="upper right")
    fig1.savefig("figQ1_qfi.pdf")
    print("Saved figQ1_qfi.pdf")
    print(f"  At k_opt = 2/h = {2/h_test:.1e} rad/m:")
    k_opt_idx = np.argmin(np.abs(k_vals - 2/h_test))
    print(f"    Coherent:  δA_min = {delta_A_coh[k_opt_idx]*1e9:.2e} nm")
    print(f"    Squeezed:  δA_min = {delta_A_sq[k_opt_idx]*1e9:.2e} nm  "
          f"(×{enhancement_sq[k_opt_idx]:.1f} better)")
    print(f"    NOON:      δA_min = {delta_A_noon[k_opt_idx]*1e9:.2e} nm  "
          f"(×{enhancement_noon[k_opt_idx]:.1f} better)")

    # ================================================================
    # Figure Q2: SVD Analysis
    # ================================================================
    nx = 128
    dx = 0.5e-6
    heights_svd = [5e-6, 10e-6, 20e-6, 40e-6, 80e-6]

    fig2, (ax2a, ax2b) = plt.subplots(1, 2, figsize=(14, 6))
    colours = plt.cm.viridis(np.linspace(0.15, 0.9, len(heights_svd)))

    for h_i, c in zip(heights_svd, colours):
        H_mat, k_vec = build_forward_matrix_1d(nx, dx, h_i, params)
        svd_result = svd_analysis(H_mat, k_vec, params)
        S = svd_result['S']
        rank_eff = svd_result['effective_rank']
        info = svd_result['info_content']

        # Plot singular value spectrum
        ax2a.semilogy(np.arange(1, len(S) + 1), S / S[0], color=c, lw=1.5,
                       label=f"$h = {h_i*1e6:.0f}$ µm (rank={rank_eff})")

    ax2a.set_xlabel("Singular value index $i$", fontsize=12)
    ax2a.set_ylabel("Normalised singular value $\\sigma_i/\\sigma_1$", fontsize=12)
    ax2a.set_title("(a) Singular value spectrum of forward operator $\\mathbf{H}$",
                   fontsize=12, fontweight="bold")
    ax2a.legend(fontsize=9, loc="upper right")
    ax2a.grid(True, alpha=0.3, lw=0.5)
    ax2a.set_ylim(1e-18, 2)

    # Add info content vs h
    h_fine = np.logspace(-6, -3, 100)
    info_fine = []
    rank_fine = []
    for h_j in h_fine:
        H_mat_j, k_vec_j = build_forward_matrix_1d(nx, dx, h_j, params)
        svd_j = svd_analysis(H_mat_j, k_vec_j, params)
        info_fine.append(svd_j['info_content'])
        rank_fine.append(svd_j['effective_rank'])

    ax2b_twin = ax2b.twinx()
    ax2b.loglog(h_fine * 1e6, info_fine, "b-", lw=2, label="Info content $I(h)$")
    ax2b_twin.loglog(h_fine * 1e6, rank_fine, "r--", lw=2, label="Effective rank")
    ax2b.set_xlabel("Ion height $h$ [µm]", fontsize=12)
    ax2b.set_ylabel("Shannon information $I(h)$ [nats]", fontsize=12, color="b")
    ax2b_twin.set_ylabel("Effective rank", fontsize=12, color="r")
    ax2b.set_title("(b) Information content and effective rank vs.~$h$\n"
                   f"$N = {nx}$, $\\Delta x = {dx*1e6:.1f}$ µm, FOV = {nx*dx*1e6:.0f} µm",
                   fontsize=11, fontweight="bold")
    ax2b.grid(True, alpha=0.3, lw=0.5, which="both")
    ax2b.tick_params(axis="y", colors="b")
    ax2b_twin.tick_params(axis="y", colors="r")

    # Combined legend for twin axes
    lines1, labels1 = ax2b.get_legend_handles_labels()
    lines2, labels2 = ax2b_twin.get_legend_handles_labels()
    ax2b.legend(lines1 + lines2, labels1 + labels2, fontsize=9, loc="lower left")

    fig2.suptitle("Inverse Problem Structure: SVD Analysis of the Forward Operator",
                  fontsize=14, fontweight="bold")
    plt.tight_layout()
    fig2.savefig("figQ2_svd.pdf")
    print("Saved figQ2_svd.pdf")
    for h_i in heights_svd:
        H_mat, k_vec = build_forward_matrix_1d(nx, dx, h_i, params)
        svd_r = svd_analysis(H_mat, k_vec, params)
        print(f"  h = {h_i*1e6:5.0f} µm:  rank = {svd_r['effective_rank']:3d},  "
              f"I = {svd_r['info_content']:.2f} nats")

    # ================================================================
    # Figure Q3: Micromotion ITF Comparison
    # ================================================================
    k_plot = np.logspace(3, 8, 500)
    h_mm = 40e-6
    channels = compare_itf_channels(k_plot, h_mm, params)

    fig3, (ax3a, ax3b) = plt.subplots(1, 2, figsize=(14, 6))

    # Panel (a): ITF comparison
    ax3a.loglog(k_plot, channels['H_dw'], "b-", lw=2,
                label="$\\Delta\\omega_x$ (curvature, $\\propto k^2 e^{-kh}$)")
    ax3a.loglog(k_plot, channels['H_beta'], "r-", lw=2,
                label="$\\beta_x$ (gradient, $\\propto |k| e^{-kh}$)")
    ax3a.set_xlabel("Wavenumber $k$ [rad/m]", fontsize=12)
    ax3a.set_ylabel("Normalised ITF magnitude", fontsize=12)
    ax3a.set_title("(a) Complementary ITF channels\n"
                   f"$h = {h_mm*1e6:.0f}$ µm",
                   fontsize=12, fontweight="bold")
    ax3a.legend(fontsize=9)
    ax3a.grid(True, alpha=0.3, lw=0.5, which="both")
    ax3a.set_ylim(1e-4, 2)

    # Panel (b): derivative order comparison
    # H_dw ∝ k² e^{-kh}: peaks at k=2/h
    # H_beta ∝ k e^{-kh}: peaks at k=1/h
    k_opt_dw = 2.0 / h_mm
    k_opt_beta = 1.0 / h_mm
    lam_opt_dw = 2 * np.pi / k_opt_dw
    lam_opt_beta = 2 * np.pi / k_opt_beta

    # Illustrate the complementary coverage
    k_band = np.logspace(3, 7, 200)
    ax3b.fill_between(k_band, 0, 1, where=(k_band > 1e3), color="blue", alpha=0.1,
                       label="$\\Delta\\omega_x$ dominant")
    ax3b.fill_between(k_band, 0, 1, where=(k_band < 1e5), color="red", alpha=0.1,
                       label="$\\beta_x$ dominant")
    ax3b.axvline(k_opt_dw, color="blue", ls="--", lw=1.2,
                  label=f"$\\Delta\\omega_x$ peak ($\\lambda={lam_opt_dw*1e6:.0f}$ µm)")
    ax3b.axvline(k_opt_beta, color="red", ls="--", lw=1.2,
                  label=f"$\\beta_x$ peak ($\\lambda={lam_opt_beta*1e6:.0f}$ µm)")
    ax3b.set_xscale("log")
    ax3b.set_xlabel("Wavenumber $k$ [rad/m]", fontsize=12)
    ax3b.set_ylabel("Relative sensitivity", fontsize=12)
    ax3b.set_title("(b) Spatial frequency coverage\n"
                   "gradient + curvature = full $k$-band",
                   fontsize=12, fontweight="bold")
    ax3b.legend(fontsize=8, loc="upper right")
    ax3b.grid(True, alpha=0.3, lw=0.5)
    ax3b.set_ylim(0, 1.1)

    fig3.suptitle("Multi-Signal Fusion: Micromotion as Complementary Observable",
                  fontsize=14, fontweight="bold")
    plt.tight_layout()
    fig3.savefig("figQ3_micromotion.pdf")
    print("Saved figQ3_micromotion.pdf")

    # ================================================================
    # Figure Q4: Two-height analysis — proof of degeneracy
    # ================================================================
    k_proof = np.logspace(3, 7, 300)
    h1_proof = 20e-6
    h2_proof = 80e-6
    ratio, kappa = two_height_condition_number(k_proof, h1_proof, h2_proof)

    fig4, ax4 = plt.subplots(1, 1, figsize=(8, 5))
    ax4.semilogx(k_proof, ratio, "k-", lw=2)
    ax4.set_xlabel("Wavenumber $k$ [rad/m]", fontsize=12)
    ax4.set_ylabel("Ratio $\\Delta\\omega(h_1)/\\Delta\\omega(h_2) = e^{-k(h_1-h_2)}$",
                   fontsize=12)
    ax4.set_title(
        "Two-Height Degeneracy: topography and charge\n"
        "share identical $h$-scaling",
        fontsize=12, fontweight="bold",
    )
    ax4.grid(True, alpha=0.3, lw=0.5)
    ax4.text(0.5, 0.5,
             "Two-height scanning alone\n"
             "CANNOT separate topography\n"
             "from surface charge.\n\n"
             "Additional information required:\n"
             "• In-situ charge elimination (UV, ion bombardment)\n"
             "• Known charge distribution (reference scan)\n"
             "• Independent material-contrast channel ($\\dot{\\bar{n}}$)",
             transform=ax4.transAxes, fontsize=10, va="center", ha="center",
             bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))

    fig4.tight_layout()
    fig4.savefig("figQ4_two_height_proof.pdf")
    print("Saved figQ4_two_height_proof.pdf")
    print()
    print("Key theoretical result: two-height scanning alone cannot separate")
    print("topography from charge. Both contributions share the e^{-kh} scaling.")
    print("The forward matrix is rank-1 → condition number is formally infinite.")
    print("Resolution strategies: in-situ charge elimination, known reference,")
    print("or independent material-contrast channel (heating rate).")
    print()

    # ================================================================
    # Numerical summary
    # ================================================================
    print("=" * 65)
    print("  Q1 Expansion Results Summary")
    print("=" * 65)
    print(f"  Quantum enhancement at k_opt, h = {h_test*1e6:.0f} µm:")
    print(f"    Coherent state:  δA_min = {delta_A_coh[k_opt_idx]*1e9:.2e} nm (baseline)")
    print(f"    Squeezed (13dB): δA_min = {delta_A_sq[k_opt_idx]*1e9:.2e} nm "
          f"(×{delta_A_coh[k_opt_idx]/delta_A_sq[k_opt_idx]:.1f})")
    print(f"    NOON (N=5):      δA_min = {delta_A_noon[k_opt_idx]*1e9:.2e} nm "
          f"(×{delta_A_coh[k_opt_idx]/delta_A_noon[k_opt_idx]:.1f})")
    print()
    print(f"  SVD effective rank vs. h (N = {nx}):")
    for h_i in heights_svd:
        H_mat, k_vec = build_forward_matrix_1d(nx, dx, h_i, params)
        svd_r = svd_analysis(H_mat, k_vec, params)
        print(f"    h = {h_i*1e6:5.0f} µm:  rank = {svd_r['effective_rank']:3d},  "
              f"I = {svd_r['info_content']:.2f} nats")
    print()

    print("✓ All Q1 expansion figures generated:")
    print("  figQ1_qfi.pdf            — Quantum Fisher Information analysis")
    print("  figQ2_svd.pdf            — SVD of forward operator")
    print("  figQ3_micromotion.pdf     — Micromotion ITF comparison")
    print("  figQ4_two_height_proof.pdf — Proof: two-height cannot separate charge")
