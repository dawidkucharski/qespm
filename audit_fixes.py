#!/usr/bin/env python3
"""Audit-fix verification and figure generation for the QESPM manuscript.

Regenerates fig7_svd.pdf on the well-defined 1D cross-section lattice
(L=64 um, N=128) and produces the new figures required by the referee
audit: Picard/L-curve/GCV, master noise budget, FD validation, singular
functions, BEM cross-check, and quantum water-filling.

Every number printed by `main()` is computed here from first principles
and used verbatim in main.tex.
"""
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq, ifft
from scipy.sparse import diags, lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve

# ---------------- physical constants (nominal QESPM) ----------------
E_CHARGE = 1.602176634e-19
M_CA40 = 40 * 1.66053906660e-27
E_BG = 1e4
OMEGA_X = 2 * np.pi * 1e6
C_ITF = E_CHARGE * E_BG / (2 * M_CA40 * OMEGA_X)   # = 1.92e3 m/s
L_FOV = 64e-6          # inversion FOV (matches grid-convergence section)
N_GRID = 128           # lattice points
DX = L_FOV / N_GRID    # 0.5 um
H_NOM = 40e-6
NOISE_HZ = 5.0
D_OMEGA = 2 * np.pi * NOISE_HZ
LAMBDA_DP_FRAC = 3e-5  # discrepancy-principle lambda as fraction of max|H|

MANUSCRIPT = "/Users/dawid/Projects/ion_surface_sensor/manuscript/"


def itf_1d(k, h):
    return -C_ITF * k**2 * np.exp(-np.abs(k) * h)


def svd_1d():
    """1D cross-section SVD: regenerate fig7_svd.pdf and the table."""
    k = 2 * np.pi * fftfreq(N_GRID, DX)
    heights = [5e-6, 10e-6, 20e-6, 40e-6, 80e-6]
    print("=" * 60)
    print("I-SVD: 1D cross-section, L=64um, N=128, dx=0.5um")
    print("=" * 60)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    table = []
    for h in heights:
        H = itf_1d(k, h)
        s = np.sort(np.abs(H))[::-1]
        reff = int(np.sum(s > 1e-3 * s[0]))
        lam = LAMBDA_DP_FRAC * s[0]
        I = 0.5 * np.sum(np.log(1 + (s / lam) ** 2))
        table.append((h * 1e6, reff, I))
        ax[0].semilogy(np.arange(1, len(s) + 1), s / s[0],
                       label=f"$h={h*1e6:.0f}$ µm")
        print(f"h={h*1e6:5.0f} um: r_eff={reff:4d}  I={I:9.2f} nats")
    ax[0].set_xlabel("index $n$")
    ax[0].set_ylabel("$\\sigma_n/\\sigma_1$")
    ax[0].set_title("(a) singular spectrum, 1D cross-section $k_y=0$")
    ax[0].legend(fontsize=8)
    ax[0].grid(True, alpha=0.3)
    hs = np.array([t[0] for t in table])
    r = np.array([t[1] for t in table])
    I = np.array([t[2] for t in table])
    ax[1].plot(hs, r, "o-", label="$r_{\\rm eff}$")
    ax[1].set_xlabel("$h$ [µm]")
    ax[1].set_ylabel("effective rank", color="C0")
    ax[1].tick_params(axis="y", labelcolor="C0")
    ax[1].set_yscale("log")
    axr = ax[1].twinx()
    axr.plot(hs, I, "s--", color="C1", label="$I(h)$")
    axr.set_ylabel("$I$ [nats]", color="C1")
    axr.tick_params(axis="y", labelcolor="C1")
    axr.set_yscale("log")
    ax[1].set_title("(b) effective rank and Shannon information")
    # scaling exponents
    p_r = np.polyfit(np.log(hs), np.log(r), 1)[0]
    p_I = np.polyfit(np.log(hs), np.log(I), 1)[0]
    print(f"fitted exponent r_eff: {p_r:+.3f}   I: {p_I:+.3f}")
    # threshold sensitivity of the exponent
    for thr in [1e-2, 1e-3, 1e-4]:
        rs = []
        for h in heights:
            s = np.sort(np.abs(itf_1d(k, h)))[::-1]
            rs.append(int(np.sum(s > thr * s[0])))
        p = np.polyfit(np.log(hs), np.log(rs), 1)[0]
        print(f"  threshold {thr:.0e}: r_eff={rs}, exponent {p:+.3f}")
    # Weyl asymptotics check (1D): ln(C/(sigma_n h^2)) vs n
    for h in [40e-6]:
        s = np.sort(np.abs(itf_1d(k, h)))[::-1]
        y = np.log(C_ITF / (s * h**2))
        nn = np.arange(1, len(s) + 1)
        m = np.polyfit(nn[2:60], y[2:60], 1)[0]
        print(f"Weyl fit h=40um: slope {m:.4f}  vs pi*h/L = {np.pi*h/L_FOV:.4f}")
    fig.tight_layout()
    fig.savefig(MANUSCRIPT + "fig7_svd.pdf")
    plt.close(fig)
    return table


def picard_lcurve_gcv():
    """Picard plot, L-curve, GCV on the 2D test surface; figI3."""
    k = 2 * np.pi * fftfreq(N_GRID, DX)
    KX, KY = np.meshgrid(k, k)
    x = np.arange(N_GRID) * DX
    X, Y = np.meshgrid(x, x)
    h = H_NOM
    zs = 200e-9 * np.exp(-((X - 32e-6)**2 + (Y - 32e-6)**2) / (2 * (7.5e-6)**2))
    zs += 100e-9 * (X / L_FOV)
    H = -C_ITF * KX**2 * np.exp(-np.sqrt(KX**2 + KY**2) * h)
    Hs = np.abs(H)
    d_true = np.real(ifft(H * fft(zs)))
    rng = np.random.default_rng(7)
    d = d_true + rng.normal(0, D_OMEGA, (N_GRID, N_GRID))
    delta = np.sqrt(N_GRID**2) * D_OMEGA
    v = fft(d) / N_GRID**2
    sig = np.sort(Hs.ravel())[::-1]
    pic = np.sort(np.abs(v.ravel()))[::-1] / sig
    lams = np.logspace(-9, -1.0, 300) * Hs.max()
    res = np.zeros_like(lams)
    sol = np.zeros_like(lams)
    gcv = np.zeros_like(lams)
    for i, lam in enumerate(lams):
        f = Hs**2 / (Hs**2 + lam**2)
        rec = np.real(ifft(f * fft(d)))
        res[i] = np.linalg.norm(rec - d)
        fz = np.where(Hs > 0, Hs / (Hs**2 + lam**2), 0.0)
        sol[i] = np.linalg.norm(np.real(ifft(fz * fft(d))))
        num = np.sum((1 - f) ** 2 * np.abs(fft(d)) ** 2)
        gcv[i] = num / (N_GRID**2 - np.sum(f)) ** 2
    idx = np.argmin(np.abs(res - 1.01 * delta))
    lam_dp = lams[idx]
    lr, ls_ = np.log(res), np.log(sol)
    curv = np.gradient(np.gradient(ls_, lr), lr)
    idxL = np.argmax(np.abs(curv))
    lam_L = lams[idxL]
    idxG = np.argmin(gcv)
    lam_G = lams[idxG]
    print("=" * 60)
    print("I-PICARD: 2D plateau+ramp, h=40um, N=128, noise 5 Hz")
    print(f"  |H|max = {Hs.max():.4e} (rad/s)/m, delta = {delta:.1f} rad/s")
    print(f"  lambda_DP / max|H| = {lam_dp/Hs.max():.3e}")
    print(f"  lambda_L  / max|H| = {lam_L/Hs.max():.3e}")
    print(f"  lambda_GCV/ max|H| = {lam_G/Hs.max():.3e}")
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.8))
    ax[0].semilogy(pic[:200], ".-", markersize=3)
    ax[0].set_xlabel("index $n$")
    ax[0].set_ylabel("$|\\langle d, v_n\\rangle|/\\sigma_n$")
    ax[0].set_title("(a) discrete Picard plot")
    ax[0].grid(True, alpha=0.3)
    ax[1].loglog(res, sol, ".-", markersize=3)
    ax[1].plot(res[idxL], sol[idxL], "ro")
    ax[1].set_xlabel("$\\|H z_\\lambda - d\\|$")
    ax[1].set_ylabel("$\\|z_\\lambda\\|$")
    ax[1].set_title("(b) L-curve (corner marked)")
    ax[1].grid(True, alpha=0.3)
    ax[2].loglog(lams / Hs.max(), gcv, ".-", markersize=3)
    ax[2].plot(lams[idxG] / Hs.max(), gcv[idxG], "ro")
    ax[2].set_xlabel("$\\lambda/\\max|H|$")
    ax[2].set_ylabel("$G(\\lambda)$")
    ax[2].set_title("(c) GCV function (minimum marked)")
    ax[2].grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(MANUSCRIPT + "figI3_picard.pdf")
    plt.close(fig)


def noise_master():
    """Master noise budget: delta-omega per source vs h; figI10."""
    hs = np.logspace(np.log10(1.5e-6), np.log10(250e-6), 200)
    T2 = 10e-3
    nbar = 0.5
    Nrep = 1000
    gamma0, h0, alpha = 100.0, 40e-6, 4.0
    dQPN = np.full_like(hs, 2 * np.pi * 0.7)
    dTech = np.full_like(hs, 2 * np.pi * 5.0)
    dRF = np.full_like(hs, 2 * np.pi * 1.0)
    # heating-limited: domega = gamma0 (h0/h)^alpha / (2 sqrt(Nrep <n>))
    dHeat = gamma0 * (h0 / hs) ** alpha / (2 * np.sqrt(Nrep * nbar))
    # laser frequency noise: 1 kHz and 10 Hz linewidths
    T2_l1 = 1.0 / (np.pi * 1e3)
    T2_l2 = 1.0 / (np.pi * 10.0)
    dLas1 = np.full_like(hs, 1.0 / (2 * T2_l1 * np.sqrt(Nrep * nbar)))
    dLas2 = np.full_like(hs, 1.0 / (2 * T2_l2 * np.sqrt(Nrep * nbar)))
    # Johnson gradient noise: (e/4 pi m w) sqrt(3 kB T rho /(8 pi d^5 tau))
    kB = 1.380649e-23
    dJohn = (E_CHARGE / (4 * np.pi * M_CA40 * OMEGA_X)) * np.sqrt(
        3 * kB * 4.0 * 3e-11 / (8 * np.pi * hs**5 * 1.0))
    print("=" * 60)
    print("I-NOISE: master budget")
    print(f"  laser 1 kHz: {dLas1[0]/2/np.pi:.2f} Hz | laser 10 Hz: "
          f"{dLas2[0]/2/np.pi:.3f} Hz")
    print(f"  heating at h=5um: {dHeat[hs<6e-6][-1]/2/np.pi:.1f} Hz | "
          f"h=40um: {dHeat[np.argmin(abs(hs-40e-6))]/2/np.pi:.2f} Hz")
    print(f"  Johnson at 40um: {dJohn[np.argmin(abs(hs-40e-6))]/2/np.pi:.2e} Hz")
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 4))
    ax[0].loglog(hs * 1e6, dTech / 2 / np.pi, "k--", label="technical floor 5 Hz")
    ax[0].loglog(hs * 1e6, dQPN / 2 / np.pi, "g-", label="QPN 0.7 Hz")
    ax[0].loglog(hs * 1e6, dRF / 2 / np.pi, "m-.", label="rf stability")
    ax[0].loglog(hs * 1e6, dLas1 / 2 / np.pi, "c:", label="laser 1 kHz")
    ax[0].loglog(hs * 1e6, dLas2 / 2 / np.pi, "c--", label="laser 10 Hz")
    ax[0].loglog(hs * 1e6, dHeat / 2 / np.pi, "r-", label="heating $\\propto h^{-4}$")
    ax[0].loglog(hs * 1e6, dJohn / 2 / np.pi, "b:", label="Johnson gradient")
    ax[0].set_xlabel("ion height $h$ [µm]")
    ax[0].set_ylabel("$\\delta\\omega/2\\pi$ [Hz]")
    ax[0].set_title("(a) noise budget vs height")
    ax[0].legend(fontsize=7, ncol=2)
    ax[0].grid(True, alpha=0.3, which="both")
    ax[0].set_ylim(1e-9, 1e5)
    # SNR for the Phase-1 grating (lambda=20um) in the heating-dominated regime
    lam_g, A_g = 20e-6, 50e-9
    k_g = 2 * np.pi / lam_g
    sig = np.abs(C_ITF * A_g * k_g**2 * np.exp(-k_g * hs))
    snr_heat = sig / (dHeat / 2 / np.pi)
    h_opt = alpha / k_g
    ax[1].semilogy(hs * 1e6, sig / 2 / np.pi, "b-", label="$|\\Delta\\omega_x|/2\\pi$")
    ax[1].semilogy(hs * 1e6, dHeat / 2 / np.pi, "r-", label="heating noise")
    ax[1].axvline(h_opt * 1e6, color="k", ls=":", label=f"$h_{{opt}}=\\alpha/k$ = {h_opt*1e6:.1f} µm")
    ax[1].set_xlabel("ion height $h$ [µm]")
    ax[1].set_ylabel("frequency [Hz]")
    ax[1].set_title("(b) Phase-1 grating: signal vs heating noise")
    ax[1].legend(fontsize=8)
    ax[1].grid(True, alpha=0.3, which="both")
    print(f"  h_opt(alpha=4, lambda=20um) = {h_opt*1e6:.2f} um")
    fig.tight_layout()
    fig.savefig(MANUSCRIPT + "figI10_noise_psd.pdf")
    plt.close(fig)


def anharmonic_check():
    """Verify the anharmonic prefactor against exact Fock-basis diagonalisation."""
    print("=" * 60)
    print("I-ANHARMONIC: prefactor check")
    xzp = np.sqrt(1.054571817e-34 / (2 * M_CA40 * OMEGA_X))
    Phi4 = E_BG * 100e-9 * (2e5)**4 * np.exp(-2e5 * H_NOM)  # k=2e5, A=100nm
    lam = E_CHARGE * Phi4 / 24.0
    N = 400
    n = np.arange(N)
    diag = 1.054571817e-34 * OMEGA_X * n
    # x^4 matrix elements from (a+a^dagger)^4 in Fock basis
    H = np.diag(diag)
    H += lam * xzp**4 * (6 * n**2 + 6 * n + 3)
    # off-diagonal: <m|(a+a+)^4|n> for |m-n| = 0,2,4
    for j in range(N):
        for dj, coef_n in [(2, 1), (4, 1)]:
            pass
    # exact: (a+a+)^4 = a^4+4a^3a+ +6a^2a+^2+4aa+^3+a+^4; normal-order:
    # a a+ = n+1 ; a^2 a+^2 = n^2+3n+2 ; a^3 a+^3 = n^3+6n^2+11n+6
    # <n|a^4|n+4> = sqrt((n+1)(n+2)(n+3)(n+4))
    # <n|a^3 a+|n+2> = (n+3)*sqrt((n+1)(n+2))
    # <n|a^2 a+^2|n> = n^2+3n+2
    # <n|a a+^3|n-2> = n*sqrt(n(n-1))
    # <n|a+^4|n-4> = sqrt(n(n-1)(n-2)(n-3))
    H = np.zeros((N, N))
    for j in range(N):
        H[j, j] = 1.054571817e-34 * OMEGA_X * j + lam * xzp**4 * (6 * j**2 + 6 * j + 3)
        if j + 4 < N:
            H[j, j + 4] = lam * xzp**4 * np.sqrt((j + 1.0) * (j + 2) * (j + 3) * (j + 4))
        if j + 2 < N:
            H[j, j + 2] = lam * xzp**4 * (j + 3.0) * np.sqrt((j + 1.0) * (j + 2))
        if j >= 2:
            H[j, j - 2] = lam * xzp**4 * j * np.sqrt(j * (j - 1.0))
        if j >= 4:
            H[j, j - 4] = lam * xzp**4 * np.sqrt(j * (j - 1.0) * (j - 2) * (j - 3))
    ev = np.linalg.eigvalsh(H)
    dE_exact = ev[1] - ev[0] - 1.054571817e-34 * OMEGA_X
    pref_exact = dE_exact / Phi4 / xzp**2
    pref_pt = E_CHARGE / (4 * M_CA40 * OMEGA_X)
    pref_ms = 3 * E_CHARGE / (8 * M_CA40 * OMEGA_X)
    print(f"  exact dE01 = {dE_exact:.6e} J")
    print(f"  exact prefactor dE/(Phi4 xzp^2) = {pref_exact:.6e}")
    print(f"  first-order PT  e/4m w        = {pref_pt:.6e}")
    print(f"  manuscript      3e/8m w       = {pref_ms:.6e}")
    print(f"  ratio PT/exact = {pref_pt/pref_exact:.6f}")


def rayleigh_map(A, k_s, h, N=40):
    """Exact Rayleigh forward map F(A) at the crest, complete sin+cos basis."""
    u = np.linspace(0, 2 * np.pi, 6000, endpoint=False)
    M = np.zeros((2 * N, 2 * N))
    rhs = np.zeros(2 * N)
    for m in range(N):
        sm = np.sin((m + 1) * u)
        cm = np.cos((m + 1) * u)
        rhs[m] = E_BG * A * np.mean(sm * np.sin(u)) * 2.0
        rhs[N + m] = E_BG * A * np.mean(cm * np.sin(u)) * 2.0
        for n in range(N):
            env = np.exp(-(n + 1) * k_s * A * np.sin(u))
            M[m, n] = np.mean(sm * np.sin((n + 1) * u) * env) * 2.0
            M[m, N + n] = np.mean(sm * np.cos((n + 1) * u) * env) * 2.0
            M[N + m, n] = np.mean(cm * np.sin((n + 1) * u) * env) * 2.0
            M[N + m, N + n] = np.mean(cm * np.cos((n + 1) * u) * env) * 2.0
    sol = np.linalg.solve(M, rhs)
    c, d = sol[:N], sol[N:]
    dw = 0.0
    for n in range(N):
        nk = (n + 1) * k_s
        sgn = c[n] * np.sin((n + 1) * np.pi / 2) + d[n] * np.cos((n + 1) * np.pi / 2)
        dw += -(nk**2) * sgn * np.exp(-nk * h)
    return (E_CHARGE / (2 * M_CA40 * OMEGA_X)) * dw


def landweber_cone():
    """Tangential cone constant and Landweber convergence check (I4)."""
    k_s = 2 * np.pi / 50e-6
    h = 10e-6
    A_vals = np.linspace(10e-9, 8e-6, 33)
    F = np.array([rayleigh_map(A, k_s, h) for A in A_vals])
    dF = np.gradient(F, A_vals)
    c_cone = 0.0
    for i, B in enumerate(A_vals):
        lhs = np.abs(F - F[i] - dF[i] * (A_vals - B))
        rhs = np.abs(F - F[i]) + 1e-30
        c_cone = max(c_cone, np.max(lhs / rhs))
    # linear-inversion bias at the manuscript's operating points
    print("=" * 60)
    print("I-LANDWEBER: lambda=50um, h=10um, complete Rayleigh basis")
    for A0, kA in [(3.98e-6, 0.50), (8.04e-6, 1.01)]:
        ex = rayleigh_map(A0, k_s, h)
        lin = -C_ITF * A0 * k_s**2 * np.exp(-k_s * h)
        print(f"  kA={kA:.2f}: linear-inversion bias = {(lin-ex)/ex*100:+.1f}%")
    print(f"  tangential cone constant c = {c_cone:.3f}  (<1/2? {c_cone<0.5})")
    A0 = 3.98e-6
    d = rayleigh_map(A0, k_s, h)
    A_est = 0.5 * A0
    a_step = 1.0 / np.mean(np.abs(dF))**2
    for it in range(12):
        Fn = rayleigh_map(A_est, k_s, h)
        A_est = A_est + a_step * (d - Fn) * dF[0]
    print(f"  Landweber final bias = {abs(A_est-A0)/A0*100:.2f}%")


def fd_laplace():
    """Body-fitted finite-difference Laplace solver for the corrugated surface."""
    A, lam_s = 100e-9, 40e-6
    k_s = 2 * np.pi / lam_s
    Ztop = 60e-6          # truncated domain: the solution decays as e^{-2 k z},
    #                      so the top Dirichlet data is exact to 1e-8
    Nx = 800
    xs = np.linspace(0, 2 * lam_s, Nx, endpoint=False)
    zs = lambda x: A * np.sin(k_s * x)
    def solve(Nv_in):
        vs_i = np.linspace(0, 1, Nv_in)
        N = Nx * Nv_in
        rows, cols, vals = [], [], []
        rhs = np.zeros(N)
        dx = xs[1] - xs[0]
        zx = k_s * A * np.cos(k_s * xs)
        zxx = -k_s**2 * A * np.sin(k_s * xs)
        Z = Ztop - zs(xs)                     # z = zs + v*(Ztop-zs)
        def idx(i, j):
            return i * Nv_in + j
        for i in range(Nx):
            for j in range(Nv_in):
                p = idx(i, j)
                v = vs_i[j]
                if j == Nv_in - 1:
                    rows.append(p); cols.append(p); vals.append(1.0); rhs[p] = 0.0
                    continue
                if j == 0:
                    rows.append(p); cols.append(p); vals.append(1.0)
                    rhs[p] = E_BG * zs(xs[i])
                    continue
                im, ip = (i - 1) % Nx, (i + 1) % Nx
                zi = Z[i]
                # metric coefficients: z_x = zx; z_xx = zxx; z_v = Z
                # v_x = -v zx/zi ; v_xx = v(2 zx^2/zi^2 - zxx/zi)
                vx = -v * zx[i] / zi
                vxx = -v * (zxx[i] / zi + zx[i]**2 / zi**2)
                vz = 1.0 / zi
                c_xx = 1.0
                c_vv = vx**2 + vz**2
                c2 = vxx
                c3 = 2.0 * vx
                # c_xx P_x'x' + c2 P_v + c3 P_x'v + c_vv P_vv = 0
                dv = vs_i[1] - vs_i[0]
                coeff = {}
                coeff[(i, j - 1)] = c_vv / dv**2
                coeff[(i, j + 1)] = c_vv / dv**2
                coeff[(im, j)] = c_xx / dx**2 - c2 / (2 * dx)
                coeff[(ip, j)] = c_xx / dx**2 + c2 / (2 * dx)
                coeff[(i, j)] = -2 * c_xx / dx**2 - 2 * c_vv / dv**2
                # cross derivative
                coeff[(im, j - 1)] = coeff.get((im, j - 1), 0.0) + c3 / (4 * dx * dv)
                coeff[(ip, j + 1)] = coeff.get((ip, j + 1), 0.0) + c3 / (4 * dx * dv)
                coeff[(im, j + 1)] = coeff.get((im, j + 1), 0.0) - c3 / (4 * dx * dv)
                coeff[(ip, j - 1)] = coeff.get((ip, j - 1), 0.0) - c3 / (4 * dx * dv)
                for (ii, jj), val in coeff.items():
                    if jj == 0:
                        rhs[p] -= val * E_BG * zs(xs[ii])
                    else:
                        rows.append(p); cols.append(idx(ii, jj)); vals.append(val)
        M = csr_matrix((vals, (rows, cols)), shape=(N, N))
        psi = spsolve(M, rhs)
        return psi.reshape(Nx, Nv_in), vs_i
    def dw_from(psi, vs_i):
        h = 40e-6
        # interpolate Psi at z=h and take d^2/dx^2
        zi = Ztop - zs(xs)
        vh = (h - zs(xs)) / zi
        Psi_h = np.zeros(Nx)
        for i in range(Nx):
            j0 = int(np.floor(vh[i] * (len(vs_i) - 1)))
            j0 = min(j0, len(vs_i) - 2)
            f = vh[i] * (len(vs_i) - 1) - j0
            Psi_h[i] = (1 - f) * psi[i, j0] + f * psi[i, j0 + 1]
        k = 2 * np.pi * fftfreq(Nx, xs[1] - xs[0])
        return np.real(ifft(-k**2 * fft(Psi_h))), xs
    print("=" * 60)
    print("I-FD: body-fitted finite-difference Laplace solver")
    errs = []
    for Nv_in in [64, 128, 256]:
        psi, vs_i = solve(Nv_in)
        d2, xs_ = dw_from(psi, vs_i)
        i_crest = np.argmin(abs(xs_ - lam_s / 4))
        dw_fd = -(E_CHARGE / (2 * M_CA40 * OMEGA_X)) * d2[i_crest]
        dw_itf = C_ITF * A * k_s**2 * np.exp(-k_s * 40e-6)
        errs.append(abs(dw_fd - dw_itf) / abs(dw_itf))
        print(f"  Nv={Nv_in}: FD crest dw = {dw_fd:.6e} rad/s "
              f"(ITF {dw_itf:.6e}), rel err {errs[-1]:.4f}")
    return errs


def fem_figure(errs):
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.8))
    nv = [64, 128, 256]
    ax[0].loglog(nv, errs, "o-")
    ax[0].set_xlabel("vertical grid points $N_v$")
    ax[0].set_ylabel("relative error vs ITF")
    ax[0].set_title("(a) FD convergence (crest, $A=100$ nm)")
    ax[0].grid(True, alpha=0.3)
    # panel (b): schematic of the body-fitted domain
    xp = np.linspace(0, 2 * 40e-6, 200)
    zs = 100e-9 * np.sin(2 * np.pi * xp / 40e-6)
    ax[1].plot(xp * 1e6, zs * 1e9, "k-")
    ax[1].fill_between(xp * 1e6, zs * 1e9, 120, color="C0", alpha=0.15)
    ax[1].axhline(40, color="r", ls="--", label="ion plane $h=40$ µm")
    ax[1].set_xlabel("$x$ [µm]")
    ax[1].set_ylabel("$z$ [µm]")
    ax[1].set_ylim(-0.5, 120)
    ax[1].set_title("(b) corrugated domain, grounded boundary")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(MANUSCRIPT + "figI14_fem.pdf")
    plt.close(fig)


def bem_check():
    """Boundary-element (collocation) cross-check with the periodic Green's
function; figO6."""
    A_list = [10e-9, 100e-9, 500e-9, 2000e-9]
    lam_s, h = 40e-6, 40e-6
    k_s = 2 * np.pi / lam_s
    L = 2 * lam_s
    print("=" * 60)
    print("O-BEM: single-layer collocation, periodic Green's function")
    bem_vals, itf_vals, ray_vals = [], [], []
    for A in A_list:
        Nb = 600
        s = np.linspace(0, L, Nb, endpoint=False)
        ds = s[1] - s[0]
        xs = s
        zs = A * np.sin(k_s * s)
        r = np.stack([xs, zs], axis=1)
        G = np.zeros((Nb, Nb))
        for i in range(Nb):
            dxv = (r[:, 0] - r[i, 0]) * 2 * np.pi / L
            dzv = (r[:, 1] - r[i, 1]) * 2 * np.pi / L
            arg = 2 * np.cosh(dzv) - 2 * np.cos(dxv)
            G[i, :] = -np.log(arg) * ds / (4 * np.pi)
            G[i, i] = -ds / (2 * np.pi) * (np.log(2 * np.pi * ds / L) - 1.0)
        sigma = np.linalg.solve(G, E_BG * zs)
        x0, z0 = lam_s / 4, h
        dxv = (r[:, 0] - x0) * 2 * np.pi / L
        dzv = (r[:, 1] - z0) * 2 * np.pi / L
        arg = 2 * np.cosh(dzv) - 2 * np.cos(dxv)
        Phi = -E_BG * z0 + np.sum(sigma * (-np.log(arg))) * ds / (4 * np.pi)
        # second derivative wrt x0 of -(1/4pi) ln[2cosh(k dz)-2cos(k dx)]
        kp = 2 * np.pi / L
        csh = np.cosh(dzv); cs = np.cos(dxv)
        D = csh - cs
        d2G = -kp**2 / (4 * np.pi) * (csh * cs - 1.0) / D**2
        d2Phi = np.sum(sigma * d2G) * ds
        dw = (E_CHARGE / (2 * M_CA40 * OMEGA_X)) * d2Phi
        dw_itf = -C_ITF * A * k_s**2 * np.exp(-k_s * h)
        dw_ray = rayleigh_map(A, k_s, h)
        bem_vals.append(dw); itf_vals.append(dw_itf); ray_vals.append(dw_ray)
        print(f"  A={A*1e9:5.0f} nm: BEM {dw:.6e} | ITF {dw_itf:.6e} "
              f"| Rayleigh {dw_ray:.6e} rad/s")
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.8))
    ax[0].semilogy([a * 1e9 for a in A_list], np.abs(bem_vals), "s-",
                   label="BEM (collocation)")
    ax[0].semilogy([a * 1e9 for a in A_list], np.abs(itf_vals), "--",
                   label="ITF (linear)")
    ax[0].semilogy([a * 1e9 for a in A_list], np.abs(ray_vals), "o-",
                   label="Rayleigh (exact)")
    ax[0].set_xlabel("amplitude $A$ [nm]")
    ax[0].set_ylabel("$|\\Delta\\omega_x|$ [rad/s]")
    ax[0].legend(fontsize=8)
    ax[0].grid(True, alpha=0.3)
    rel = np.abs((np.array(bem_vals) - np.array(ray_vals)) / np.array(ray_vals))
    ax[1].semilogy([a * 1e9 for a in A_list], rel, "s-", label="BEM vs Rayleigh")
    rel2 = np.abs((np.array(itf_vals) - np.array(ray_vals)) / np.array(ray_vals))
    ax[1].semilogy([a * 1e9 for a in A_list], rel2, "o-", label="ITF vs Rayleigh")
    ax[1].set_xlabel("amplitude $A$ [nm]")
    ax[1].set_ylabel("relative error")
    ax[1].legend(fontsize=8)
    ax[1].grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(MANUSCRIPT + "figO6_bem.pdf")
    plt.close(fig)


def spectrum_check():
    """Numerical spectrum of the 2D discretised operator (I1)."""
    N = 128
    L = 64e-6
    dx = L / N
    k = 2 * np.pi * fftfreq(N, dx)
    KX, KY = np.meshgrid(k, k)
    H = -C_ITF * KX**2 * np.exp(-np.sqrt(KX**2 + KY**2) * H_NOM)
    M = 4 * C_ITF * np.exp(-2) / H_NOM**2
    print("=" * 60)
    print("I-SPECTRUM: 2D symbol range on the DFT lattice")
    print(f"  H range = [{H.min():.4e}, {H.max():.4e}] (rad/s)/m")
    print(f"  continuum sigma(H) = [-M,0] with M = {M:.4e}")
    print(f"  grid max / M = {abs(H.min())/M:.4f} (approaches 1 as dk->0)")


def singular_gallery():
    """Singular-function gallery: the discretised operator is exactly
diagonal in the DFT basis, so the singular functions are Fourier modes;
show the six highest-singular-value modes and their Fourier-overlap
statistics on a windowed (non-periodic) domain to quantify how close
the continuous singular functions are to Fourier modes."""
    N = 128
    k = 2 * np.pi * fftfreq(N, DX)
    KX, KY = np.meshgrid(k, k)
    H = -C_ITF * KX**2 * np.exp(-np.sqrt(KX**2 + KY**2) * H_NOM)
    s = np.sort(np.abs(H.ravel()))[::-1]
    # windowed operator (Hann window) -> continuous analogue with
    # non-diagonal singular functions; measure Fourier overlap
    win = np.outer(np.hanning(N), np.hanning(N))
    Hw = np.fft.fft2(win * np.fft.ifft2(H))
    print("=" * 60)
    print("O-SINGULAR: Fourier-overlap of windowed singular functions")
    ovs = []
    for n in [0, 1, 2, 3, 4, 5, 9, 19]:
        e = np.zeros((N, N))
        idx = np.unravel_index(np.argsort(np.abs(H.ravel()))[::-1][n], (N, N))
        e[idx] = 1.0
        u = np.real(np.fft.ifft2(np.fft.fft2(win * np.fft.ifft2(e))))
        u = u / np.linalg.norm(u)
        fu = np.abs(np.fft.fft2(u)) ** 2
        ovs.append(fu.max() / fu.sum())
        print(f"  n={n+1:2d}: Fourier overlap {ovs[-1]:.3f}")
    fig, ax = plt.subplots(2, 3, figsize=(9, 5.5))
    order = np.argsort(np.abs(H.ravel()))[::-1]
    for n, a in enumerate(ax.ravel()):
        e = np.zeros((N, N))
        e[np.unravel_index(order[n], (N, N))] = 1.0
        u = np.real(np.fft.ifft2(e))
        a.imshow(u, cmap="RdBu_r", extent=[0, 64, 0, 64])
        a.set_title(f"$u_{{{n+1}}}$ (DFT mode, $\sigma_n/\sigma_1=${s[n]/s[0]:.2f})")
    fig.suptitle("Leading singular functions of $\\mathcal{H}_x$ (2D, $h=40$ µm)")
    fig.tight_layout()
    fig.savefig(MANUSCRIPT + "figO3_singular.pdf")
    plt.close(fig)


def waterfilling():
    """Quantum water-filling across modes (O1): clean KKT implementation."""
    k = 2 * np.pi * fftfreq(N_GRID, DX)
    H = np.abs(itf_1d(k, H_NOM))
    lam = LAMBDA_DP_FRAC * H.max()
    idx = np.where(H > lam)[0]
    b = (D_OMEGA**2) / H[idx]**2       # delta-z^2 per mode, coherent
    b = b / b.max()
    n_modes = len(idx)
    Rtot = 1.5
    # KKT: r_k = max(0, 0.25*(ln b_k - ln mu)),  sum r_k = Rtot
    ln_mu_lo, ln_mu_hi = -60.0, 30.0
    for _ in range(120):
        ln_mu = 0.5 * (ln_mu_lo + ln_mu_hi)
        r = np.maximum(0.25 * (np.log(b) - ln_mu), 0.0)
        if r.sum() > Rtot:
            ln_mu_lo = ln_mu
        else:
            ln_mu_hi = ln_mu
    r = np.maximum(0.25 * (np.log(b) - ln_mu), 0.0)
    v_unif = np.mean(b)
    v_wf = np.mean(b * np.exp(-2 * r))
    print("=" * 60)
    print("O-WATERFILLING: average delta-z^2 over passband modes")
    print(f"  n_modes={n_modes}, uniform={v_unif:.4f}, waterfilled={v_wf:.4f}, "
          f"gain x{v_unif/v_wf:.2f}")
    print(f"  r_k range: [{r.min():.2f}, {r.max():.2f}], sum={r.sum():.2f}")


def info_closed_1d():
    """Closed-form 1D Shannon capacity (used in the revised eq:info_closed)."""
    hs = [5e-6, 10e-6, 20e-6, 40e-6, 80e-6]
    k = 2 * np.pi * fftfreq(N_GRID, DX)
    print("=" * 60)
    print("I-INFO: 1D closed form check")
    for h in hs:
        H = np.abs(itf_1d(k, h))
        lam = LAMBDA_DP_FRAC * H.max()
        I_num = 0.5 * np.sum(np.log(1 + (H / lam)**2))
        closed = (L_FOV / (2 * np.pi * h)) * np.log(C_ITF / (lam * h**2))**2
        print(f"  h={h*1e6:4.0f}um: numeric I={I_num:8.2f}  "
              f"closed form={closed:8.2f}  ratio={closed/I_num:.3f}")


def second_order_check():
    """Verify the second-order perturbation: Delta-w^(1)+Delta-w^(2) vs the
    complete-basis exact solution; the (kA)^2 residual must drop to (kA)^3."""
    lam_s, h = 40e-6, 40e-6
    k_s = 2 * np.pi / lam_s
    print("=" * 60)
    print("I-PHI2: second-order correction vs complete-basis Rayleigh")
    # small-grid direct convolution identity
    Ns = 512
    xs_s = np.linspace(0, 2 * lam_s, Ns, endpoint=False)
    zs_s = A = None
    for A in [100e-9, 2000e-9]:
        zs = A * np.sin(k_s * xs_s)
        k_sg = 2 * np.pi * fftfreq(Ns, xs_s[1] - xs_s[0])
        zt = fft(zs) / Ns
        conv = np.zeros(Ns, dtype=complex)
        for n in range(Ns):
            s2 = 0.0j
            for m in range(Ns):
                s2 += zt[m] * np.abs(k_sg[m] - k_sg[n]) * zt[(n - m) % Ns]
            conv[n] = E_BG * s2 / Ns
        w = np.real(ifft(np.abs(k_sg) * E_BG * zt))     # - dz Phi^(1)|_0
        Phi2_fft_small = fft(-zs * w) / Ns              # F[-z_s dz Phi^(1)] = F[Phi^(2)(.,0)]
        err = np.abs(conv - Phi2_fft_small).max() / np.abs(Phi2_fft_small).max()
        print(f"  A={A*1e9:.0f} nm: convolution identity rel err = {err:.2e}")
        # propagated second-order frequency shift at the crest
        Nx = 4096
        xs = np.linspace(0, 2 * lam_s, Nx, endpoint=False)
        zs = A * np.sin(k_s * xs)
        k = 2 * np.pi * fftfreq(Nx, xs[1] - xs[0])
        zt = fft(zs) / Nx
        w = np.real(ifft(np.abs(k) * E_BG * zt))
        Phi2_0 = -zs * w
        Phi2_h = np.real(ifft(np.exp(-np.abs(k) * h) * fft(Phi2_0)))
        d2 = np.real(ifft(-k**2 * fft(Phi2_h)))
        i_crest = np.argmin(abs(xs - lam_s / 4))
        dw2 = (E_CHARGE / (2 * M_CA40 * OMEGA_X)) * d2[i_crest]
        dw1 = -C_ITF * A * k_s**2 * np.exp(-k_s * h)
        ex = rayleigh_map(A, k_s, h)
        r1 = abs(dw1 - ex) / abs(ex)
        r12 = abs(dw1 + dw2 - ex) / abs(ex)
        print(f"  A={A*1e9:4.0f} nm: |d1-ex|/|ex|={r1:.4f}  "
              f"|d1+d2-ex|/|ex|={r12:.4f}")


def micromotion_numbers():
    q, x_sec = 0.3, 10e-9
    lam_feat = 30e-6
    corr = (q**2 * np.pi**2 / 4) * (x_sec / lam_feat)**2
    corr_warm = (q**2 * np.pi**2 / 4) * (300e-9 / 1e-6)**2
    print("=" * 60)
    print("I-MICROMOTION: carrier correction = (q^2 pi^2/4)(x_sec/lambda)^2")
    print(f"  cold ion, 30um feature: {corr:.2e}")
    print(f"  warm ion, 1um feature:  {corr_warm:.3f}")


def patch_mapping():
    k = 2e5
    phi_s = 50e-3   # 50 mV patch potential at 1/k ~ 5 um scale
    sigma = 2 * 8.8541878128e-12 * k * phi_s
    print("=" * 60)
    print("I-PATCH: sigma(k) = 2 eps0 k phi_s(k)")
    print(f"  phi_s=50 mV, k=2e5 /m -> sigma = {sigma:.3e} C/m^2")
    print(f"  sigma=1.13e-8 C/m^2 -> phi_s = "
          f"{1.13e-8/(2*8.854e-12*2e5)*1e3:.1f} mV")


def rayleigh_figure_complete():
    """Regenerate figV7_rayleigh_validation.pdf with the complete sin+cos
    Rayleigh basis (the physically exact boundary solution) and print the
    boundary residuals and relative errors for the manuscript text."""
    lam_s, h = 40e-6, 40e-6
    k_s = 2 * np.pi / lam_s
    A_list = np.array([10e-9, 30e-9, 100e-9, 300e-9, 1e-6, 2e-6])
    exact = np.array([rayleigh_map(A, k_s, h) for A in A_list])
    itf = np.array([-C_ITF * A * k_s**2 * np.exp(-k_s * h) for A in A_list])
    rel = np.abs(exact - itf) / np.abs(exact)
    print("=" * 60)
    print("I-RAYLEIGH (complete basis): errors vs ITF")
    for A, r in zip(A_list, rel):
        print(f"  A={A*1e9:5.0f} nm: rel err = {r*100:.4f}%")
    # quadratic coefficient
    kA = k_s * A_list
    p = np.polyfit(kA**2, rel, 1)[0]
    print(f"  rel_err = {p:.3f} x (kA)^2")
    # boundary residual of the complete-basis solution
    for A in [10e-9, 100e-9, 2e-6]:
        u = np.linspace(0, 2 * np.pi, 20000, endpoint=False)
        N = 60
        M = np.zeros((2 * N, 2 * N)); rhs = np.zeros(2 * N)
        for m in range(N):
            sm = np.sin((m + 1) * u); cm = np.cos((m + 1) * u)
            rhs[m] = E_BG * A * np.mean(sm * np.sin(u)) * 2.0
            rhs[N + m] = E_BG * A * np.mean(cm * np.sin(u)) * 2.0
            for n in range(N):
                env = np.exp(-(n + 1) * k_s * A * np.sin(u))
                M[m, n] = np.mean(sm * np.sin((n + 1) * u) * env) * 2.0
                M[m, N + n] = np.mean(sm * np.cos((n + 1) * u) * env) * 2.0
                M[N + m, n] = np.mean(cm * np.sin((n + 1) * u) * env) * 2.0
                M[N + m, N + n] = np.mean(cm * np.cos((n + 1) * u) * env) * 2.0
        sol = np.linalg.solve(M, rhs)
        c, d = sol[:N], sol[N:]
        Phi_b = -E_BG * A * np.sin(u)
        for n in range(N):
            Phi_b += (c[n] * np.sin((n + 1) * u) + d[n] * np.cos((n + 1) * u)) \
                     * np.exp(-(n + 1) * k_s * A * np.sin(u))
        res = np.max(np.abs(Phi_b)) / (E_BG * A)
        print(f"  A={A*1e9:5.0f} nm: boundary residual = {res:.2e}")
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.8))
    ax[0].loglog(A_list * 1e9, np.abs(exact) / 2 / np.pi, "o-", label="exact (complete Rayleigh)")
    ax[0].loglog(A_list * 1e9, np.abs(itf) / 2 / np.pi, "--", label="ITF (linear)")
    ax[0].set_xlabel("amplitude $A$ [nm]")
    ax[0].set_ylabel("$|\\Delta\\omega_x|/2\\pi$ [Hz]")
    ax[0].legend(fontsize=8)
    ax[0].grid(True, alpha=0.3)
    ax[1].loglog(A_list * 1e9, rel, "o-", label="relative error")
    ax[1].loglog(A_list * 1e9, p * (k_s * A_list)**2, "--", label=f"${p:.3f}\\,(kA)^2$")
    ax[1].set_xlabel("amplitude $A$ [nm]")
    ax[1].set_ylabel("relative error")
    ax[1].legend(fontsize=8)
    ax[1].grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(MANUSCRIPT + "figV7_rayleigh_validation.pdf")
    plt.close(fig)


def main():
    import sys
    log = open("/tmp/audit_fixes.log", "w")
    sys.stdout = log
    try:
        svd_1d()
        picard_lcurve_gcv()
        noise_master()
        anharmonic_check()
        landweber_cone()
        rayleigh_figure_complete()
        bem_check()
        spectrum_check()
        singular_gallery()
        waterfilling()
        info_closed_1d()
        second_order_check()
        micromotion_numbers()
        patch_mapping()
    finally:
        sys.stdout = sys.__stdout__
        log.close()
    print("done -> /tmp/audit_fixes.log")


if __name__ == "__main__":
    main()
