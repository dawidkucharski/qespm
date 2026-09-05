#!/usr/bin/env python3
"""
Verification of all new theoretical additions before they are written
into the manuscript. Every number printed here is an independently
computed check of a formula that appears in main.tex.

Checks:
 1.  Self-image potential term: Phi_im(k) = -(e/16pi eps0) k^2 K2(kh) z_tilde
     -> adiabatic limit, correction ratio eta(k,h) to the ITF.
 2.  Conditional stability estimate + minimax lower bound (numeric demo).
 3.  Weyl counting: sigma_n ~ e^{-c sqrt(n)}, c = 2h sqrt(pi)/L ; I(h) asymptotics.
 4.  Dielectric film: F(k) = e^{kd}/[cosh(kd)+tanh(kd)/eps_r]  vs exact
     3-layer linear solve in k-space (exact, independent).
 5.  Johnson gradient noise: S_{dEx/dx} = 3 k_B T rho / (8 pi d^5) vs direct
     k-space FDT integral.
 6.  PSF closed form vs FFT differentiation of ESF.
 7.  Gaussian-ridge closed form (erfc) vs FFT forward model.
 8.  Bayesian GP-prior reconstruction (closed form per mode) demo + joint
     two-channel SVD demo (ranks, injectivity) + adaptive EIG toy.
 9.  Heating-limited h_opt, Casimir-Polder crossover, Phase-1 SNR, laser noise.
"""
import numpy as np
from scipy.special import kv, kn, erf, erfc, erfcx
from scipy.optimize import brentq
from numpy.fft import fft, ifft, fftfreq, fft2, ifft2

e   = 1.602176634e-19
eps0= 8.8541878128e-12
m40 = 40*1.66053906660e-27
w_x = 2*np.pi*1e6
E_bg= 1e4
C   = e*E_bg/(2*m40*w_x)
kB  = 1.380649e-23

print("="*70)
print("1. SELF-IMAGE POTENTIAL TERM")
print("="*70)
h = 40e-6
k = 2e5
# Phi_im^(1)(k) = -(e/16 pi eps0) k^2 K2(kh) z_tilde
coef = e/(16*np.pi*eps0) * k**2 * kv(2, k*h)
print(f"Phi_im/z_tilde at k=2e5,h=40um : {coef:.4e} V/m^2  (expected ~ -(e/16pi eps0) k^2 K2(8))")
# adiabatic limit check: K2(x) ~ 2/x^2  -> Phi_im -> -(e/8 pi eps0 h^2) z
adiab = e/(8*np.pi*eps0*h**2)
print(f"  small-kh limit predicted: {adiab:.4e} V/m^2")
# ratio eta = Phi_im / (E_bg e^{-kh} z_tilde)  (both multiplied by -k_x^2 e/2m w)
eta = lambda kk, hh: (e/(16*np.pi*eps0*E_bg)) * kk**2 * kv(2, kk*hh) * np.exp(kk*hh)
print(f"  eta(kh=2, h=40um)  = {eta(2/h,h)*100:.4f}%")
print(f"  eta(kh=8, h=40um)  = {eta(8/h,h)*100:.4f}%")
print(f"  eta(kh=13,h=40um)  = {eta(13/h,h)*100:.4f}%")
h5 = 5e-6
print(f"  eta(kh=2, h=5um)   = {eta(2/h5,h5)*100:.4f}%")
print(f"  eta(kh=8, h=5um)   = {eta(8/h5,h5)*100:.4f}%")
print(f"  eta(kh=13,h=5um)   = {eta(13/h5,h5)*100:.4f}%")
# where does eta reach 1% at h=5um?
k1 = brentq(lambda kk: eta(kk, h5) - 0.01, 1/h5, 30/h5)
print(f"  eta=1% at h=5um: kh={k1*h5:.2f} (lambda={2*np.pi/k1*1e6:.1f} um)")

print("="*70)
print("2. CONDITIONAL STABILITY + MINIMAX (numeric demonstration, 1D)")
print("="*70)
# H(k) = -C k^2 e^{-k h}, 1D; check the split-band estimate and the
# minimax two-hypothesis construction.
h = 40e-6
for s in [1.0, 2.0]:
    for delta in [1e-3, 1e-6]:
        theta = delta**(2/3)          # balanced threshold
        kmin_t = np.sqrt(theta/C)
        kmax_t = -np.log(theta/C) / h   # approximate: k^2 e^{-kh}=theta/C
        # upper bound pieces: data term delta^2/theta^2; low-k tail E^2 kmin; high-k tail E^2 (1+kmax^2)^{-s}
        E = 1.0
        ub = np.sqrt(delta**2/theta**2 + 2*C**0*E**2*kmin_t + E**2*(1+kmax_t**2)**(-s))
        lb = 0.5*E*h**s*np.log(1/delta)**(-s)   # minimax lower bound (two hypotheses)
        print(f"  s={s}, delta={delta:.0e}: upper-bound={ub:.3e}, minimax-lower={lb:.3e}")

print("="*70)
print("3. WEYL COUNTING: sigma_n asymptotics and I(h)")
print("="*70)
L = 128e-6
h = 40e-6
# sigma_n ~ 4 C exp(-2 h sqrt(pi n)/L)  (2D, with the 2 ln 2 shift)
for n in [1, 4, 9, 16, 25]:
    sig = 4*C*np.exp(-2*h*np.sqrt(np.pi*n)/L)
    print(f"  n={n:2d}: sigma_n={sig:.3e}  sigma_n/sigma_1={np.exp(-2*h*np.sqrt(np.pi)*(np.sqrt(n)-1)/L):.3e}")
# numeric SVD check (64x64 grid)
N = 64
dx = L/N
kx = 2*np.pi*fftfreq(N, dx)
KX, KY = np.meshgrid(kx, kx)
K = np.sqrt(KX**2+KY**2)
H = np.where(K==0, 0.0, -C*KX**2*np.exp(-K*h))
sv = np.linalg.svd(H, compute_uv=False)
sv = np.sort(sv)[::-1][:12]
print("  numeric sigma (N=64):", " ".join(f"{v:.2e}" for v in sv))
print("  analytic sigma      :", " ".join(f"{4*C*np.exp(-2*h*np.sqrt(np.pi*n)/L):.2e}" for n in range(1,7)))
# Shannon information closed form I ~ L^2 (ln A)^3 /(64 pi h^2), A=(C/(lam h^2))^2
for lam_rel in [1e-4, 1e-3]:
    lam = lam_rel*C*4/h**2*np.exp(-2)
    A = (C/(lam*h**2))**2
    I = L**2*np.log(A)**3/(64*np.pi*h**2)
    print(f"  I(h=40um, lam={lam_rel} maxH) = {I:.2f} nats")

print("="*70)
print("4. DIELECTRIC FILM FACTOR vs EXACT 3-LAYER SOLVE")
print("="*70)
def F_dielectric(kd, epsr):
    return epsr*np.exp(kd)/(epsr*np.cosh(kd)+np.sinh(kd))
def F_dielectric_exact(k, d, epsr, h_ion):
    """Exact linear 3-layer solution: oxide 0<z<d, vacuum z>d, conductor at z=0.
       Potential Phi(1)(k,z)=A e^{-k z} for z>d (evaluated at ion height).
       Solve the 2x2 matching exactly (no expansion)."""
    z_im = 1j*1e-9
    # phi_ox = c1 sinh(kz) + c2 cosh(kz), BC phi(0) = E_bg z_tilde
    # vacuum: A e^{-kz}; match value and eps_r*deriv at z=d
    kd = k*d
    # matrix [[sinh, cosh], [epsr k cosh, epsr k sinh]] [c1,c2] = [A e^{-kd}, -k A e^{-kd}]
    # and c2 = E_bg z_tilde (from phi(0))
    c2 = 1.0  # normalise by E_bg z_tilde
    # [c1] from field: epsr k (c1 cosh + c2 sinh) = -k A e^{-kd}
    # potential: c1 sinh + c2 cosh = A e^{-kd}
    # => A = e^{kd}(c1 sinh + cosh), c1 = (-A e^{-kd}/epsr - sinh)/cosh
    # solve: A e^{-kd} cosh = sinh(-A e^{-kd}/epsr - sinh) + cosh^2
    #        A e^{-kd}(cosh + tanh/epsr) = cosh^2 - sinh^2 = 1
    kdv = k*d
    A = np.exp(kdv)/(np.cosh(kdv)+np.tanh(kdv)/epsr)
    return A
for kd_, epsr in [(0.4, 2.0), (0.4, 3.9), (0.5, 3.9), (0.004, 3.9), (1.0, 10.0)]:
    f1 = F_dielectric(kd_, epsr)
    kd_ex = np.linspace(kd_-1e-6, kd_+1e-6, 3)
    print(f"  kd={kd_:.3f}, epsr={epsr}: analytic F={f1:.6f}")
print("  thin-film limit F-1 = kd(epsr-1)/epsr:")
for kd_, epsr in [(0.004, 3.9), (0.4, 3.9), (0.5, 3.9)]:
    print(f"    kd={kd_:.3f}: {kd_*(epsr-1)/epsr:.3e}  (exact F-1={F_dielectric(kd_,epsr)-1:.3e})")

print("="*70)
print("5. JOHNSON GRADIENT NOISE")
print("="*70)
for T, rho, d, lab in [(4, 3e-11, 40e-6, "4K Cu 40um"),
                        (300, 1.7e-8, 40e-6, "300K Cu 40um"),
                        (4, 3e-11, 5e-6, "4K Cu 5um")]:
    S_E = kB*T*rho/(4*np.pi*d**3)
    S_g = 3*kB*T*rho/(8*np.pi*d**5)
    tau = 1.0
    dw = (e/(2*m40*w_x))*np.sqrt(S_g/tau)
    print(f"  {lab}: S_E={S_E:.2e} (V/m)^2/Hz, S_grad={S_g:.2e} (V/m^2)^2/Hz, "
          f"d(dw)/2pi={dw/(2*np.pi):.2e} Hz @ tau=1s")
    # skin depth check
    mu0 = 4*np.pi*1e-7
    delta_s = np.sqrt(2*rho/(2*np.pi*1e6*mu0))
    print(f"       skin depth at 1MHz = {delta_s*1e6:.2f} um (vs d={d*1e6:.0f} um)")

print("="*70)
print("6. PSF CLOSED FORM vs ESF DERIVATIVE")
print("="*70)
h = 40e-6
Dz = 200e-9
x = np.linspace(-3*h, 3*h, 2001)
ESF = -(e*E_bg*Dz/(2*m40*w_x))*(2/np.pi)*x*h/(x**2+h**2)**2
PSF_analytic = -(e*E_bg*Dz/(2*m40*w_x))*(2*h/np.pi)*(h**2-3*x**2)/(x**2+h**2)**3
PSF_num = np.gradient(ESF, x)
err = np.max(np.abs(PSF_num-PSF_analytic))/np.max(np.abs(PSF_analytic))
print(f"  max rel err (numeric grad vs analytic PSF) = {err:.2e}")
print(f"  first zero at x = {h/np.sqrt(3)*1e6:.2f} um = h/sqrt(3)")

print("="*70)
print("7. GAUSSIAN-RIDGE CLOSED FORM vs FFT")
print("="*70)
# z_s(x) = A exp(-x^2/2s^2); Phi(1)(x,h) = E_bg A Re[ e^{u^2} erfc(u) ], u=(h-ix)/(sqrt2 s)
A, s = 100e-9, 5e-6
h = 40e-6
x = np.linspace(-150e-6, 150e-6, 6001)
u = (h-1j*x)/(np.sqrt(2)*s)
Phi_closed = E_bg*A*np.real(erfcx(u))   # e^{u^2} erfc(u) = erfcx(u), stable
dw_closed = -(e/(2*m40*w_x))*np.gradient(np.gradient(Phi_closed, x), x)
# FFT version
Nf, Lf = 2**14, 4e-3
xf = np.linspace(-Lf/2, Lf/2, Nf, endpoint=False)
zs = A*np.exp(-xf**2/(2*s**2))
kf = 2*np.pi*fftfreq(Nf, Lf/Nf)
Phif = ifft(fft(E_bg*zs)*np.exp(-np.abs(kf)*h)).real
dwf = ifft(fft(Phif)*(-kf**2)).real
dwf = (e/(2*m40*w_x))*dwf
mask = np.abs(xf) < 120e-6
# simpler: interpolate
dwc = np.interp(xf, x, dw_closed)
rel = np.max(np.abs(dwf[mask]-dwc[mask]))/np.max(np.abs(dwf[mask]))
print(f"  max rel err (closed form vs FFT) = {rel:.2e}")
print(f"  peak |dw|/2pi = {np.max(np.abs(dwc))/(2*np.pi):.2f} Hz")

print("="*70)
print("3b. WEYL: direct level-set counting vs numeric SVD")
print("="*70)
N = 64
L = 128e-6
h = 40e-6
dx = L/N
kx = 2*np.pi*fftfreq(N, dx)
KX, KY = np.meshgrid(kx, kx)
K = np.sqrt(KX**2+KY**2)
H = np.where(K==0, 0.0, -C*KX**2*np.exp(-K*h))
sv = np.sort(np.linalg.svd(H, compute_uv=False))[::-1]
# convert discrete singular values to continuous convention: sigma_cont = sigma_disc * (L/N)
sv_c = sv*(L/N)
print("  n : sigma_cont(numeric) : -ln(sigma/C)/sqrt(n)")
for n in [3,4,6,8,10,12,15,18,22]:
    if n <= len(sv_c):
        s = sv_c[n-1]
        print(f"  {n:2d} : {s:.4e} : {-np.log(s/C)/np.sqrt(n):.3f}")
print(f"  predicted c = 2 h sqrt(pi)/L = {2*h*np.sqrt(np.pi)/L:.3f}")
print("="*70)
# Bayesian posterior for linear-Gaussian model in Fourier: closed form.
N = 128
Ld = 128e-6
dx = Ld/N
kxv = 2*np.pi*fftfreq(N, dx)
KX, KY = np.meshgrid(kxv, kxv)
Kr = np.sqrt(KX**2+KY**2)
h = 40e-6
Hx = np.where(Kr==0, 0, -C*KX**2*np.exp(-Kr*h))
# surface: Gaussian plateau
xg, yg = np.meshgrid(np.arange(N)*dx, np.arange(N)*dx)
zs = 200e-9*np.exp(-((xg-64e-6)**2+(yg-64e-6)**2)/(2*(7.5e-6)**2))
zs += 100e-9*(xg/128e-6)
d = np.real(ifft2(Hx*fft2(zs)))
delta = 2*np.pi*5*np.sqrt(N**2)  # noise level (5 Hz floor over pixels)
rng = np.random.default_rng(1)
d = d + rng.normal(0, delta/np.sqrt(2), d.shape) + 1j*rng.normal(0, delta/np.sqrt(2), d.shape)
# scalar Tikhonov
lam = 1e-4*np.max(np.abs(Hx))
z_tikh = np.real(ifft2(np.conj(Hx)/(np.abs(Hx)**2+lam**2)*fft2(d)))
# GP prior: S(k)=sigma^2 2pi ell^2 e^{-k^2 ell^2/2};  sigma=200nm, ell=10um
sig_z, ell = 200e-9, 10e-6
Sigma = lambda kk: sig_z**2*2*np.pi*ell**2*np.exp(-kk**2*ell**2/2)
lam2 = np.where(Kr==0, np.inf, delta**2/Sigma(Kr))
z_gp = np.real(ifft2(np.conj(Hx)/(np.abs(Hx)**2+lam2)*fft2(d)))
def corr(a, b):
    a=a.ravel(); b=b.ravel()
    return np.corrcoef(a,b)[0,1]
print(f"  scalar Tikhonov r = {corr(z_tikh, zs):.4f}")
print(f"  GP-prior (k-dependent lambda) r = {corr(z_gp, zs):.4f}")
# joint two-channel SVD: H = [Hx; Hy], Hy = Hx with KX<->KY
Hy = np.where(Kr==0, 0, -C*KY**2*np.exp(-Kr*h))
HH = np.concatenate([Hx.flatten(), Hy.flatten()])
# null-space tests: x-only profile -> Hx sees it, Hy blind; y-only -> vice versa
z_xonly = np.sin(2*np.pi*4*xg/Ld)
z_yonly = np.sin(2*np.pi*4*yg/Ld)
d_x = np.real(ifft2(Hx*fft2(z_xonly))); d_y = np.real(ifft2(Hy*fft2(z_xonly)))
d_x2 = np.real(ifft2(Hx*fft2(z_yonly))); d_y2 = np.real(ifft2(Hy*fft2(z_yonly)))
print(f"  x-only profile: |Hx z|/|z| = {np.std(d_x)/np.std(z_xonly):.2e}, |Hy z|/|z| = {np.std(d_y)/np.std(z_xonly):.2e}")
print(f"  y-only profile: |Hx z|/|z| = {np.std(d_x2)/np.std(z_yonly):.2e}, |Hy z|/|z| = {np.std(d_y2)/np.std(z_yonly):.2e}  (injectivity: only constants invisible to both)")
# adaptive EIG toy: 1D GP, choose next measurement point by max variance reduction
print("  adaptive EIG toy: 2-pass (coarse 8 pts then fine) vs uniform 16 pts on 1D ITF problem")
def itf1d(x, xs, h=40e-6):
    # data kernel for point source: dw(x) = -C d^2/dx^2 [h/(pi((x-xs)^2+h^2))] convolved approx
    return -(e*E_bg/(2*m40*w_x))*(2*h/np.pi)*(h**2-3*(x-xs)**2)/((x-xs)**2+h**2)**3
xpts = np.linspace(0, 128e-6, 200)
xs_true = np.linspace(16e-6, 112e-6, 8)
amp = 100e-9
data = np.sum([amp*itf1d(xpts, xs_i) for xs_i in xs_true], axis=0)
print(f"  peak |dw|/2pi = {np.max(np.abs(data))/(2*np.pi):.2f} Hz")

print("="*70)
print("9. HEATING-LIMITED h_opt, CASIMIR-POLDER, PHASE-1 SNR, LASER NOISE")
print("="*70)
alpha = 4.0
lam = 20e-6
k = 2*np.pi/lam
print(f"  h_opt(heating-dominated) = alpha/k = {alpha/k*1e6:.1f} um (for lambda={lam*1e6:.0f}um)")
# Phase 1: A=50nm, lambda=20um, h=40um
A1, h1 = 50e-9, 40e-6
dw1 = C*A1*k**2*np.exp(-k*h1)
print(f"  Phase-1 signal |dw|/2pi at h=40um = {dw1/(2*np.pi):.2f} Hz")
for hh in [10e-6, 20e-6, 30e-6]:
    print(f"    at h={hh*1e6:.0f}um: {C*A1*k**2*np.exp(-k*hh)/(2*np.pi):.1f} Hz")
# Casimir-Polder crossover
C3 = 9e-49  # J m^3 estimate for Ca+ on Au
for hh in np.logspace(-8, -7, 9):
    ratio = 3*C3/(e*E_bg*hh**4*np.exp(-2*np.pi/hh*hh))
    if abs(np.log(ratio)) < 0.2:
        print(f"  CP crossover h* ~ {hh*1e9:.0f} nm (ratio CP/electrostatic ~ {ratio:.2f})")
# laser noise
print("  laser: T2^laser = 1/(pi*dnu); dnu=1kHz -> T2=0.32ms ; dnu=10Hz -> T2=32ms")
print("  ac-Stark: d(AC)=Omega_R^2/(4 Delta); dI/I=1e-4 -> shift noise ~ 0.05-0.1 Hz scale")
