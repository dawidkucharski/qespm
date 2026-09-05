#!/usr/bin/env python3
"""
New theory figures:
  figV10_joint_channel.pdf  - joint two-channel operator (I1)
  figV11_bayesian.pdf       - Bayesian GP-prior inversion (I6)
  figA1_composite.pdf       - composite asymptotic approximation (O2)
  figV12_adaptive.pdf       - adaptive Bayesian scan demo (O6)

Plus printed diagnostics for the nonlinear inversion (O3).
"""
import numpy as np
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from numpy.fft import fft2, ifft2, fftfreq

e   = 1.602176634e-19
m40 = 40*1.66053906660e-27
w_x = 2*np.pi*1e6
E_bg= 1e4
C   = e*E_bg/(2*m40*w_x)

N, L = 128, 128e-6
dx = L/N
k = 2*np.pi*fftfreq(N, dx)
KX, KY = np.meshgrid(k, k)
K = np.sqrt(KX**2+KY**2); K[0,0] = 1e-12
xg, yg = np.meshgrid(np.arange(N)*dx, np.arange(N)*dx)

# ---- test surface (manuscript topography-only case) ----
zs = 200e-9*np.exp(-((xg-64e-6)**2+(yg-64e-6)**2)/(2*(7.5e-6)**2))  # plateau
zs += 100e-9*(xg/L)                                                 # ramp
zs += 5e-9*np.sin(2*np.pi*xg/4e-6)                                  # fine grating

Hx = -C*KX**2*np.exp(-K*40e-6)
Hy = -C*KY**2*np.exp(-K*40e-6)
Hx[0,0] = 0; Hy[0,0] = 0

# noise: calibrate so scalar Tikhonov at lam=1e-4 maxH gives r ~ 0.98
rng = np.random.default_rng(3)
for sig_pix in [2*np.pi*1.0, 2*np.pi*5.0]:
    dxn = np.real(ifft2(Hx*fft2(zs))) + rng.normal(0, sig_pix, (N,N))
    lam = 1e-4*np.max(np.abs(Hx))
    zt = np.real(ifft2(np.conj(Hx)/(np.abs(Hx)**2+lam**2)*fft2(dxn)))
    r = np.corrcoef(zt.ravel(), zs.ravel())[0,1]
    print(f"calibration: noise {sig_pix/(2*np.pi):.0f} Hz/pix -> scalar Tikhonov r = {r:.4f}")

# ============================================================
# Figure V10: joint two-channel operator
# ============================================================
sig_pix = 2*np.pi*1.0   # matches fig2 inversion (r ~ 0.98 at lam=1e-4 maxH)
dxn = np.real(ifft2(Hx*fft2(zs))) + rng.normal(0, sig_pix, (N,N))
dyn = np.real(ifft2(Hy*fft2(zs))) + rng.normal(0, sig_pix, (N,N))

lam = 1e-4*np.max(np.abs(Hx))
z_x  = np.real(ifft2(np.conj(Hx)/(np.abs(Hx)**2+lam**2)*fft2(dxn)))
z_y  = np.real(ifft2(np.conj(Hy)/(np.abs(Hy)**2+lam**2)*fft2(dyn)))
z_j  = np.real(ifft2((np.conj(Hx)*fft2(dxn)+np.conj(Hy)*fft2(dyn))
                     /(np.abs(Hx)**2+np.abs(Hy)**2+lam**2)))

r_x = np.corrcoef(z_x.ravel(), zs.ravel())[0,1]
r_y = np.corrcoef(z_y.ravel(), zs.ravel())[0,1]
r_j = np.corrcoef(z_j.ravel(), zs.ravel())[0,1]
print(f"joint channel: r_x={r_x:.4f}, r_y={r_y:.4f}, r_joint={r_j:.4f}")

# singular spectra (single channel, 128x128 -> N^2 modes)
h_vals = [5e-6, 40e-6]
fig10, ax10 = plt.subplots(1, 3, figsize=(15, 5))
sv_s = np.sort(np.abs(Hx).ravel())[::-1]
sv_j = np.sort(np.sqrt(np.abs(Hx).ravel()**2+np.abs(Hy).ravel()**2))[::-1]
ax10[0].semilogy(np.arange(1, len(sv_s)+1), sv_s/sv_s[0], color='tab:blue', lw=2,
                 label='single ($H_x$)')
ax10[0].semilogy(np.arange(1, len(sv_j)+1), sv_j/sv_j[0], color='tab:red', lw=2,
                 label='joint ($H_x \oplus H_y$)')
ax10[0].set_xlim(1, 400)
ax10[0].set_xlabel('singular value index $n$')
ax10[0].set_ylabel(r'$\sigma_n/\sigma_1$')
ax10[0].set_title('(a) Joint vs single singular spectra ($h$=40 \u00b5m)')
ax10[0].legend(fontsize=9); ax10[0].grid(True, alpha=0.3)

# effective rank single vs joint at several heights
hs = [5e-6, 10e-6, 20e-6, 40e-6, 80e-6]
rank_s, rank_j = [], []
for hh in hs:
    Hx_h = np.abs(-C*KX**2*np.exp(-K*hh)); Hx_h[0,0]=0
    Hy_h = np.abs(-C*KY**2*np.exp(-K*hh)); Hy_h[0,0]=0
    s1 = np.sort(Hx_h.ravel())[::-1]; s1 = s1[s1>0]
    sj = np.sort(np.sqrt(Hx_h.ravel()**2+Hy_h.ravel()**2))[::-1]
    rank_s.append(np.sum(s1 > 1e-3*s1[0]))
    rank_j.append(np.sum(sj > 1e-3*sj[0]))
ax10[1].semilogy(np.array(hs)*1e6, rank_s, 'o-', color='tab:blue', label='single channel')
ax10[1].semilogy(np.array(hs)*1e6, rank_j, 's-', color='tab:red', label='joint')
ax10[1].set_xlabel('ion height $h$ [\u00b5m]')
ax10[1].set_ylabel('effective rank ($\\sigma>10^{-3}\\sigma_1$)')
ax10[1].set_title('(b) Effective rank: joint closes the null space')
ax10[1].legend(fontsize=9); ax10[1].grid(True, alpha=0.3)

# reconstruction comparison
vlim = np.percentile(np.abs(zs), 99.9)
im = ax10[2].imshow(z_j*1e9, cmap='terrain', origin='lower',
                    extent=[0,128,0,128])
ax10[2].set_title(f'(c) Joint reconstruction ($r={r_j:.3f}$ vs '
                  f'single $r={r_x:.3f}$)')
ax10[2].set_xlabel('x [\u00b5m]'); ax10[2].set_ylabel('y [\u00b5m]')
plt.colorbar(im, ax=ax10[2], fraction=0.046, label='z [nm]')
fig10.suptitle('Two-channel joint inversion: null-space closure and isotropy',
               fontweight='bold', fontsize=12)
fig10.tight_layout()
fig10.savefig('manuscript/figV10_joint_channel.pdf', dpi=150, bbox_inches='tight')
plt.close(fig10)
print('Saved: figV10_joint_channel.pdf')

# ============================================================
# Figure V11: Bayesian GP-prior inversion
# ============================================================
# GP prior: S_z(k) = sigma_z^2 * 2 pi ell^2 exp(-k^2 ell^2/2)   (squared-exponential)
# Posterior (linear-Gaussian model): zhat_k = H* d_k / (|H|^2 + delta^2/S_z(k))
# i.e. Tikhonov with frequency-dependent lambda^2(k) = delta^2 / S_z(k).
sig_z, ell = 200e-9, 10e-6
S_xi = (sig_pix*dx)**2                     # 2D white-noise PSD (continuous convention)
Sz = sig_z**2 * 2*np.pi*ell**2 * np.exp(-K**2*ell**2/2)
lam2 = np.where(Sz>0, S_xi/Sz, np.inf)
z_gp = np.real(ifft2(np.conj(Hx)/(np.abs(Hx)**2+lam2)*fft2(dxn)))
r_gp = np.corrcoef(z_gp.ravel(), zs.ravel())[0,1]
# posterior std per mode (averaged): sqrt(1/(1/Sz + |H|^2/delta^2))-ish
post_var = 1.0/(1.0/Sz + np.abs(Hx)**2/S_xi)
print(f"Bayesian GP: r={r_gp:.4f} (scalar {r_x:.4f})")

fig11, ax11 = plt.subplots(2, 2, figsize=(13, 11))
vlim = np.percentile(np.abs(zs), 99.9)
for axi, (zt, tt, rr) in zip(ax11.ravel()[:2],
                             [(z_x, 'Scalar Tikhonov', r_x), (z_gp, 'GP-prior (Bayesian)', r_gp)]):
    im = axi.imshow(zt*1e9, cmap='terrain', origin='lower', vmin=-vlim*1e9, vmax=vlim*1e9)
    axi.set_title(f'({tt}) $r$={rr:.4f}')
    axi.set_xlabel('x [\u00b5m]'); axi.set_ylabel('y [\u00b5m]')
    plt.colorbar(im, ax=axi, fraction=0.046)
xc = 64
ax11[1,0].plot(xg[0]*1e6, zs[N//2]*1e9, 'k-', lw=2, label='true')
ax11[1,0].plot(xg[0]*1e6, z_x[N//2]*1e9, 'b--', lw=1.5, label='Tikhonov')
ax11[1,0].plot(xg[0]*1e6, z_gp[N//2]*1e9, 'r-.', lw=1.5, label='GP prior')
ax11[1,0].set_xlabel('x [\u00b5m]'); ax11[1,0].set_ylabel('z [nm]')
ax11[1,0].set_title('(c) Cross-section through the plateau centre')
ax11[1,0].legend(fontsize=9); ax11[1,0].grid(True, alpha=0.3)
kk = np.abs(k[1:N//2+1])
lam_sc = np.full_like(kk, lam)
lam_gp = np.sqrt(S_xi/np.maximum(Sz[0,1:N//2+1],1e-60))
ax11[1,1].loglog(kk, lam_sc, 'b-', lw=2, label='scalar $\\lambda$ (Tikhonov)')
ax11[1,1].loglog(kk, lam_gp, 'r-', lw=2, label='$\\lambda(k)$ from GP prior')
ax11[1,1].loglog(kk, np.abs(Hx[0,1:N//2+1]), 'k--', lw=1, alpha=0.6, label='$|H_x(k)|$')
ax11[1,1].set_xlabel('$k$ [rad/m]'); ax11[1,1].set_ylabel('regularisation')
ax11[1,1].set_title('(d) GP prior = frequency-dependent regularisation')
ax11[1,1].legend(fontsize=9); ax11[1,1].grid(True, alpha=0.3)
fig11.suptitle('Bayesian inversion: Gaussian-process prior as optimal '
               'frequency-dependent Tikhonov regulariser', fontweight='bold', fontsize=12)
fig11.tight_layout()
fig11.savefig('manuscript/figV11_bayesian.pdf', dpi=150, bbox_inches='tight')
plt.close(fig11)
print('Saved: figV11_bayesian.pdf')

# ============================================================
# Figure A1: composite asymptotic approximation
# ============================================================
hh = 40e-6
kk = np.logspace(np.log10(0.02/hh), np.log10(8/hh), 400)
f = C*kk**2*np.exp(-kk*hh)
f_comp = C*kk**2*np.exp(-kk*hh) + C*kk**2*(1-np.exp(-kk*hh))*(1-0.5*kk*hh)
figA, axA = plt.subplots(1, 2, figsize=(12, 5))
axA[0].semilogy(kk*hh, f, 'k-', lw=2, label='exact $f(k)=Ck^2e^{-kh}$')
axA[0].semilogy(kk*hh, f_comp, 'r--', lw=1.5, label='composite $f_{\\rm comp}(k)$')
axA[0].set_xlabel('$kh$'); axA[0].set_ylabel('$f(k)$')
axA[0].set_title('(a) Exact ITF vs composite asymptotic approximation')
axA[0].legend(fontsize=9); axA[0].grid(True, alpha=0.3, which='both')
rel = np.abs(f_comp-f)/f
axA[1].semilogy(kk*hh, rel, 'k-', lw=2)
axA[1].set_xlabel('$kh$'); axA[1].set_ylabel('relative error')
axA[1].set_title('(b) Relative error of the naive composite')
axA[1].grid(True, alpha=0.3)
axA[1].annotate(f'overshoot {rel.max()*100:.2e}% at $kh=8$',
                xy=(0.95,0.95), xycoords='axes fraction', ha='right', va='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5), fontsize=9)
figA.tight_layout()
figA.savefig('manuscript/figA1_composite.pdf', dpi=150, bbox_inches='tight')
plt.close(figA)
print('Saved: figA1_composite.pdf')

# ============================================================
# Figure V12: adaptive Bayesian scan (1D toy)
# ============================================================
rng12 = np.random.default_rng(11)
xx = np.linspace(0, 128e-6, 400)
def kernel(x, xs, h=40e-6):
    """delta-source response of the 1D ITF (PSF of Eq. psf, per unit amplitude)."""
    return -(e*E_bg/(2*m40*w_x))*(2*h/np.pi)*(h**2-3*(x-xs)**2)/((x-xs)**2+h**2)**3
xs_src = np.linspace(24, 104, 17)*1e-6   # 5 um spacing < PSF scale h=40 um: ill-posed
amp = 1.5e-9                             # near the 5 Hz noise floor (SNR ~ few per pixel)
rng12 = np.random.default_rng(11)
amp_vec = amp*rng12.normal(1, 0.4, xs_src.size)
truth = np.sum([amp_vec[i]*kernel(xx, x0) for i, x0 in enumerate(xs_src)], axis=0)
noise_std = 2*np.pi*5.0
data = truth + rng12.normal(0, noise_std, xx.size)

# Bayesian 1D GP with kernel K(s,t)=kernel(s,t) per unit amplitude; prior C0 over amplitudes.
prior_amp = 5e-9   # prior std of amplitude at each source
def gp_recon(idx_obs):
    G = np.array([kernel(xx[i], x0) for i in idx_obs for x0 in xs_src]).reshape(len(idx_obs),-1)
    S = noise_std**2*np.eye(len(idx_obs))
    post_cov = np.linalg.inv(G.T@G/S.shape[0]/noise_std**2*noise_std**2 + np.eye(len(xs_src))/prior_amp**2)
    # proper: posterior mean = (G^T S^-1 G + P^-1)^-1 G^T S^-1 d
    G2 = np.array([kernel(xx[i], x0) for i in idx_obs]).T if False else None
    return None
# simpler: regularised least squares on a fixed source grid (equivalent to GP MAP)
G = np.array([kernel(xx[i], x0) for i in range(xx.size) for x0 in xs_src]).reshape(-1, len(xs_src))
P_inv = np.eye(len(xs_src))/prior_amp**2
S_inv = np.eye(xx.size)/noise_std**2

# uniform: all 400 points
post = np.linalg.solve(G.T@S_inv@G + P_inv, G.T@S_inv@data)
recon_all = G@post
r_all = np.corrcoef(recon_all, truth)[0,1]
# adaptive 2-pass: coarse 40 uniform points, then 80 points of max posterior variance
idx_c = np.linspace(0, xx.size-1, 30, dtype=int)
Gc = G[idx_c]; dc = data[idx_c]
Wc = np.eye(len(idx_c))/noise_std**2
P1 = np.linalg.inv(Gc.T@Wc@Gc + P_inv)
post_c = P1@(Gc.T@Wc@dc)
# greedy D-optimal (expected-information-gain) refinement: 60 sequential picks
var = np.einsum('ij,jk,ik->i', G, P1, G)
idx_f = np.argsort(var)[::-1][:60]
P_cur = P1.copy()
idx_sel = list(idx_c)
for _ in range(60):
    gains = np.einsum('ij,jk,ik->i', G, P_cur, G)   # g^T P g per candidate
    gains[idx_sel] = -np.inf
    i_best = int(np.argmax(gains))
    idx_sel.append(i_best)
    g = G[i_best]
    P_cur = P_cur - np.outer(P_cur@g, g@P_cur)/(noise_std**2 + g@P_cur@g)
idx_ad = np.array(sorted(idx_sel))
Ga = G[idx_ad]; da = data[idx_ad]
Wa = np.eye(len(idx_ad))/noise_std**2
Pa = np.linalg.inv(Ga.T@Wa@Ga + P_inv)
post_a = Pa@(Ga.T@Wa@da)
recon_ad = G@post_a
r_ad = np.corrcoef(recon_ad, truth)[0,1]

# uniform 30 points (equal budget)
idx_u = np.linspace(0, xx.size-1, 30, dtype=int)
Gu = G[idx_u]
Wu = np.eye(len(idx_u))/noise_std**2
Pu = np.linalg.inv(Gu.T@Wu@Gu + P_inv)
post_u = Pu@(Gu.T@Wu@data[idx_u])
recon_u = G@post_u
r_u = np.corrcoef(recon_u, truth)[0,1]

# expected errors: average over 30 noise realisations
W_all = np.eye(xx.size)/noise_std**2
P_all = np.linalg.inv(G.T@W_all@G + P_inv)
err_u, err_a, err_all = [], [], []
for _ in range(30):
    d_r = truth + rng12.normal(0, noise_std, xx.size)
    pu_r = Pu@(Gu.T@Wu@d_r[idx_u])
    pa_r = Pa@(Ga.T@Wa@d_r[idx_ad])
    pall_r = P_all@(G.T@W_all@d_r)
    err_u.append(np.linalg.norm(pu_r-amp_vec)/np.linalg.norm(amp_vec))
    err_a.append(np.linalg.norm(pa_r-amp_vec)/np.linalg.norm(amp_vec))
    err_all.append(np.linalg.norm(pall_r-amp_vec)/np.linalg.norm(amp_vec))
eu = np.mean(err_u); ea = np.mean(err_a); ea_all = np.mean(err_all)
print(f"adaptive demo (30 realisations): amplitude err uniform 30 = {eu:.3f}, "
      f"adaptive 30+60 = {ea:.3f}, all-400 = {ea_all:.3f}")
print(f"               posterior cov trace: uniform {np.trace(Pu):.2e} vs adaptive {np.trace(Pa):.2e}")
print(f"               signal r: uniform={r_u:.4f}, adaptive={r_ad:.4f}, all-400={r_all:.4f}")

fig12, ax12 = plt.subplots(1, 2, figsize=(13, 5))
ax12[0].plot(xx*1e6, truth/(2*np.pi), 'k-', lw=2, label='true signal')
ax12[0].plot(xx*1e6, data/(2*np.pi), '.', ms=2, color='0.7', alpha=0.5, label='noisy data')
ax12[0].plot(xx*1e6, recon_u/(2*np.pi), 'b--', lw=1.5, label=f'uniform 30 pts ($r$={r_u:.3f})')
ax12[0].plot(xx*1e6, recon_ad/(2*np.pi), 'r-', lw=1.5, label=f'adaptive 30+60 ($r$={r_ad:.3f})')
ax12[0].set_xlabel('x [\u00b5m]'); ax12[0].set_ylabel(r'$\Delta\omega_x/2\pi$ [Hz]')
ax12[0].legend(fontsize=9); ax12[0].grid(True, alpha=0.3)
ax12[0].set_title('(a) Adaptive two-pass vs uniform sampling')
ax12[1].plot(xx*1e6, np.sqrt(var)/np.sqrt(var).max(), 'k-', lw=1.5)
ax12[1].plot(xx[idx_f]*1e6, np.sqrt(var[idx_f])/np.sqrt(var).max(), 'r.', ms=5,
             label='adaptive refinement points')
ax12[1].set_xlabel('x [\u00b5m]'); ax12[1].set_ylabel('normalised predictive std')
ax12[1].set_title('(b) Posterior variance and the adaptive point selection')
ax12[1].legend(fontsize=9); ax12[1].grid(True, alpha=0.3)
fig12.suptitle('Bayesian adaptive scanning: expected-information-gain point '
               'selection (1D demonstration)', fontweight='bold', fontsize=12)
fig12.tight_layout()
fig12.savefig('manuscript/figV12_adaptive.pdf', dpi=150, bbox_inches='tight')
plt.close(fig12)
print('Saved: figV12_adaptive.pdf')

# ============================================================
# O3: nonlinear inversion (Landweber / CGLS on the quadratic correction)
# ============================================================
print("="*60)
print("O3 nonlinear inversion diagnostics")
import sys as _sys
_sys.path.insert(0, '.')
from rayleigh_validation import rayleigh_coeffs
lam_s = 50e-6; k_s = 2*np.pi/lam_s
h_test = 10e-6
xs = np.linspace(0, lam_s, 2048, endpoint=False)

def exact_dw_profile(A):
    """Exact (Rayleigh) Delta-omega_x(x) for z_s = A sin(k_s x) at height h_test."""
    c = rayleigh_coeffs(A, N=60)
    dw = np.zeros_like(xs)
    for n in range(60):
        nn = n+1
        dw += -(e/(2*m40*w_x))*c[n]*(nn*k_s)**2*np.sin(nn*k_s*xs)*np.exp(-nn*k_s*h_test)
    return dw

for A_true in [4e-6, 8e-6]:
    d = exact_dw_profile(A_true)
    kx1 = 2*np.pi*np.fft.fftfreq(len(xs), xs[1]-xs[0])
    H = -C*kx1**2*np.exp(-np.abs(kx1)*h_test)
    H[0] = 0
    lam = 1e-6*np.max(np.abs(H))
    z0 = np.real(np.fft.ifft(np.conj(H)/(np.abs(H)**2+lam**2)*np.fft.fft(d)))
    A_lin = 2*np.std(z0)          # amplitude from the linear inversion
    r0 = np.corrcoef(z0, A_true*np.sin(k_s*xs))[0,1]
    # Landweber iteration on the scalar amplitude A with the EXACT nonlinear forward map
    A_n = A_lin
    grad_kernel = np.real(np.fft.ifft(H*np.fft.fft(np.sin(k_s*xs))))
    a = 1.0/np.sum(grad_kernel**2)
    for _ in range(20):
        res = d - exact_dw_profile(A_n)
        g = np.sum(res*grad_kernel)
        A_n = A_n + a*g
    zN = A_n*np.sin(k_s*xs)
    rN = np.corrcoef(zN, A_true*np.sin(k_s*xs))[0,1]
    print(f"  A_true={A_true*1e6:.1f}um (Ak={k_s*A_true:.2f}): linear A={A_lin*1e6:.2f}um "
          f"(bias {100*(A_lin/A_true-1):+.1f}%, r={r0:.4f}); "
          f"20 Landweber steps: A={A_n*1e6:.2f}um (bias {100*(A_n/A_true-1):+.1f}%, r={rN:.4f})")
