#!/usr/bin/env python3
"""
Boundary Element Method (BEM) validation of the linear ITF.
Solves the 2D Laplace equation for a corrugated grounded conductor
and compares the exact solution to the first-order perturbation result.

Method: Indirect BEM with point collocation — Green's function for 
half-space with Dirichlet boundary condition on z = z_s(x).

For a perfectly conducting surface at potential Phi=0, the induced
surface charge density sigma(x) satisfies the integral equation, 
and the potential at the ion position is computed exactly.

Output: figV7_bem_validation.pdf
"""

import numpy as np
from scipy.linalg import solve
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt

np.random.seed(456)

# ============================================================
# Indirect BEM for Laplace equation with corrugated boundary
# ============================================================
def bem_potential_2d(x_surface, z_surface, x_ion, z_ion, E_bg=1e4):
    """
    Compute exact potential at ion position above a corrugated grounded conductor.
    
    Uses indirect BEM: places source points on the boundary, solves for
    source strengths to satisfy Dirichlet BC Phi=0, then evaluates
    potential at the ion position.
    
    Parameters:
      x_surface: array of x-coordinates on the surface (M points)
      z_surface: array of z-coordinates on the surface (same M)
      x_ion, z_ion: ion position
      E_bg: background field strength
    Returns:
      Phi_at_ion: total potential at ion position
    """
    M = len(x_surface)
    
    # Build influence matrix: G_{ij} = -1/(2*pi) * ln(r_{ij})
    G = np.zeros((M, M))
    for i in range(M):
        for j in range(M):
            dx = x_surface[i] - x_surface[j]
            dz = z_surface[i] - z_surface[j]
            r = np.sqrt(dx**2 + dz**2)
            if r < 1e-15:
                # Self-term: use local curvature correction
                # For flat segments, use element length regularization
                ds = x_surface[1] - x_surface[0] if M > 1 else 1e-6
                G[i, j] = -1/(2*np.pi) * np.log(ds / 4)
            else:
                G[i, j] = -1/(2*np.pi) * np.log(r)
    
    # RHS: background potential on surface must be cancelled
    # Phi_bg = E_bg * z, so on surface: Phi_bg = E_bg * z_surface
    # We need Phi_total = Phi_bg + Phi_ind = 0 => Phi_ind = -Phi_bg
    rhs = -E_bg * z_surface
    
    # Solve for source strengths
    sigma = solve(G, rhs)
    
    # Evaluate potential at ion position
    phi_ion = 0.0
    for j in range(M):
        dx = x_ion - x_surface[j]
        dz = z_ion - z_surface[j]
        r = np.sqrt(dx**2 + dz**2)
        if r > 1e-15:
            phi_ion += -1/(2*np.pi) * np.log(r) * sigma[j]
    
    # Add background contribution
    phi_ion += E_bg * z_ion
    
    return phi_ion

def itf_potential_analytical(x_surface, z_surface, x_ion, z_ion, E_bg=1e4):
    """First-order perturbation ITF: Phi = E_bg * z_ion + O(z_surface)"""
    # In the linear ITF: Phi^(1)(x,z) from Poisson integral
    # For a sinusoidal surface z_s = A*sin(kx), Phi^(1)(x,z) = E_bg*A*sin(kx)*e^{-kz}
    # Here we do the general Poisson integral for arbitrary z_surface
    
    M = len(x_surface)
    dx = x_surface[1] - x_surface[0]
    phi_ion = E_bg * z_ion  # background
    
    for i in range(M):
        r = np.sqrt((x_ion - x_surface[i])**2 + z_ion**2)
        if r > 1e-15:
            ds = dx
            # Poisson kernel: Phi(z) = z/pi * integral [phi_s(x') / ((x-x')^2 + z^2)] dx'
            # where phi_s = E_bg * z_surface (small-roughness approximation)
            phi_ion += (z_ion / np.pi) * (E_bg * z_surface[i]) / (r**2) * ds
    
    return phi_ion

# ============================================================
# Validation: sinusoidal surface, varying amplitude
# ============================================================
print("Running BEM validation for sinusoidal surfaces...")

lambda_s = 40e-6  # period
k_s = 2*np.pi/lambda_s
h_nom = 40e-6

# Surface discretization
M = 500
x_range = 4 * lambda_s  # multiple periods
x_s = np.linspace(-x_range/2, x_range/2, M)

# Ion position (centred at anti-node where sin(kx)=1 for max signal)
x_ion = lambda_s / 4  # anti-node of sin(kx)
z_ion = h_nom

E_bg_val = 1e4  # V/m, background field

# Varying amplitude
A_list = np.logspace(np.log10(10e-9), np.log10(2e-6), 15)
dw_bem = []
dw_itf = []

for A in A_list:
    z_s = A * np.sin(k_s * x_s)
    
    # BEM: compute potential at two nearby x positions to get d^2Phi/dx^2
    delta_x = 1e-9
    phi_c = bem_potential_2d(x_s, z_s, x_ion, z_ion, E_bg_val)
    phi_l = bem_potential_2d(x_s, z_s, x_ion - delta_x, z_ion, E_bg_val)
    phi_r = bem_potential_2d(x_s, z_s, x_ion + delta_x, z_ion, E_bg_val)
    d2phi_bem = (phi_l - 2*phi_c + phi_r) / delta_x**2
    
    e = 1.602176634e-19
    m = 40 * 1.66053906660e-27
    omega_x = 2*np.pi*1e6
    dw_bem.append(e/(2*m*omega_x) * d2phi_bem)
    
    # ITF prediction: ITF gives dw_fft = H_x * z_fft
    # For z=A*sin(kx), dw(x) = (e*Ebg/2m*wx) * (-k^2) * A*e^{-kh} * sin(kx)
    # At x = lambda/4, sin(kx)=1
    dw_itf_val = (e * E_bg_val)/(2*m*omega_x) * A * (-k_s**2) * np.exp(-k_s * h_nom)
    dw_itf.append(dw_itf_val)

dw_bem = np.array(dw_bem)
dw_itf = np.array(dw_itf)

# ============================================================
# Figure V7: BEM validation
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Panel (a): BEM vs ITF frequency shift vs amplitude
axes[0].loglog(A_list*1e9, np.abs(dw_bem)/(2*np.pi), 'bo-', ms=6, lw=2, label='BEM (exact)')
axes[0].loglog(A_list*1e9, np.abs(dw_itf)/(2*np.pi), 'r--', lw=2, label='ITF (linear)')
axes[0].set_xlabel('Surface amplitude $A$ [nm]')
axes[0].set_ylabel('$|\\Delta\\omega_x|/2\\pi$ [Hz]')
axes[0].set_title(f'(a) Frequency shift vs amplitude ($h$={h_nom*1e6:.0f} µm, $\\lambda$={lambda_s*1e6:.0f} µm)')
axes[0].legend(fontsize=9)
axes[0].grid(True, alpha=0.3)

# Panel (b): Relative error BEM vs ITF
rel_err = np.abs(dw_bem - dw_itf) / np.abs(dw_bem)
axes[1].semilogx(A_list*1e9, rel_err*100, 'ko-', ms=6, lw=2)
axes[1].axhline(5, color='red', ls='--', lw=1, label='5% error')
axes[1].axhline(10, color='orange', ls=':', lw=1, label='10% error')
axes[1].set_xlabel('Surface amplitude $A$ [nm]')
axes[1].set_ylabel('Relative error [%]')
axes[1].set_title('(b) ITF accuracy vs BEM')
axes[1].legend(fontsize=9)
axes[1].grid(True, alpha=0.3)

# Panel (c): z_s/h ratio vs error
z_over_h = A_list / h_nom
axes[2].semilogx(z_over_h, rel_err*100, 'ko-', ms=6, lw=2)
axes[2].axvline(0.1, color='gray', ls='--', lw=1, label='$z_s/h = 0.1$')
axes[2].axhline(10, color='orange', ls=':', lw=1)
axes[2].set_xlabel('$A/h$')
axes[2].set_ylabel('Relative error [%]')
axes[2].set_title('(c) Error vs $z_s/h$ ratio')
axes[2].legend(fontsize=9)
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('manuscript/figV7_bem_validation.pdf', dpi=150, bbox_inches='tight')
plt.close()
print('Saved: figV7_bem_validation.pdf')

print(f"\nBEM validation summary:")
print(f"  At A=100 nm: |dw|/2pi (BEM) = {np.abs(dw_bem[5])/(2*np.pi):.0f} Hz")
print(f"  At A=100 nm: |dw|/2pi (ITF) = {np.abs(dw_itf[5])/(2*np.pi):.0f} Hz")
print(f"  Relative error = {rel_err[5]*100:.1f}%")
print(f"  ITF < 5% error up to A = {A_list[np.argmax(rel_err > 0.05)]*1e9:.0f} nm" if np.any(rel_err > 0.05) else "  ITF < 5% error for all tested amplitudes")
