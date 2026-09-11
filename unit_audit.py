#!/usr/bin/env python3
"""Unit-and-factor audit of the key numbers quoted in the QESPM manuscript.

Recomputes every headline number in main.tex from first principles
(constants: e, m(40Ca), E_bg, omega_x) and compares with the quoted value.
Run before every submission round.
"""
import numpy as np

E = 1.602176634e-19          # elementary charge
M_CA = 40 * 1.66053906660e-27  # 40Ca+ mass [kg]
E_BG = 1e4                   # background field [V/m]
OMEGA_X = 2 * np.pi * 1e6    # secular frequency [rad/s]
C = E * E_BG / (2 * M_CA * OMEGA_X)   # ITF prefactor [m/s]

ok = []


def check(name, computed, quoted, tol=0.05):
    lo, hi = quoted * (1 - tol), quoted * (1 + tol)
    good = lo <= computed <= hi
    ok.append(good)
    print(f"{'PASS' if good else 'FAIL'}  {name:44s} computed={computed:.4g}  quoted={quoted:.4g}")
    return good


print(f"C = e E_bg/(2 m omega_x) = {C:.4g} m/s")
check("C [m/s]", C, 1.92e3, 0.01)

# ITF peak and bandwidth
check("k_opt = 2/h (h=40um) [rad/m]", 2 / 40e-6, 5e4, 0.02)
check("lambda_opt = pi h (h=40um) [um]", np.pi * 40, 125.66, 0.02)
x_lo, x_hi = 0.7613, 4.1559          # roots of x^2 e^-x = 2 e^-2
check("Delta k_FWHM = (x_hi-x_lo)/h (h=40um) [rad/m]", (x_hi - x_lo) / 40e-6, 3.395 / 40e-6, 0.01)
check("Q = 2/(x_hi-x_lo)", 2 / (x_hi - x_lo), 0.589, 0.01)
check("MTF 50% cutoff (root k^2 e^-kh = half max) [x units]", x_hi / 2, 2.0779, 0.01)

# Ramsey SQL / QPN / Type-A
uA = lambda Nrep: 1 / (4 * np.pi * 5e-3 * np.sqrt(Nrep * 0.5))
check("QPN floor, N_rep=1e3, <n>=0.5 [Hz]", uA(1e3), 0.7, 0.05)
check("Type-A u_A(f), N_rep=5.6e3 [Hz]", uA(5.6e3), 0.30, 0.05)

# static-demonstration signal levels (A=50 nm, lambda=20 um grating)
k = 2 * np.pi / 20e-6
A0 = 50e-9
for hq, Hzq in [(40e-6, 5.3), (30e-6, 122), (20e-6, 2.8e3), (10e-6, 65e3)]:
    f = C * k**2 * A0 * np.exp(-k * hq) / (2 * np.pi)
    check(f"signal A=50nm lam=20um h={hq*1e6:.0f}um [Hz]", f, Hzq, 0.06)

# A_min at the ITF peak vs 5 Hz floor (h=40um, k=k_opt)
dw = 2 * np.pi * 5.0
Amin_peak = dw / (C * (2 / 40e-6) ** 2 * np.exp(-2))
check("A_min at ITF peak, 5 Hz floor [nm]", Amin_peak * 1e9, 0.05, 0.3)

# uncertainty budget: charge row (u=1.13e-8 C/m^2, k=2e5, E_bg=1e4)
c_sigma = 1 / (2 * 8.8541878128e-12 * 2e5 * E_BG)
check("charge row u_c [nm]", 1.13e-8 * c_sigma * 1e9, 319, 0.02)

# (2m omega_x/e) prefactor in the sensitivity table
check("2 m omega_x / e [kg rad/(C s)]", 2 * M_CA * OMEGA_X / E, 5.21, 0.02)
check("QESPM curvature sensitivity at QPN [V/m^2/sqrt(Hz)]", 5.21 * 2 * np.pi * 0.7, 23, 0.05)
check("QESPM curvature sensitivity at 5 Hz [V/m^2/sqrt(Hz)]", 5.21 * 2 * np.pi * 5, 160, 0.05)

# dephasing penalty at tau = T2/2
check("sqrt(e) penalty", np.sqrt(np.e), 1.65, 0.02)

# Rayleigh ITF error law 0.14 (kA)^2 at A=2um, lambda=40um
check("ITF rel error at A=2um, lam=40um [%]", 0.1397 * (2 * np.pi * 2e-6 / 40e-6) ** 2 * 100, 1.4, 0.15)

print("\nALL PASS" if all(ok) else "\nSOME CHECKS FAILED")
