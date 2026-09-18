# STMP Restructure — Move Map and Completeness Check

Date: 2026-09-18. Main manuscript reduced from ~72 pages to **32 pages**
(0 errors / 0 undefined / 0 `??` / 0 warnings); Supplementary Material is
**30 pages** (same clean status). No scientific content was deleted: all
material removed from the main text was relocated to the Supplementary
Material (Secs. S1–S16, grouped thematically) or compressed with the
central equations retained.

## A. New main-paper structure

1. Introduction
2. Measurement concept and metrological measurand
3. QESPM forward model and spatial transfer function
4. Spatial information limits and identifiability
5. Quantum precision and metrological uncertainty
6. Numerical validation and realistic surface reconstruction
7. Operating envelope, limitations and experimental implementation
8. Discussion
9. Conclusions

## B. Supplementary organisation

- S1 Supplementary figure index and AFM evaluation details
- S2 Floquet analysis
- S3 Exact (Rayleigh–Floquet) cross-validation (incl. BEM/FD hierarchy figure)
- S4 Robustness to model violation (incl. dielectric factor, validity diagram)
- S5 Monte Carlo uncertainty propagation (incl. Type-B/Type-A derivations, correlations, falsifiable predictions, sensitivity-comparison Table S2)
- S6 Analytical benchmark: step-edge response
- S7 Calibration and traceability considerations (incl. traceability-chain figure)
- S8 Symbols and constants
- S9 Experimental roadmap
- S10 Forward-model derivations (incl. higher-order corrections, convergence majorant, equilibrium shift, self-image, asymptotics, sensitivity coefficients, minimum detectable amplitude, Figs. S17–S19)
- S11 Quantum resource allocation (incl. multiparameter QCRB and optimal measurement)
- S12 Stability and operating regime
- S13 Response to broadband surface roughness
- S14 Sequential and adaptive measurement strategies (incl. sampling theorem)
- S15 Johnson–Nyquist noise derivations
- S16 Regularisation-theory proofs (incl. operator-theoretic details, Weyl asymptotics, joint operator Fig. S23, Picard/L-curve/GCV Fig. S24, Backus–Gilbert Fig. S25, Bayesian GP and joint posterior Figs. S26–S28, nonlinear inversion)

## C. Moved-to-Supplement map (section by section)

| Moved content (old main location) | Now in | Main-text cross-reference |
|---|---|---|
| Second-order perturbation propositions + proofs | S10 | Sec. 3.1 (summary retained) |
| Perturbation-series convergence lemma + geometric majorant | S10 | Sec. 3.1 |
| Equilibrium-position shift analysis | S10 | Sec. 3.1 |
| Self-image lemma, η(k,h) values, two-term ITF | S10 | Sec. 3.1 |
| Zero-phase filter proof | S10 | Sec. 3.2 |
| Low-k/high-k expansions, composite bounds, Q-factor derivation, inflection lemma | S10 (Fig. S17) | Sec. 3.2 (Q, cutoffs retained) |
| Sensitivity-coefficient/conditioning lemma | S10 | Sec. 5.4 (budget link) |
| Extended L² operator-theory analysis, prolate eigenvalues, singular gallery | S16 (Fig. S18) | Sec. 4.2 |
| Floquet details | S2 | Sec. 3.1 |
| Numerical implementation chain | main Sec. 3.4 (compressed; code in Zenodo) | Data Availability |
| Microscopic noise models (heating, laser, Johnson) | S2, S15 | Sec. 3.6 |
| Heating-channel spatial kernel derivation | S15 | Sec. 4.4 |
| Dielectric-film lemma and lossy-film analysis | S4 | Sec. 2.3 (F(k) retained) |
| Validity diagram | S4 (Fig. S13) | Sec. 2.3 |
| SVD details, Weyl lemmas and proofs, information capacity | S16 | Sec. 4.2 (law + table retained) |
| Joint-operator spectrum/corollary proofs + figure | S16 (Fig. S23) | Sec. 4.3 |
| Theorem proof details (rank-1) | S16 | Sec. 4.4 (proof sketch retained) |
| Picard condition, L-curve/GCV/DP selection + figure | S16 (Fig. S24) | Sec. 5.3 |
| Source set, conditional-stability and minimax theorem proofs | S16 | Sec. 5.3 (rates retained) |
| Backus–Gilbert analysis + figure | S16 (Fig. S25) | Sec. 5.3 (result retained) |
| Bayesian GP posterior, joint posterior + figures | S16 (Figs. S26–S28) | Sec. 5.3 |
| Nonlinear (Landweber) inversion | S16 | Sec. 5.3 |
| Reconstruction conditions, noise-cutoff lemma, A_min band + figure | S10 (Fig. S19) | Sec. 4.1 |
| Sampling theorem and throughput | S14 | Sec. 7.2–7.3 |
| Type-B justifications, Type-A derivation, correlations | S5 | Sec. 5.4 |
| Falsifiable experimental predictions | S5 | Sec. 6 (summary) |
| Sensitivity comparison table | S5 (Table S2) | Sec. 8.1 |
| Traceability-chain figure | S7 (Fig. S16) | Sec. 5.4 |
| Solver-hierarchy figure | S3 (Fig. S11) | Sec. 6.1 |
| Extended Monte Carlo details | S5 | Sec. 6.5 |
| Step-edge derivation, Gaussian-bump benchmark | S6 | Sec. 6.5 |
| Multiparameter QCRB, optimal measurement, resource allocation | S11 | Sec. 5.1 |
| Model-violation details | S4 | Sec. 6.5 |
| Roadmap details | S9 | Sec. 7.3 |
| Adaptive-scanning details | S14 | Sec. 7.3 |
| Figures moved from main to supplement | Figs. S11, S13, S16–S19, S23–S28 | as above |

Main-paper figures kept (11): geometry (Fig. 1), ITF (Fig. 2), forward chain
(Fig. 3), noise budget (Fig. 4), SVD (Fig. 5), QFI (Fig. 6), uncertainty
decomposition (Fig. 7), inversion (Fig. 8), AFM evaluation (Fig. 9),
published comparison (Fig. 10), sensitivity–standoff landscape (Fig. 11).
Main-paper tables kept (4): quantity hierarchy, four-regime comparison,
uncertainty budget, technique comparison.

## D. Main-paper completeness check

Every item on the reviewer reading list is present in the main manuscript:
1. What QESPM is and the problem it addresses — Sec. 1.
2. Measurand/indication/influence quantities — Sec. 2.1, Table 1.
3. Physical principle — Secs. 2.2–3.1 (boundary perturbation, Poisson integral, Eq. (1)–(4)).
4. Forward operator and ITF — Eq. (5)–(6), Sec. 3.2.
5. Information loss with standoff — Sec. 4.1–4.2 (effective-rank table, Weyl law, Fig. 5).
6. Topography–charge identifiability — Theorem 1 + Corollaries 1–2, Sec. 4.4.
7. Quantum precision limit and what enhancement can/cannot do — Eq. (11), Table 2, Sec. 5.1–5.2.
8. Reconstruction/inversion strategy — Eq. (12)–(14), Sec. 5.3.
9. Metrological uncertainty framework — Table 3, Sec. 5.4, traceability summary.
10. Principal numerical validation — Sec. 6 (Rayleigh 0.01%→0.14(kA)², BEM ≤0.6%, FD convergence, synthetic r=0.981).
11. Realistic uncertainty result — 319 nm, explicitly conditional (Sec. 5.4).
12. Computational vs experimental status — stated in abstract, Secs. 6.3, 7, Data Availability.
13. Limitations and implementation — Sec. 7.
14. Metrological significance and conclusions — Secs. 8–9.

All central equations, all four theorems/lemmas with their statements, all
principal numerical values and every major conclusion appear in the main
paper; no major conclusion depends on supplement-only material.

## E. Final consistency audit (automated)

- `main.tex`: 11 figures, 4 tables; all figure files exist; every float label referenced; no broken `\ref`; compile 0/0/0/0; references begin on page 30 of 32 (no floats after bibliography).
- `supplementary.tex`: 28 figures, 4 tables; all figure files exist; every float label referenced; no broken `\ref`; compile 0/0/0/0; own bibliography added (citations resolve).
- Every hard-coded "Supplementary Sec./Fig./Table Sx" reference in `main.tex` was aligned with the actual numbering extracted from the compiled `.aux`.
- `unit_audit.py`: 19/19 PASS after the restructure (all quoted numbers unchanged).
- Resolution terminology audit retained: ESF FWHM 0.34h, PSF first zero h/√3, Rayleigh-style 1.15h, MTF-50 1.51h, λ_opt=πh used as distinct metrics; the L² operator claim states the Fourier-multiplier/non-compact structure correctly.
