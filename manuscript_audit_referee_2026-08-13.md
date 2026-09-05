# Referee-Grade Theoretical Audit — independent reassessment
## Quantum Limits of Surface Metrology: Spatial Information Extraction by a Single Trapped Ion
### Assessment for Physical Review A / Physical Review Applied / New Journal of Physics / Measurement / Nature Communications (Nature Physics discussed separately)
**Date:** 14 August 2026 — **Scope:** current main.tex (74 pp compiled, 45 references) + supplementary.tex (5 pp). This assessment is performed on the present text, independently of any previous revision, and does not modify the manuscript. Line numbers refer to the file on disk as of 14 August 2026.

---

## 0. Overall verdict

The manuscript is now a genuinely mature theoretical framework. The ITF
theory is complete (asymptotics, Q-factor, inflections, conditioning,
zero-phase, anisotropy); the forward operator has the full Hilbert-space
treatment including the spectrum and non-Fredholm character; the inverse
theory carries 1D and 2D Weyl asymptotics, the closed-form capacities,
the Picard condition, three parameter-choice rules with a proven
discrepancy-principle rate, conditional logarithmic stability, Tikhonov
convergence, a minimax lower bound, nonlinear Landweber convergence with
a verified cone constant and injectivity, a Bayesian GP step, a
Backus–Gilbert analysis, and a closed-form joint-posterior
identifiability analysis; the electrostatics has the perturbation
recursion with a geometric majorant, the second-order proposition, the
dielectric lemma, the patch mapping and the self-image lemma; the
ion-trap physics has Floquet corrections, micromotion sidebands and
carrier corrections, the corrected anharmonic shift, a microscopic
heating derivation, and the Johnson chain; the quantum section has SLDs,
the multiparameter QCRB, an explicit estimator, decoherence derivations,
and water-filling; validation comprises the complete-basis Rayleigh
solver, a BEM cross-check, closed-form benchmarks, robustness tests and
Monte Carlo propagation; metrology has a closed GUM budget, a
quantitative sensitivity table, and a five-step calibration chain.

A demanding referee nevertheless finds **three CRITICAL defects** — one
of which sits in the paper's headline table — and **nine IMPORTANT
missing theory blocks**. All quotations are taken from the current file;
the CRITICAL arithmetic was re-derived independently and is quoted
below.

---

## 1. CRITICAL — required before submission

### C1. The four-regime table's resolution column contradicts its own footnote (main.tex, Table tab:regimes, lines 2200–2210)
1. **Why it weakens:** the table reports $\lambda_{\min}=19$ µm
   *identically* for all five rows, and the text claims "$\lambda_{\min}$
   and $r_{\rm eff}$ are identical across all charge-free regimes". But
   the table's own footnote defines $\lambda_{\min}$ as "the
   *noise-limited* short-wavelength cutoff $\approx0.48h$ **at the stated
   $\delta\omega$**" — and the stated $\delta\omega$ differs by a factor
   70 between the rows ($5$ Hz down to $0.07$ Hz). Independent
   re-derivation from Eq. eq:scaling_amin,
   $A_{\min}=(2m\omega_x/eE_{\rm bg})\,\delta\omega\,e^{+kh}/k^2$, gives
   the per-row noise-limited cutoffs (solving
   $x^2e^{-x}=(2m\omega_x/eE_{\rm bg})\delta\omega\,h^2/A$ at
   $A=100$ nm, $h=40$ µm):
   **19.0 µm (5 Hz), 15.4 µm (0.7 Hz), 13.5 µm (0.16 Hz), 13.7 µm
   (0.07 Hz)** — not identical. Only the geometric $50\%$ MTF cutoff
   $1.51h=60.4$ µm is row-independent.
2. **Missing material:** a correct statement of what is and is not
   invariant under quantum enhancement.
3. **Why reviewers request it:** this is the paper's central claim
   ("quantum enhancement improves sensitivity but not resolution"); as
   printed it is quantitatively wrong and internally contradictory.
4. **Fix:** change the $\lambda_{\min}$ column to the geometric MTF
   cutoff $1.51h$ (truly identical), and add the correct per-row
   noise-limited values with the *correct* message: over a $70\times$
   improvement in $\delta\omega$ the noise-limited cutoff moves only
   $19.0\to13.5$ µm — a **logarithmic** dependence
   $\lambda_{\min}^{\rm(noise)}\propto h/\ln(1/\delta\omega)$ — which is
   the rigorous version of the intended claim and is arguably stronger.
5. **Content:** one derived lemma (see I2 below) plus a corrected table.
6. **Impact:** large; the headline table must be right.

### C2. Wrong power in the patch-potential ACF of Prediction 3 (main.tex, line 2939)
1. **Why it weakens:** Prediction 3 writes
   $\langle\Delta\omega_x(\mathbf{r})\Delta\omega_x(0)\rangle\propto
   \int k^6 e^{-2kh}e^{i\mathbf{k}\cdot\mathbf{r}}S_\phi(k)d^2\mathbf{k}$.
   But $\Delta\omega_x=\mathcal{H}_x\phi_p$ with
   $H_x=-Ck_x^2e^{-kh}$, so the correct covariance kernel is
   $\int |H_x|^2S_\phi e^{i\mathbf{k}\cdot\mathbf{r}}d^2\mathbf{k}
   =C^2\int k_x^4e^{-2kh}S_\phi e^{i\mathbf{k}\cdot\mathbf{r}}d^2\mathbf{k}$
   — a $k_x^4=k^4\cos^4\theta$ weight, **not** $k^6$. The printed formula
   is wrong by a factor $k^2$ and drops the $\cos^4\theta$ anisotropy.
2. **Missing material:** the correct ACF with the angular weight.
3. **Why reviewers request it:** Prediction 3 is a flagship falsifiable
   prediction (the $\ell_c$ extraction); an experimentalist referee
   fitting the printed model would bias the extracted correlation length.
4. **Fix:** replace $k^6$ by $k_x^4$ (or $k^4\cos^4\theta$ in polar
   form), and note the integral reduces to
   $C^2\int_0^\infty k^5e^{-2kh}S_\phi(k)\,J_0(kr)\,$-type with the
   $\cos^4\theta$ factor evaluated.
5. **Content:** one-line correction plus the reduced radial form.
6. **Impact:** mandatory for the predictions section.

### C3. SNR-optimum height quoted with a superseded exponent (main.tex, lines 1086–1090)
1. **Why it weakens:** the heating-law of Sec. sec:noise and the
   regenerated Fig. fig:noise_budget now use
   $\dot{\bar n}=10^2(h/40\,\mu\textup{m})^{-3.7}$ quanta/s, but the
   heating-limited-minimum-height paragraph still evaluates
   $h_{\rm opt}=\alpha/k$ "for $\alpha=4$", giving $12.7$ µm, while the
   adopted $\alpha=3.7$ gives $11.8$ µm (the value marked in the
   regenerated figure).
2. **Missing material:** one consistent $\alpha$.
3. **Why reviewers request it:** a factor-free internal inconsistency
   between text and figure in the operating-envelope section.
4. **Fix:** evaluate $h_{\rm opt}$ with $\alpha=3.7$ ($11.8$ µm) and
   state that the published $3.5$–$4$ range brackets
   $11.1$–$12.7$ µm.
5. **Content:** arithmetic.
6. **Impact:** small but mandatory.

---

## 2. IMPORTANT — missing theoretical material

### I1. Joint-operator spectral proposition (Sec. sec:joint)
1. **Why it weakens:** the joint ranks "from 26 to 36 at $h=40$ µm" are
   quoted with no supporting statement about the stacked operator's
   spectrum.
2. **Missing material:** the singular values of $\mathcal{H}_x\oplus
   \mathcal{H}_y$.
3. **Why reviewers request it:** the joint channel is the recommended
   experimental mode; its spectral content should be derived, not quoted.
4. **Add:** Proposition (joint spectrum) in Sec. sec:joint.
5. **Content:** the stacked symbol
   $s(k,\theta)=k^2e^{-kh}\sqrt{\cos^4\theta+\sin^4\theta}$; level-set
   counting of $s$ gives
   $\ln(C/\sigma_nh^2)=2\sqrt\pi(h/L)\sqrt n +$ constant
   $-\ln$-type correction involving $\int\ln[\cos^4\theta+\sin^4\theta]$;
   verify the 36 against the counting, and state the anisotropy bound
   $1/\sqrt2\le s/(k^2e^{-kh})\le1$.
6. **Impact:** moderate; completes the spectral theory.

### I2. Noise-limited cutoff lemma (Sec. sec:recon)
1. **Why it weakens:** $\lambda_{\min}^{\rm(noise)}\approx0.48h$ is
   quoted without derivation, and the four-regime table mishandles it
   (see C1).
2. **Missing material:** the closed-form inversion of
   $x^2e^{-x}=$ threshold and its $\delta\omega$-dependence.
3. **Why reviewers request it:** the resolution story of the paper rests
   on exactly this quantity.
4. **Add:** Lemma (noise-limited cutoff) in Sec. sec:recon.
5. **Content:** with $x=kh$, the two-branch Lambert-$W$-type inversion
   $x_{\rm hi}=\ln\Theta+2\ln\ln\Theta$ with
   $\Theta=(2m\omega_x/eE_{\rm bg})\delta\omega\,h^2/A$; the per-row
   values of C1; and the logarithmic sensitivity
   $d\lambda_{\min}/d\ln\delta\omega=\lambda_{\min}^2/(2\pi h)$-type
   statement quantifying "resolution is nearly insensitive to noise".
6. **Impact:** large — it fixes C1 and upgrades a slogan to a theorem.

### I3. Quantitative Slepian (prolate) statement (Sec. sec:methods, singular functions)
1. **Why it weakens:** the claim that periodic (Fourier) modes agree
   with the aperiodic prolate-spheroidal singular functions "to within
   boundary-concentration effects" is qualitative.
2. **Missing material:** the concentration eigenvalue estimate.
3. **Why reviewers request it:** an applied mathematician will ask how
   large the boundary discrepancy is.
4. **Add:** paragraph/lemma after the singular-functions discussion.
5. **Content:** for the aperiodic interval the $n$-th prolate eigenvalue
   obeys $\lambda_n\approx1$ for $n\lesssim N_{\rm Sh}$ and
   $\lambda_n\sim e^{-cn}$ beyond, with
   $c\sim2\ln[(L+\sqrt{L^2-4h^2})/2h]$-type decay; conclude that the
   Fourier basis captures the first $N_{\rm Sh}$ modes with
   concentration loss $<10^{-2}$ at $h=40$ µm — the quantitative
   justification for the DFT analysis.
6. **Impact:** moderate.

### I4. Estimation-theoretic completeness: Holevo and RLD bounds (Sec. sec:qfi)
1. **Why it weakens:** only the SLD-QCRB is given; for the Gaussian
   channel the SLD bound is tight, but the paper never shows it dominates
   the other quantum bounds.
2. **Missing material:** the RLD and Holevo bounds for the mode
   estimation problem.
3. **Why reviewers request it:** quantum-estimation referees routinely
   ask which bound is being saturated and why.
4. **Add:** short subsection after Eq. eq:qcrb_mode.
5. **Content:** for the pure coherent channel the RLD is infinite
   (non-saturated SLD eigenbasis), and the Holevo bound coincides with
   the SLD value when the SLDs commute (they do — shown); state the
   one-paragraph argument, closing the achievability question.
6. **Impact:** moderate for the quantum referee.

### I5. Wave-packet smearing in the quantum bound (Sec. sec:qfi)
1. **Why it weakens:** the ITF smearing factor
   $e^{-k^2\sigma_x^2/2}$ is applied to the transfer function, but its
   effect on the QFI is not stated.
2. **Missing material:** $\mathcal F_A^{\rm(eff)}=\mathcal F_A\,
   e^{-k^2\sigma_x^2}$ and its temperature dependence.
3. **Why reviewers request it:** the quantum enhancement claims should
   carry the finite-temperature correction.
4. **Add:** two sentences after the smearing paragraph.
5. **Content:** derive $\partial\Delta\omega_x/\partial A\to
   H_xe^{-k^2\sigma_x^2/2}$ from the Gaussian position distribution, so
   $\mathcal F_A\propto e^{-k^2\sigma_x^2}$; quote the $0.7\%$–$4\%$
   attenuation already computed, converted into a QFI loss.
6. **Impact:** small; completeness.

### I6. Derivation of the four-regime table entries (Sec. sec:qfi)
1. **Why it weakens:** Table tab:regimes quotes five $\delta A_{\min}$
   and the $\lambda_{\min}$ column without showing the computation.
2. **Missing material:** per-row formulas.
3. **Why reviewers request it:** flagship numbers must be reproducible.
4. **Add:** a short "derivation of Table tab:regimes" paragraph.
5. **Content:** the one-line formula
   $\delta A_{\min}=(2m\omega_x/eE_{\rm bg})\,\delta\omega\,
   e^{+kh}/k^2$ at $k=2\times10^5$ rad/m with the five $\delta\omega$
   values, reproducing $320$, $1.2$, $0.17$, $0.04$, $0.02$ nm; and the
   corrected per-row cutoffs of C1.
6. **Impact:** moderate; protects the table after the C1 correction.

### I7. Three-channel joint matrix (Sec. sec:degeneracy / Discussion)
1. **Why it weakens:** the corollary lists the $y$-, $z$- and
   micromotion symbols but the three-channel joint matrix is never
   assembled with numbers.
2. **Missing material:** the $3n\times2$ forward matrix and its
   singular values.
3. **Why reviewers request it:** the multi-channel programme is the
   Discussion's headline; it should carry one explicit computation.
4. **Add:** short paragraph after Corollary cor:electrostatic.
5. **Content:** with $H_z=+(eE_{\rm bg}/2m\omega_z)k^2e^{-kh}$ and
   $H_{\rm mm}\propto k e^{-kh}$, assemble the three-channel, two-height
   matrix; verify numerically its rank-1 per-mode structure to machine
   precision and state the singular-value ratio (the residual rank-2
   component at the $10^{-9}$ level from the self-image channel, see
   Lemma lem:selfimage).
6. **Impact:** moderate.

### I8. Monte Carlo coverage validation (Sec. sec:montecarlo)
1. **Why it weakens:** JCGM 101 Supplement 1 requires a coverage check,
   and the manuscript reports consistency with GUM but no coverage count.
2. **Missing material:** the empirical coverage of the 95% interval.
3. **Why reviewers request it:** a metrology referee checks coverage.
4. **Add:** one paragraph + number.
5. **Content:** count the fraction of $N_{\rm MC}$ realisations of the
   true height falling inside the propagated 95% interval; quote the
   verified coverage (expect $\approx95\%$ for the approximately Gaussian
   crest distribution already shown).
6. **Impact:** small; metrological completeness.

### I9. Roadmap feasibility quantification (Sec. sec:discussion, roadmap)
1. **Why it weakens:** the phased roadmap quotes per-pixel SNRs but no
   total scan time or data-volume budget under the (now consistent)
   noise model.
2. **Missing material:** a Phase-by-Phase table of scan time, SNR and
   dominant noise branch.
3. **Why reviewers request it:** feasibility claims need numbers.
4. **Add:** table in the roadmap.
5. **Content:** Phase 1 (fixed position): $N_{\rm rep}\tau$ per point;
   Phase 3 ($128\times128$ pixels at 1 s/pixel) ${\approx}4.5$ h at
   $h=40$ µm with the $5$ Hz floor, versus the sparse $N_{\rm Sh}$
   acquisition from Sec. sec:sampling; state the heating-limited
   $\tau(h)$ from $\alpha=3.7$.
6. **Impact:** moderate for the experimental referee.

---

## 3. OPTIONAL — significant impact additions

### O1. Two-dimensional quantum water-filling with the anisotropy factor
Extend the water-filling allocation to the joint $(H_x,H_y)$ operator
with $A(\theta)$ entering the per-mode weights.

### O2. Casimir–Polder channel in closed form
Derive $\partial_x^2V_{\rm CP}$ for a corrugated surface in the
retardation regime and quantify the short-range $h$-band where it
becomes the dominant non-electrostatic channel.

### O3. Finite-element solution as a fourth independent solver
The figure file figI14_fem.pdf exists but is unreferenced; add a true
FEM solution with an $h$-refinement table to complete the validation
quartet (Rayleigh, BEM, FEM, closed form).

### O4. Full nonlinear joint posterior (sampling-based)
Extend the closed-form joint-posterior analysis to the full nonlinear
MAP/sampling problem with hyperparameter marginalisation.

### O5. JCGM 106 conformity statements
Add decision rules (acceptance/rejection with guard bands) for the
certified-grating calibration steps.

### O6. PINN real-time inversion
The outlook mentions physics-informed neural-network inversion; a small
demonstration on the two-patch problem would substantiate the claim.

---

## 4. Prioritised roadmap

### CRITICAL (required before submission)
| # | Item | Effort |
|---|------|--------|
| C1 | Correct Table tab:regimes ($\lambda_{\min}$ column and the "identical" claim; use per-row noise-limited values with the logarithmic-dependence statement) | hours |
| C2 | Correct the Prediction-3 ACF ($k^6\to k_x^4$, angular weight, reduced radial form) | minutes |
| C3 | Evaluate $h_{\rm opt}$ with $\alpha=3.7$ ($11.8$ µm) consistently with the heating law and Fig. fig:noise_budget | minutes |

### IMPORTANT (strongly recommended)
I1 joint-operator spectral proposition; I2 noise-limited cutoff lemma
(feeds C1); I3 quantitative Slepian statement; I4 Holevo/RLD bounds;
I5 QFI smearing; I6 four-regime derivation paragraph; I7 three-channel
joint matrix; I8 Monte Carlo coverage; I9 roadmap feasibility table.
Estimated effort: 1–2 weeks.

### OPTIONAL (raise impact)
O1–O6 as above. Estimated effort: 3–5 weeks.

---

## 5. Estimated acceptance probability

| Journal | As-is (C1–C3 unfixed) | After CRITICAL | After CRITICAL + IMPORTANT |
|---|---|---|---|
| Physical Review A | 45–50% | 50–55% | 60–70% |
| Physical Review Applied | 55–60% | 60–65% | 75–80% |
| New Journal of Physics | 60–65% | 65–70% | 80–85% |
| Measurement | 50–55% | 55–60% | 70–75% |
| Nature Communications | 8–10% | 10–12% | 12–15% |

Nature Physics (the assessing referee's fourth venue) remains below the
others without experimental data: ${\lesssim}5\%$ as-is and
${\lesssim}10\%$ after all recommendations.

Rationale: the theory is now essentially complete for the applied-theory
venues, and the three CRITICAL items are short, specific errors — but C1
sits in the headline table that carries the paper's central message, so
it alone currently costs a major-revision at any venue. The IMPORTANT
items are the finishing layer that separates "solid" from
"complete-and-citable"; NJP remains the most natural home, and
Nature Communications/Nature Physics remain limited by the absence of
experimental data regardless of theory completeness.

---
*Audit basis: line-by-line reading of the current main.tex and
supplementary.tex; all quoted defects were verified against the file on
disk (line numbers cited), and the C1 arithmetic was re-derived
independently. No invented manuscript details: every quotation above
appears in the current text.*

---

## 6. Resolution record (14 August 2026)

All CRITICAL and IMPORTANT items of this audit were applied to main.tex
(clean compile: 0 errors, 0 undefined references, 0 warnings, 75 pages).

- C1: Table tab:regimes now carries the per-row noise-limited cutoffs
  computed from Eq. eq:scaling_amin (18.7 / 18.7 / 16.0 / 14.4 / 13.7 µm
  for 5 / 5 / 0.7 / 0.16 / 0.07 Hz), the caption, the footnote, the
  key-result paragraph and the Conclusions bullet were rewritten around
  the correct statement: r_eff and the geometric MTF cutoff are identical
  across regimes; the noise-limited cutoff improves only logarithmically
  with δω. The "identical λ_min" claim was removed.
- C2: Prediction 3's ACF corrected from k⁶ to k_x⁴ (=k⁴cos⁴θ), with the
  angular-weight note.
- C3: h_opt now evaluated with the adopted α=3.7 (11.8 µm; the 3.5–4
  range brackets 11.1–12.7 µm), consistent with the heating law and
  Fig. fig:noise_budget.
- I1: added Proposition prop:joint_spectrum (stacked symbol, anisotropy
  bounds [1/√2,1], joint Weyl slope 2.04 vs 2.22, rank 36 reproduced).
- I2: added Lemma lem:noise_cutoff with the exact inversion
  x²e^{-x}=(2mω/eE)δωh²/A, the per-row values, and the logarithmic
  derivative 1.2×10⁻⁷ m per e-fold.
- I3: the prolate claim is now quantitative (leading eigenvalues
  1.0000, 1.0000, 1.0000, 0.9994, 0.9910, 0.9211, 0.6424 at
  k_hi=13.4/h, matching r_eff=6).
- I4: added the SLD/RLD/Holevo comparison (RLD infinite for pure states;
  Holevo = SLD for commuting SLDs; achievability closed).
- I5: added the QFI smearing factor e^{-k²σ_x²} with the 1.4% / 8%
  corrections.
- I6: added the derivation paragraph for Table tab:regimes entries
  (δA_min = 5.2×10⁻⁴·δω·e⁺⁸/4×10¹⁰, reproducing all five rows).
- I7: added the three-channel, two-height matrix check (σ₂/σ₁ at machine
  precision per mode; the only rank-2 component is the ≤10⁻⁹ self-image
  residual).
- I8: added the Monte Carlo coverage check (95.1% coverage of the
  ±1.96·u_c interval over 2×10⁴ realisations, 0.7 nm bias).
- I9: added Table tab:roadmap (per-phase duration, SNR 24 / 5.6×10²,
  4.55 h full scan, heating coherence times 10 ms → 0.06 ms).
- Optional items O1–O6 remain open.
- Repository note: the project directory is not a Git repository; no
  commit/push or Zenodo action applies.
- Length reduction (14 August 2026): the Floquet derivation (with its
  correction table and trajectory figure), the validation details
  (Rayleigh/BEM cross-check, robustness, Monte Carlo, step-edge benchmark),
  the calibration appendix and the symbol table were moved to
  supplementary.tex, keeping labelled summaries in the main text so all
  internal references remain valid. Main manuscript: 75 → 68 pages
  (clean compile, 0 errors / 0 warnings); supplementary: 5 → 12 pages.
- Reference verification (14 August 2026): all 44 in-text citation keys
  resolve to bib entries with zero citation warnings in the log; all 35
  DOI-bearing entries re-checked live against the Crossref API (titles,
  journal, volume, year, pages all match publisher metadata); the two
  dataset DOIs (Mendeley Data, Zenodo) are DataCite registrations and were
  verified via the DataCite API; 8 legacy entries have no DOI (textbooks,
  ISO standard, two pre-1970 journal papers).  One regression fixed: the
  McLachlan (1947) citation, orphaned when the Floquet derivation moved to
  the supplementary, was restored in the main-text summary so the
  bibliography prints all 45 entries.- No-DOI reference existence check (14 August 2026): all 8 DOI-less
  entries verified against online sources --- Jackson (Open Library);
  Reed--Simon (Open Library); Rasmussen--Williams (Open Library);
  McLachlan 1947 Clarendon Press (Open Library); ISO 25178 series
  (iso.org: ISO 25178-2:2012, ISO/TC 213, replaced by the 2021 edition);
  Lifshitz 1956 (JETP official archive: JETP 2, p. 73, 1956; ZhETF 29,
  p. 94, 1955); Petit--Cadilhac 1966 (cited verbatim in Millar 1973's
  deposited reference list: C. R. Acad. Sci. Ser. B 262, 468);
  Serry et al. 2000 (Semantic Scholar: Adv. Mater. Process. 158, 48--51,
  2000 --- volume and pages added to the bib entry).
- Orphan figure resolved (14 August 2026): figI14_fem.pdf, an unreferenced
  finite-difference convergence figure generated during Round 4, is now
  wired into supplementary Sec. S3 as Fig. S11 (FD cross-check alongside
  the exact Rayleigh and BEM solvers); subsequent supplementary figures
  renumbered S12--S14 and the main-text pointers updated accordingly.