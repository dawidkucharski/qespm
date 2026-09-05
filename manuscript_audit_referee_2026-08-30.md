# Referee Audit — Round 6 (from scratch, independent of previous rounds)

**Date:** 30 August 2026
**Auditor stance:** senior theoretical physicist / applied mathematician; refereeing for Physical Review A, Physical Review Applied, New Journal of Physics, Measurement, Nature Physics.
**Scope:** identify everything *missing* from the theoretical development. The scientific idea is assumed correct and unchanged. This audit does **not** edit the manuscript.
**Basis:** complete reading of the current `main.tex` (70 pp, 45 refs) and `supplementary.tex` (13 pp, S1–S8) as they exist on disk today, plus independent recomputation of every quantitative claim quoted below. Quotations are verbatim; my own arithmetic is labelled "referee's check". All findings below were located by line in the current files (line numbers cited).

---

## Part A — Findings

### CRITICAL (required before submission)

#### C1. The sensitivity-comparison table's QESPM entries are irreproducible from their stated definition

**Location:** Table `tab:sensitivity`, row "QESPM (this work)", lines 3217–3218. The caption states the QESPM entry is the "gradient sensitivity $(2m\omega_x/e)\,\delta\omega$ from this work"; the row quotes **72 V/m²/√Hz (QPN)** and **5.1×10² V/m²/√Hz (5 Hz floor)**.

**Referee's check.** $(2m\omega_x/e)\cdot\delta\omega$ with $m=40$ u, $\omega_x/2\pi=1$ MHz, $\delta\omega=2\pi\times0.7$ rad/s gives $5.21\times10^{-9}\times4.40=2.3\times10^{-8}$ V/m²/√Hz — nine orders of magnitude below the quoted 72. (With the 5 Hz floor: $1.6\times10^{-7}$ vs $5.1\times10^2$ — the same factor $3.1\times10^9$.) The quoted values are instead consistent (to $\sim$5%) with $k^2E_{\rm bg}\,\delta A_{\min}$ at the evaluation point $k=2\times10^5$ rad/m — i.e. the noise-equivalent **surface-potential curvature** ($68$ and $488$ V/m²), not the caption's definition.
1. **Why it weakens the paper.** A flagship comparison table against KPFM/EFM/NV/electrometry contains two numbers that cannot be produced from the paper's own equations as captioned. A metrology referee will recompute the definition and flag the table as unreliable.
2. **What is missing.** Either the correct definition (noise-equivalent curvature $k^2E_{\rm bg}\delta A_{\min}$, with $k$ and $h$ stated), or the values recomputed from the caption's $(2m\omega_x/e)\delta\omega$.
3. **Why reviewers request it.** The table is the paper's quantitative positioning against competitor techniques; a 10⁹ discrepancy in the headline entries undermines it.
4. **Proposed subsection.** Fix the caption/values in Table `tab:sensitivity`; add one sentence deriving the QESPM entry from Eq. `eq:scaling_amin`.
5. **Content.** One-line derivation: $\delta(\partial_x^2\Phi)=k^2E_{\rm bg}\delta A_{\min}=k^2E_{\rm bg}(2m\omega_x/eE_{\rm bg})\delta\omega\,e^{+kh}/k^2 = (2m\omega_x/e)\,\delta\omega\,e^{+kh}$, evaluated at $kh=8$.
6. **Impact.** +6–8% across all target journals.

---

#### C2. The comparison table still advertises "${\sim}10^{-2}$ nm vertical resolution" without the sensitivity/uncertainty qualifier

**Location:** Table `tab:comparison`, QESPM row, line 3188: "QESPM & ${\sim}1.5h$ (MTF) & ${\sim}10^{-2}$\,nm & No & UHV native".

1. **Why it weakens the paper.** The manuscript now correctly states (Abstract, §3.1, §3.14) that quantum-limited *detection sensitivity* differs from *reconstruction uncertainty*, and that the realistic combined standard uncertainty is 319 nm. A table cell reading "Vertical res. ${\sim}10^{-2}$ nm" with no qualifier is precisely the "picometre surface microscope" claim the paper elsewhere disavows; a referee comparing Table `tab:comparison` with Table `tab:uncertainty` sees a contradiction.
2. **What is missing.** The qualifier: the $10^{-2}$ nm figure is the QPN detection sensitivity at $k_{\rm opt}=2/h$ in the charge-free regime; the realistic reconstruction uncertainty is 319 nm.
3. **Why reviewers request it.** Internal consistency between the positioning table and the uncertainty section.
4. **Proposed subsection.** Amend the cell to "$\sim10^{-2}$ nm (QPN sensitivity, charge-free; realistic $u_c=319$ nm)" or add a table footnote.
5. **Content.** One line + footnote.
6. **Impact.** +4–6%.

---

#### C3. False injectivity statement and garbled null-space characterisation (§2.1, closing paragraph, lines 765–769)

**Quote.** "The null space of $\mathcal{H}_h$ (for the $x$-channel alone) consists of all surface profiles with $\tilde{z}_s(k_x=0, k_y)=0$, i.e., profiles that depend only on $y$. … $\mathcal{H}_h$ is injective on the subspace of zero-mean functions; the non-injectivity is limited to the physically trivial constant mode."

1. **Why it weakens the paper.** Two errors in one paragraph of the operator-theory section: (i) the null space is profiles whose transform is **supported on** $\{k_x=0\}$ (zero off the line), not profiles with $\tilde z_s(k_x=0,k_y)=0$ — the clause that follows ("profiles that depend only on $y$") is correct, the formula is not; (ii) the single-channel operator is **not** injective on zero-mean functions — any zero-mean y-only profile is annihilated. Injectivity holds only for the joint operator $\mathcal{H}_x\oplus\mathcal{H}_y$ (proved correctly in Proposition `prop:joint`).
2. **What is missing.** A correct Lemma: $\ker\mathcal{H}_h=\{\tilde z_s:\operatorname{supp}\tilde z_s\subseteq\{k_x=0\}\}$; and the injectivity statement restricted to the joint operator.
3. **Why reviewers request it.** A self-contradiction in the most theorem-dense section.
4. **Proposed subsection.** Replace the closing paragraph of §2.1 with `Lemma [Null space and injectivity]`.
5. **Content.** Statement + one-line proof from the symbol.
6. **Impact.** +8–10% at PRA/PRApplied.

---

#### C4. False "exactly once" claim in the proof of Proposition `prop:spectrum` (line 723)

**Quote.** "The symbol $H_x=-Ck_x^2e^{-kh}$ is real, continuous, and takes every value in $[-M,0]$ exactly once along the line $k_y=0$ ($k_x\in[0,\infty)$)…"

1. **Why it weakens the paper.** $k_x\mapsto Ck_x^2e^{-k_xh}$ on $[0,\infty)$ rises from 0 to $M$ at $k_x=2/h$ and decays back to 0: every interior value is attained **twice**, not once. The conclusion (essential range $=[-M,0]$) is correct; the asserted reason is false.
2. **What is missing.** The correct continuity argument; optionally the positive remark that the two preimages are exactly the roots $x_{\rm lo},x_{\rm hi}$ of Corollary `cor:qfactor`.
3. **Why reviewers request it.** A false sentence inside a proof.
4. **Proposed subsection.** Correct the one sentence in the proof.
5. **Content.** One-sentence fix.
6. **Impact.** +4–6% at PRA/PRApplied.

---

#### C5. Equilibrium-position-shift numbers are inconsistent with the manuscript's own parameters by 3–4 orders of magnitude (line 325)

**Quote.** "For the nominal parameters ($h = 40$\,\textmu m, $|\Phi^{(1)}| \sim 100$\,\textmu V, feature scale $\sim 30$\,\textmu m), I find $|\Delta\mathbf{r}| \sim 10^{-10}$\,m and $|\Delta\mathbf{r}|/\lambda_{\rm feature} \sim 3 \times 10^{-6}$. The resulting correction to $\Delta\omega_x$ … is $\mathcal{O}(|\Delta\mathbf{r}|^2/\lambda_{\rm feature}^2) \sim 10^{-11}$ relative."

**Referee's check.** From the paper's own formula $\Delta\mathbf{r}=-(e/m\omega_x^2)\nabla\Phi^{(1)}$: $e/m\omega_x^2=6.1\times10^{-8}$ m per (V/m); with $\nabla\Phi^{(1)}\sim2\pi|\Phi^{(1)}|/\lambda_{\rm feature}\approx21$ V/m one obtains $|\Delta\mathbf{r}|\approx1.3\times10^{-6}$ m (four orders above the quoted value; even the cruder $\Phi/\lambda$ gradient gives $2.0\times10^{-7}$ m). Consequently $(\Delta r/\lambda_{\rm feature})^2\sim10^{-3}$, not $10^{-11}$ — the stated correction is eight orders too small. The qualitative conclusion (correction small, $\sim0.1\%$) survives; the quoted numbers do not.
1–4. As a quantitative Methods claim, this is the first arithmetic a referee repeats. **Proposed subsection:** rewrite the paragraph with the corrected displacement ($\sim10^{-6}$ m), corrected relative correction ($\sim10^{-3}$), and the honest statement that the residual is 0.1%-level and absorbed into the Type-B budget.
5. **Content.** Two-line computation.
6. **Impact.** +5–8% across journals.

---

#### C6. The printed Ramsey frequency-uncertainty formula is missing a factor $2\pi$ (four locations)

**Locations:** line 846 (QPN floor definition), line 942 (heating-limited minimum-height paragraph), line 2830 (Type-A derivation), line 2936 (Prediction 5).

**Quote (line 2830).** "$u_A(f) = 1/(2\tau\sqrt{N_{\rm rep}\langle n\rangle})$; with $\tau = T_2/2 = 5$\,ms, $\langle n\rangle=0.5$ and $N_{\rm rep}=5.6\times10^3$…, $u_A(f) = 0.30$\,Hz".

**Referee's check.** The printed formula gives $1/(2\cdot5\,\text{ms}\cdot\sqrt{2800})=1.89$ Hz, not 0.30 Hz; the standard Ramsey SQL formula $\delta f=1/(4\pi\tau\sqrt{N_{\rm rep}\langle n\rangle})$ reproduces 0.30 Hz exactly (and the 0.7 Hz QPN floor at $N_{\rm rep}=10^3$). The numerical values throughout (including the laser-noise entries 11 Hz / 0.1 Hz in the noise budget) were computed with the $4\pi\tau$ form; the printed formula has lost the $2\pi$ in every occurrence.
1. **Why it weakens the paper.** The QPN floor drives the four-regime table, roadmap SNRs, and the quantum-enhancement claims; substituting the printed formula gives 4.5 Hz and a 6.4× inconsistency.
2. **What is missing.** The correct formula with its one-line derivation ($\delta\phi_{\rm SQL}=1/(2\sqrt{\langle n\rangle})$ per shot; $\delta\omega=\delta\phi/(2\pi\tau\sqrt{N_{\rm rep}})$).
3. **Why reviewers request it.** Formula–value mismatch in the flagship table.
4. **Proposed subsection.** Correct all four occurrences; state the formula once.
5. **Content.** One equation + one sentence.
6. **Impact.** +8–10% at PRA/PRApplied/Measurement.

---

### IMPORTANT (strongly recommended)

#### I1. The second-order frequency-shift kernel is never derived; only its numerical coefficient is quoted

$\Phi^{(2)}$ is derived (Eq. `eq:phi2`), but $\Delta\omega_x^{(2)}=(e/2m\omega_x)\,\partial_x^2\Phi^{(2)}|_{(x,y,h)}$ is never evaluated analytically; the coefficient $0.14\,(kA)^2$ (Lemma `lem:perturbation`) appears as a numerical accident.
1. **Why it weakens the paper.** The nonlinear electrostatic correction is not carried through to the observable.
2. **What is missing.** The Fourier form $\Delta\tilde\omega_x^{(2)}(k)=-\frac{eE_{\rm bg}}{2m\omega_x}e^{-kh}k_x^2\int\tilde z_s(k')|\mathbf k-\mathbf k'|\tilde z_s(\mathbf k-\mathbf k')\frac{d^2k'}{(2\pi)^2}$; its reduction for a sinusoid (explicit $2k$-harmonic coefficient); proof that the Rayleigh $c_2$ coefficient equals this closed form.
3. **Why reviewers request it.** Electrostatics referees expect the analytic structure (sign, $h$-dependence) of the nonlinear term.
4. **Proposed subsection.** `Proposition [Second-order frequency shift]` after `prop:secondorder`.
5. **Content.** Two displayed equations + cross-check against the existing exact solver.
6. **Impact.** +4–6% at PRA/PRApplied.

#### I2. The topography–charge degeneracy is not stated in quantum estimation form

Theorem `thm:degeneracy` is classical; §3.10 builds the multiparameter Fisher matrix only for topography modes $\{z_k\}$.
1. **Why it weakens the paper.** The claim "quantum enhancement cannot resolve the degeneracy" is never proved at the level the paper uses elsewhere.
2. **What is missing.** Corollary: for data generated by $H_x[\tilde z_s+\tilde\sigma/2\varepsilon_0kE_{\rm bg}]$, the two-parameter QFI matrix for $(z_s,\sigma)$ at fixed height has rank 1 per mode; the quantum CRB for the charge direction diverges; only a non-electrostatic channel or priors restore finiteness.
3. **Why reviewers request it.** It closes the loop between the identifiability and QFI sections in three lines.
4. **Proposed subsection.** `Corollary [Quantum Cramér–Rao bound for the degenerate pair]` after `thm:degeneracy` or in §3.10.
5. **Content.** Rank-1 argument + reference to `eq:qfi_matrix`.
6. **Impact.** +4–6%.

#### I3. The heating-rate channel has no quantitative spatial model (the degeneracy-breaker is unquantified)

The heating channel is invoked as the non-electrostatic channel (Theorem `thm:degeneracy`, strategy (ii), §4.3); the text now correctly concedes "whether the heating channel actually restores the rank requires an explicit microscopic model of the noise source, and the quantitative rank-2 joint analysis is deferred to future work."
1. **Why it weakens the paper.** The principal escape route from the central limitation remains a promissory note.
2. **What is missing.** The spatial kernel of the heating map: from $S_E\propto d^{-4}$ for a surface patch at lateral distance $\rho$, $\dot{\bar n}(x,y)\propto\int(h^2+\rho^2)^{-2}\,\mathcal S_{\rm patch}(x',y')\,dx'dy'$ — a convolution of width $\sim h$; hence (i) heating-channel spatial resolution $\sim h$; (ii) the two-channel rank-2 condition number as a function of $k$ and patch correlation length $\ell_c$; (iii) the crossover $\ell_c$ for statistically significant separation.
3. **Why reviewers request it.** Without it, the key resolution strategy is qualitative.
4. **Proposed subsection.** New subsection `Heating-rate channel: spatial kernel and rank-2 identifiability` in Methods or Results.
5. **Content.** Derivation from Eq. `eq:turchette` + near-field scaling; numerical rank-2 demonstration; link to Prediction 3's $\ell_c$ extraction.
6. **Impact.** +6–10% across journals.

#### I4. The logarithmic-derivative formula in Lemma `lem:noise_cutoff` carries an extra factor of $x$ (lines 2676–2677)

**Quote.** "$d\lambda_{\min}/d\ln\Theta=2\pi h/[x_{\rm hi}^2(x_{\rm hi}-2)]\approx 1.2\times10^{-7}$ m per $e$-fold … at $x_{\rm hi}=13.4$."

**Referee's check.** From $\lambda_{\min}=2\pi h/x_{\rm hi}$ with $x^2e^{-x}=\Theta$: $dx/d\ln\Theta=x/(2-x)$, whence $d\lambda_{\min}/d\ln\Theta=2\pi h/[x(x-2)]$ — one power of $x$, not two. At $x=13.4$, $h=40$ µm: $2.51\times10^{-4}/(13.4\times11.4)=1.65\times10^{-6}$ m per e-fold, $\approx13.7\times$ the quoted value (and the sign convention: the derivative is negative as $\Theta$ grows). The logarithmic-dependence conclusion is unaffected.
1–4. **Proposed subsection:** correct the formula and value in `lem:noise_cutoff` (one line + optional remark).
5. **Content.** One-line fix.
6. **Impact.** +3–5%.

#### I5. The two-parameter (topography + charge) inversion is never demonstrated numerically

The joint Bayesian formalism (Eq. `eq:joint_post_var`, Fig. `fig:joint_bayes`) analyses a single Fourier mode; no surface carrying both $z_s$ and $\sigma$ is ever reconstructed, even synthetically.
1. **Why it weakens the paper.** Strategy (iii) — the paper's principal software-level answer to the degeneracy — is never shown to work on any example.
2. **What is missing.** A synthetic two-parameter demonstration: known plateau + localised charge patch; joint GP(Laplace) inversion; separation quality vs. prior strength; the residual floor $\sqrt{p_\varphi}/E_{\rm bg}$.
3. **Why reviewers request it.** A claimed resolution strategy must be exercised.
4. **Proposed subsection.** New figure + paragraph in §3.11 (`Joint topography–charge inversion (synthetic)`).
5. **Content.** Simulation with the existing FFT chain + posterior MAP; correlation table vs. prior strength.
6. **Impact.** +5–8%.

#### I6. Two quoted discrepancy-principle values differ by a factor of 300 without an explicit dimensionality statement

**Locations.** Line 2195: "$\lambda_{\rm DP} \approx 3 \times 10^{-5}\max|H_x|$ for the synthetic test case"; line 2219 (figure paragraph): "the discrepancy principle … gives $\lambda_{\rm DP} \approx 9.0\times10^{-3}\max|H_x|$" on the 2D test surface.
1. **Why it weakens the paper.** A referee will notice two values for the same quantity differing by 300× in the same section. The reconciliation (the first is the 1D cross-section with noise norm $\propto\sqrt{N_{\rm pix}}$, the second the full 2D surface with norm $\propto N_{\rm pix}$) is never stated.
2. **What is missing.** One sentence: "the $3\times10^{-5}$ value refers to the 1D cross-section and the $9\times10^{-3}$ value to the full $N\times N$ surface; the factor $\sim\sqrt{N_{\rm pix}}$ reflects the pixel-number scaling of the noise norm in the discrepancy equation."
3. **Why reviewers request it.** Apparent internal contradiction.
4. **Proposed subsection.** Amend the two paragraphs of §3.11.
5. **Content.** One clarifying sentence at each location.
6. **Impact.** +2–4%.

#### I7. The azimuthal average $\langle A(\theta)\rangle=2\sqrt2/\pi\approx0.90$ is stated without derivation (§2.1)

Correct result, but the elliptic integral $\frac{1}{2\pi}\int_0^{2\pi}\sqrt{\cos^4\theta+\sin^4\theta}\,d\theta=2\sqrt2/\pi$ should appear as a one-line Lemma. **Impact:** +2–3%.

---

### OPTIONAL (would significantly increase impact)

- **O1. Two-dimensional exact Rayleigh solution** (bi-sinusoidal grating) — would extend the nonlinear Landweber/uniqueness results (`thm:landweber`, `prop:injective`) from the scalar amplitude to 2D.
- **O2. Finite-element validation** — FD (Fig. S11), BEM and Rayleigh are present; a FEM Laplace solve is the canonical additional solver family.
- **O3. Image-charge density formula** — §2.1 states the equivalence to a continuous image-charge distribution "proportional to the local curvature" without the formula; a Proposition giving the density to $\mathcal O(z_s)$ converts a heuristic into a theorem.
- **O4. Micromotion-channel sensitivity and QFI** — $\delta A_{\min}^{\rm(mm)}(k;h)$ and the two-channel (gradient + curvature) joint QFI are never computed.
- **O5. Uncertainty budget for roughness parameters** — the budget covers the point height $z$; a propagated $R_a/R_q$ budget would align the metrology section with the ISO discussion.
- **O6. Abstract harmonisation** — "$A_{\min}\sim1$ nm" vs 1.2 nm in text; "$r=0.62\to0.76$" vs the more precise $0.617\to0.776$.
- **O7. Submodularity proof sketch** for the adaptive-scanning $(1-e^{-1})$ guarantee (currently a parenthetical citation).

---

## Part B — Prioritised roadmap

**CRITICAL (required before submission)**
1. C6 — Ramsey formula missing $2\pi$ (four locations).
2. C1 — `tab:sensitivity` QESPM entries: definition vs. values (10⁹ mismatch).
3. C3 — null-space/injectivity paragraph (§2.1).
4. C4 — "exactly once" in the spectrum proof.
5. C5 — equilibrium-shift numbers.
6. C2 — `tab:comparison` "$10^{-2}$ nm" cell unqualified.

**IMPORTANT (strongly recommended)**
7. I1 — second-order frequency-shift kernel.
8. I2 — quantum CRB form of the degeneracy.
9. I3 — heating-channel spatial kernel.
10. I4 — `lem:noise_cutoff` derivative.
11. I5 — two-parameter synthetic inversion demo.
12. I6 — $\lambda_{\rm DP}$ 1D/2D clarification.
13. I7 — $\langle A(\theta)\rangle$ derivation.

**OPTIONAL (would significantly increase impact)**
14–20. O1–O7.

---

## Part C — Estimated acceptance probability (per cent)

| Journal | Before (current manuscript) | After CRITICAL+IMPORTANT | After ALL |
|---|---|---|---|
| Physical Review A | 45 | 58 | 62 |
| Physical Review Applied | 55 | 68 | 72 |
| New Journal of Physics | 60 | 72 | 76 |
| Measurement | 48 | 62 | 66 |
| Nature Communications | 9 | 13 | 15 |

*Basis.* The framing/overclaiming fixes of 30 August (sensitivity-vs-uncertainty distinction in the Abstract and §3.14, retitled ISO subsection, geometry paragraph, hedged heating and NOON language, corrected step-edge metrics, resolution vocabulary) are in place and correct — they raise the baseline above the previous round. What remains are the six CRITICAL arithmetic/rigor items (all independently recomputed today: C6 formula 2π, C1 table definition 10⁹, C5 equilibrium shift 10⁴/10⁸, C4 twice-attained values, C3 self-contradictory null-space paragraph, C2 unqualified table cell) and the seven IMPORTANT completeness items. Nature Communications additionally requires a broader-audience synthesis of the transfer-function-limited sensing principle; I2/I3 are the enablers.

---

## Part D — Verification note

Recomputed independently in this pass and found **correct**: ITF constants $1.51h$, $3.395/h$, $Q=0.589$, $C=1.92\times10^3$ m/s; four-regime cutoffs 18.7/16.0/14.4/13.7 µm and $\delta A_{\min}$ rows via $5.2\times10^{-4}\delta\omega e^8/4\times10^{10}$; $h_{\rm opt}=11.8$ µm; joint rank 36, Weyl slope 2.04 vs 2.22; prolate eigenvalues; MC coverage 95.1% / bias 0.7 nm; roadmap SNRs 24.3/563, $T_2$ 10/3.4/0.77/0.06 ms; uncertainty-table coefficients $c_h=kz$, $c_{E_{\rm bg}}=z/E_{\rm bg}$, $c_\omega=z/\omega$, $c_m=z/m$, $c_\sigma=1/(2\varepsilon_0kE_{\rm bg})$, Type-A $c_A=3.9\times10^{-11}$ m·s, and the Bayesian floor $\sqrt{p_\varphi}/E_{\rm bg}=319$ nm; step-edge metrics (RMS $2.3\times10^{-3}$, max $3.9\times10^{-3}$); Johnson-noise spectral density $k_BT\rho/4\pi d^3=2.1\times10^{-21}$ (V/m)²/Hz; patch-potential mapping $\phi=\sigma/2\varepsilon_0k=3.2$ mV at $k=2\times10^5$ rad/m.

---

## Part E — Resolution record (30 August 2026, applied)

All CRITICAL and IMPORTANT items were applied to main.tex (clean compile:
0 errors / 0 warnings, 72 pages):

- **C1 (tab:sensitivity)** — corrected values to $2.3\times10^{1}$ and
  $1.6\times10^{2}$ V/m²/√Hz with the definition $(2m\omega_x/e)\,\delta\omega$
  (rad/s), and "gradient" renamed "curvature" to match the observable
  $\partial_x^2\Phi$. *(Audit arithmetic correction: the Round-6 finding's
  quoted 2.3×10⁻⁸ was in error by 10⁹ — $(2m\omega_x/e)=5.21$, giving
  22.9 V/m²/√Hz; the published table value 72 was exactly $\pi$ times the
  definitional value, which is what made the values irreproducible.)*
- **C2 (tab:comparison)** — caption now states that the QESPM vertical
  entry is the QPN detection sensitivity at $k_{\rm opt}=2/h$ in the
  charge-free regime, with the realistic $u_c=319$ nm quoted.
- **C3 (§2.1)** — replaced the erroneous null-space/injectivity paragraph
  with Lemma `lem:nullspace` (correct kernel characterisation and
  injectivity of the joint operator only), with proof.
- **C4 (prop:spectrum)** — "exactly once" corrected to "each interior value
  twice, at the two roots $x_{\rm lo},x_{\rm hi}$ of Corollary `cor:qfactor`".
- **C5 (equilibrium shift)** — recomputed: $|\Delta\mathbf{r}|\approx
  1.3\times10^{-6}$ m and the relative correction $\sim2\times10^{-3}$,
  with the correction re-classified as a bounded Type-B systematic.
- **C6 (Ramsey formula)** — all four occurrences corrected to
  $1/(4\pi\tau\sqrt{N_{\rm rep}\langle n\rangle})$, with the one-line
  phase-to-frequency derivation added at the first occurrence.
- **I1** — new Proposition `prop:secondorder_dw` deriving the second-order
  frequency-shift kernel and the ratio $2(k_sA)e^{-k_sh}$ for a sinusoid,
  verified against the exact Rayleigh solver to ≤0.3% at $k_sh=2,4,8$.
- **I2** — new Corollary `cor:qcrb_degenerate`: the two-parameter QFI
  matrix of $(z_s,\sigma)$ at fixed height is rank-1 per mode and the
  quantum CRB diverges in the charge direction (with proof).
- **I3** — new "Heating-channel spatial kernel" paragraph in §2.6:
  $\dot{\bar n}\propto\int S_{\rm patch}(h^2+\rho^2)^{-3}\,d^2\rho$,
  total weight $\pi/2h^4$ (recovering $\alpha=4$), half-width $0.51h$,
  and the quadratic-in-charge character under the fluctuating-patch model
  (kernel integrals verified analytically and numerically).
- **I4 (lem:noise_cutoff)** — derivative corrected to
  $2\pi h/[x_{\rm hi}(x_{\rm hi}-2)]\approx1.65\times10^{-6}$ m per
  $e$-fold, with the sign convention stated.
- **I5** — new synthetic two-parameter joint inversion demonstration
  (Fig. `fig:joint_tc`, generated by `joint_topo_charge_demo.py`):
  naive single-channel $r=0.784$ (RMSE 184 nm) vs joint MAP $r_z=0.791$
  (RMSE 172.6 nm) with charge recovery $r_\sigma=0.926$; charge-equivalent
  peak 319 nm = the budget value.
- **I6 (λ_DP)** — both quoted values now carry explicit 1D-vs-2D
  statements (the $\sqrt{N_{\rm pix}}$ noise-norm factor).
- **I7** — $\langle A(\theta)\rangle$ corrected from $2\sqrt2/\pi\approx
  0.90$ (wrong) to $\frac{1}{2\pi}\int_0^{2\pi}\sqrt{\cos^4\theta+\sin^4\theta}
  \,d\theta=\frac{2}{\pi}E(1/2)\approx0.860$, verified by direct quadrature
  (integral $5.4026$).

Optional items O1–O7 remain open.

---

## Part F — L-curve figure bug (30 August 2026, follow-up)

The L-curve panel of Fig. `fig:picard` was rendering **empty**: the
solution norm was computed as `f/Hs` with `Hs` exactly zero on the
$k_x=0$ column of the lattice, producing `0/0 = NaN` in every point of the
curve (the published claim "$\lambda_{\rm L}\lesssim10^{-9}\max|H_x|$" was
an artefact of `argmax` over NaN). Fixed in `audit_fixes.py`
(`fz = Hs/(Hs²+λ²)` with zero guard) and the figure regenerated. Corrected
values: $\lambda_{\rm L}\approx5.1\times10^{-5}\max|H_x|$ (2D), with
$\lambda_{\rm GCV}\approx9.3\times10^{-4}$ and
$\lambda_{\rm DP}\approx9.0\times10^{-3}$; the empirical
$\lambda=10^{-4}\max|H_x|$ lies between the L-curve and GCV values. Both
text locations updated; compile clean (72 pp, 0/0/0/0/0).
