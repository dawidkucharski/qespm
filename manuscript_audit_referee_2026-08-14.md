# Referee Audit — Round 5 (from scratch, independent of previous rounds)

**Date:** 14 August 2026
**Auditor stance:** senior theoretical physicist / applied mathematician; refereeing for Physical Review A, Physical Review Applied, New Journal of Physics, Measurement, Nature Physics.
**Scope:** identify everything *missing* from the theoretical development; the scientific idea is assumed correct and unchanged. This audit does **not** edit the manuscript.
**Basis:** complete line-by-line reading of the current `main.tex` (68 pp, 45 refs) and `supplementary.tex` (13 pp, sections S1–S8) as they exist on disk today, plus independent numerical checks of every quantitative claim quoted below. All quotations are verbatim from the current files. Where I quote my own arithmetic, it is explicitly labelled "referee's check".

---

## Part A — Findings

### CRITICAL (required before submission)

#### C1. False injectivity statement and garbled null-space characterisation (§2.1, "Forward model", closing paragraph)

**Quote.** "The null space of $\mathcal{H}_h$ (for the $x$-channel alone) consists of all surface profiles with $\tilde{z}_s(k_x=0, k_y)=0$, i.e., profiles that depend only on $y$. Combining both channels reduces the null space to the constant offset (unobservable in any scanning measurement). $\mathcal{H}_h$ is injective on the subspace of zero-mean functions; the non-injectivity is limited to the physically trivial constant mode."

1. **Why it weakens the paper.** It contains two mathematical errors in a single paragraph of the operator-theory section, the part of the paper that will be read most carefully by the applied-mathematics referee. (i) The characterisation "profiles with $\tilde{z}_s(k_x=0,k_y)=0$" is the opposite of the correct statement: the symbol $H_x=-Ck_x^2e^{-kh}$ vanishes *on* the line $k_x=0$, so the null space is the set of profiles whose transform is **supported on** $\{k_x=0\}$ (zero off the line), i.e. y-only profiles — which is what the same paragraph itself says in the next clause. (ii) The final sentence then claims $\mathcal{H}_h$ is injective on zero-mean functions, directly contradicting the preceding sentence: a y-only profile with zero mean is zero-mean, yet lies in $\ker\mathcal{H}_h$. The correct statement (injectivity of the **joint** operator $\mathcal{H}_x\oplus\mathcal{H}_y$ on zero-mean functions) is proved later in Proposition `prop:joint`; the single-channel claim is simply false.
2. **What is missing.** A correct, numbered statement of the single-channel kernel (Lemma: $\ker\mathcal{H}_h=\{\tilde z_s:\,\operatorname{supp}\tilde z_s\subseteq\{k_x=0\}\}$), and a correct injectivity statement restricted to the joint operator.
3. **Why reviewers request it.** It is a self-contradiction in a theorem-strewn section; an applied-mathematics referee must reject or demand correction because it undermines confidence that the operator theory was checked.
4. **Proposed subsection.** Replace the closing paragraph of §2.1 with a `Lemma [Null space and injectivity]` carrying the exact statements above, with a one-line proof from the symbol.
5. **Content.** Statement + proof only (no new machinery).
6. **Impact.** Removing a demonstrable error: probability gain +8–10% at PRA/PRApplied; +5% elsewhere.

---

#### C2. False claim in the proof of the spectrum proposition ("takes every value … exactly once")

**Quote** (proof of Proposition `prop:spectrum`): "The symbol $H_x=-Ck_x^2e^{-kh}$ is real, continuous, and takes every value in $[-M,0]$ exactly once along the line $k_y=0$ ($k_x\in[0,\infty)$), so the essential range is the whole interval."

1. **Why it weakens the paper.** The function $k_x\mapsto Ck_x^2e^{-k_xh}$ on $[0,\infty)$ rises from 0 to its maximum $M$ at $k_x=2/h$ and decays back to 0: every value in $(0,M)$ is taken **twice**, not once (0 and $M$ once). The stated sentence is false as written, even though the *conclusion* (essential range $=[-M,0]$) is correct.
2. **What is missing.** The correct argument (the symbol is continuous with $\lim_{k_x\to0}f=\lim_{k_x\to\infty}f=0$ and attains $M$; hence by continuity the range is $[0,M]$; the word "once" should be deleted).
3. **Why reviewers request it.** A false assertion inside a proof is exactly what a theory referee marks as "proof not rigorous as written".
4. **Proposed subsection.** Correct the one sentence in the proof of `prop:spectrum`.
5. **Content.** One-sentence correction; optionally add "each interior value is attained twice, corresponding to the two roots $x_{\rm lo}, x_{\rm hi}$ of Corollary `cor:qfactor`" — which turns the error into a positive link to the passband analysis.
6. **Impact.** +5–8% at PRA/PRApplied; low elsewhere.

---

#### C3. Equilibrium-position-shift numbers are internally inconsistent by 3–4 orders of magnitude

**Quote** (§2.1, "Equilibrium position shift"): "For the nominal parameters ($h = 40$\,\textmu m, $|\Phi^{(1)}| \sim 100$\,\textmu V, feature scale $\sim 30$\,\textmu m), I find $|\Delta\mathbf{r}| \sim 10^{-10}$\,m and $|\Delta\mathbf{r}|/\lambda_{\rm feature} \sim 3 \times 10^{-6}$. The resulting correction to $\Delta\omega_x$ … is $\mathcal{O}(|\Delta\mathbf{r}|^2/\lambda_{\rm feature}^2) \sim 10^{-11}$ relative."

**Referee's check (independent arithmetic, labelled as such).** From the paper's own formula $\Delta\mathbf{r}=-(e/m\omega_x^2)\nabla\Phi^{(1)}$ with $m=40$ u, $\omega_x/2\pi=1$ MHz: $e/m\omega_x^2 = 6.1\times10^{-8}$ m per (V/m). For a feature of scale $\lambda_{\rm feature}=30$ µm carrying $|\Phi^{(1)}|\sim100$ µV at the ion, the gradient scale is $2\pi|\Phi^{(1)}|/\lambda_{\rm feature}\approx21$ V/m, giving $|\Delta\mathbf{r}|\approx1.3\times10^{-6}$ m — four orders of magnitude above the quoted $10^{-10}$ m. Even the cruder gradient $|\Phi^{(1)}|/\lambda_{\rm feature}\approx3.3$ V/m gives $2.0\times10^{-7}$ m (three orders). Consequently $(\Delta r/\lambda_{\rm feature})^2\sim10^{-3}$ (using $1.3\times10^{-6}$ m), not $10^{-11}$ as claimed: the stated correction is eight orders of magnitude too small.
1. **Why it weakens the paper.** Two quoted quantitative results in the Methods section are wrong by orders of magnitude under the manuscript's own nominal parameters. The *conclusion* ("correction negligible") survives — at the correct magnitude the correction is $\sim0.1\%$, still small — but a referee who repeats this three-line calculation will conclude the numbers were not checked.
2. **What is missing.** Correct values with the derivation shown (two lines), or a restated nominal configuration under which $10^{-10}$ m would actually hold (that would require $|\Phi^{(1)}|\sim10$ nV-scale gradients, contradicting the quoted 100 µV).
3. **Why reviewers request it.** Quantitative Methods claims are the referee's first arithmetic target.
4. **Proposed subsection.** Rewrite the "Equilibrium position shift" paragraph with the corrected displacement ($\sim10^{-6}$ m), the corrected relative correction ($\sim10^{-3}$), and the honest statement that the residual is $0.1\%$-level and absorbed into the Type-B budget (or negligible at the stated precision).
5. **Content.** Two-line computation; no new framework.
6. **Impact.** +5–8% across all target journals; prevents a "numbers not checked" desk impression.

---

#### C4. Ramsey frequency-uncertainty formula and its numerical values disagree by a factor of $2\pi$

**Quotes.** Type-A derivation (§3.14): "$u_A(f) = 1/(2\tau\sqrt{N_{\rm rep}\langle n\rangle})$; with $\tau = T_2/2 = 5$\,ms, $\langle n\rangle=0.5$ and $N_{\rm rep}=5.6\times10^3$ repetitions (28\,s of integration), $u_A(f) = 0.30$\,Hz". Prediction 5: "$\delta\omega_{\rm QPN} = 1/(2\tau\sqrt{N_{\rm rep}\langle n\rangle})$ … $\sim 0.7$\,Hz for ${}^{40}$Ca$^+$ … $N_{\rm rep}\approx10^3$".

**Referee's check.** The printed formula gives $1/(2\cdot5\,{\rm ms}\cdot\sqrt{5600\cdot0.5})=1.89$ Hz, not 0.30 Hz; and $1/(2\cdot5\,{\rm ms}\cdot\sqrt{1000\cdot0.5})=4.47$ Hz, not 0.7 Hz. The quoted values are reproduced exactly by the standard Ramsey SQL formula with the correct $2\pi$: $\delta f = 1/(4\pi\tau\sqrt{N_{\rm rep}\langle n\rangle})$ — i.e. $1.89/(2\pi)=0.301$ Hz and $4.47/(2\pi)=0.712$ Hz. The numbers were computed with the $4\pi\tau$ form; the printed formula has lost the $2\pi$ in both places where it appears.
1. **Why it weakens the paper.** The QPN floor feeds the entire four-regime table (`tab:regimes` rows 3–5), the roadmap SNRs, and the central "quantum enhancement" claims. A metrology referee will substitute the printed formula, get 4.5 Hz, and find the table inconsistent by 6.4×.
2. **What is missing.** The correct formula $1/(4\pi\tau\sqrt{N_{\rm rep}\langle n\rangle})$ (equivalently $\delta\phi=1/(2\sqrt{\langle n\rangle})$ per shot and $\delta\omega=\delta\phi/(2\pi\tau\sqrt{N_{\rm rep}})$), stated once and reused.
3. **Why reviewers request it.** Formula–value mismatch in the flagship quantitative table.
4. **Proposed subsection.** Correct both occurrences; add the one-line derivation ($\delta\phi_{\rm SQL}=1/2\sqrt{\langle n\rangle}$; phase–frequency conversion $2\pi\tau$).
5. **Content.** One equation + one sentence.
6. **Impact.** +8–10% at PRA/PRApplied/Measurement (metrology credibility); lower elsewhere.

---

### IMPORTANT (strongly recommended)

#### I1. The second-order frequency-shift kernel is never derived; only its numerical coefficient is quoted

The paper derives $\Phi^{(2)}$ (Eq. `eq:phi2`) and quotes the relative error law $0.14\,(kA)^2$ from the exact Rayleigh solver (§2.1, Lemma `lem:perturbation`), but never evaluates $\Delta\omega_x^{(2)}=(e/2m\omega_x)\,\partial_x^2\Phi^{(2)}|_{(x,y,h)}$ analytically.
1. **Why it weakens the paper.** A first-principles electrostatic treatment must carry the nonlinear correction through to the *observable*: the coefficient $0.14$ appears as a numerical accident rather than a derived constant.
2. **What is missing.** The Fourier form $\Delta\tilde\omega_x^{(2)}(k) = -\frac{eE_{\rm bg}}{2m\omega_x}\,e^{-kh}\,k_x^2\int \tilde z_s(k')\,|\mathbf k-\mathbf k'|\,\tilde z_s(\mathbf k-\mathbf k')\,\tfrac{d^2k'}{(2\pi)^2}$, its reduction for a single sinusoid (the $2k$ harmonic with explicit coefficient), and the proof that the Rayleigh $c_2$ coefficient equals this closed form — turning the $0.14$ into a *derived* number.
3. **Why reviewers request it.** Electrostatics referees expect the nonlinear term's analytic structure, including its sign and $h$-dependence.
4. **Proposed subsection.** `Proposition [Second-order frequency shift]` immediately after `prop:secondorder`, with corollary matching the Rayleigh coefficient.
5. **Content.** Two displayed equations; numerical cross-check vs the exact solver already available.
6. **Impact.** +4–6% at PRA/PRApplied.

---

#### I2. The topography–charge degeneracy is not stated in estimation-theoretic (quantum) form

Theorem `thm:degeneracy` proves classical rank-1 non-identifiability. The QFI section (§3.10) builds the multiparameter Fisher matrix only for the topography modes $\{z_k\}$.
1. **Why it weakens the paper.** The paper's central message is that quantum enhancement cannot help with the charge degeneracy — but the manuscript never proves it at the level it uses elsewhere: the two-parameter QFI matrix for $(z_s,\sigma)$ at a fixed height is singular (rank 1 per mode), hence the quantum Cramér–Rao bound is infinite for any unbiased joint estimator. Stated as a corollary, this makes the "no quantum fix" conclusion exact rather than rhetorical.
2. **What is missing.** Corollary: for data generated by $H_x\,[\tilde z_s+\tilde\sigma/2\varepsilon_0kE_{\rm bg}]$, the QFI matrix of $(z_s,\sigma)$ at one height has rank 1, its null direction is the charge direction; consequently $\delta A_{\rm charge}$ (CRB) diverges, and only a non-electrostatic channel or priors restore finiteness.
3. **Why reviewers request it.** It closes the loop between the two flagship sections (identifiability and QFI) and is three lines given the existing machinery.
4. **Proposed subsection.** `Corollary [Quantum Cramér–Rao bound for the degenerate pair]` after `thm:degeneracy` or in §3.10.
5. **Content.** Rank-1 argument + reference to `eq:qfi_matrix`; optional two-parameter SLD computation showing $L_{z_s}+L_{\sigma}$ degeneracy.
6. **Impact.** +4–6%.

---

#### I3. The heating-rate channel — the proposed degeneracy breaker — has no quantitative spatial model

The paper repeatedly invokes $\dot{\bar n}(x,y)$ as the non-electrostatic channel (Theorem `thm:degeneracy`, Corollary `cor:electrostatic`, strategy (ii), §4.3), and the Discussion concedes: "the quantitative rank-2 joint analysis with the heating channel is deferred to future work."
1. **Why it weakens the paper.** The entire escape route from the central limitation is asserted, not established: there is no formula relating the *local* heating rate to local surface parameters, so the claim that heating "provides a second, quasi-independent measurement" is unquantified.
2. **What is missing.** The spatial kernel of the heating channel: from $S_E(\omega_x)$ of a surface patch at lateral distance $\rho$, $\dot{\bar n}\propto(h^2+\rho^2)^{-2}$ (the $d^{-4}$ near-field law with $d=\sqrt{h^2+\rho^2}$), i.e. the heating map is the patch-fluctuation spectrum convolved with a kernel of width $\sim h$. This yields: (i) the spatial resolution of the heating channel $\sim h$; (ii) a two-channel rank-2 statement with the joint condition number of $(H_x\text{-map}, \text{heating-map})$ as a function of $k$ and of the patch correlation length $\ell_c$; (iii) the crossover $\ell_c$ at which separation becomes statistically possible at given noise.
3. **Why reviewers request it.** Without it, the key resolution strategy of the paper is a promissory note; the Discussion section's own admission makes the gap visible.
4. **Proposed subsection.** New subsection in Methods or Results: `Heating-rate channel: spatial kernel and rank-2 identifiability`, with a corollary-level statement of when the $(z_s,\sigma)$ pair becomes jointly identifiable.
5. **Content.** Derivation of the convolution kernel from Eq. `eq:turchette` + near-field scaling; numerical rank-2 demonstration on the lattice; connection to Prediction 3's $\ell_c$ extraction.
6. **Impact.** +6–10% across journals (turns the central practical claim from qualitative to quantitative).

---

#### I4. The logarithmic-derivative formula in Lemma `lem:noise_cutoff` carries an extra factor of $x$ and a value 14× too small

**Quote:** "the logarithmic derivative $d\lambda_{\min}/d\ln\Theta=2\pi h/[x_{\rm hi}^2(x_{\rm hi}-2)]\approx 1.2\times10^{-7}$ m per $e$-fold … at $x_{\rm hi}=13.4$."

**Referee's check.** From $\lambda_{\min}=2\pi h/x_{\rm hi}$ with $x_{\rm hi}$ solving $x^2e^{-x}=\Theta$: $dx/d\ln\Theta = x/(2-x)$ (negative for $x>2$, as $\lambda_{\min}$ decreases with increasing noise), whence $d\lambda_{\min}/d\ln\Theta = 2\pi h/[x(x-2)]$ — one power of $x$ in the denominator, not two. At $x=13.4$, $h=40$ µm: $2\pi h/(13.4\times11.4)=1.65\times10^{-6}$ m per e-fold, ≈13.7× the quoted $1.2\times10^{-7}$ m. The qualitative conclusion (logarithmic dependence) is unaffected; the quantitative derivative and its formula are wrong.
1. **Why it weakens the paper.** A quoted closed-form derivative is wrong by a factor $x_{\rm hi}$; referees of the regularisation-theory sections will differentiate and find the discrepancy.
2. **What is missing.** The correct formula and value ($1.65\times10^{-6}$ m per e-fold), and the sign convention (derivative negative as $\Theta$ grows).
3. **Why reviewers request it.** Derivative claims in lemmas are checkable in one line.
4. **Proposed subsection.** Correct the formula/value in `lem:noise_cutoff`.
5. **Content.** One-line fix; optional remark that $d\lambda_{\min}/d\ln\delta\omega$ inherits the same value.
6. **Impact.** +3–5%.

---

#### I5. The two-parameter (topography + charge) inversion is never demonstrated numerically

The Bayesian joint-inversion formalism (Eq. `eq:joint_post_var`, Fig. `fig:joint_bayes`) analyses a *single* Fourier mode; nowhere is a surface with both $z_s$ and $\sigma$ reconstructed, even synthetically. The forward-simulation section only notes that the mixed case "is nearly indistinguishable from the charge-only case."
1. **Why it weakens the paper.** Strategy (iii) — Bayesian separation with smooth/sparse priors — is the paper's principal software-level answer to the degeneracy, yet it is never shown to work on any example, not even a synthetic one.
2. **What is missing.** A synthetic two-parameter demonstration: a surface with a known Gaussian plateau + a localised charge patch; joint GP(Laplace) inversion; reported separation quality vs. prior strength; a map showing where the two contributions are (and are not) separable.
3. **Why reviewers request it.** A claimed resolution strategy must be exercised; referees will ask "does it actually work?".
4. **Proposed subsection.** New figure + paragraph in §3.11 (or a Results subsection `Joint topography–charge inversion (synthetic)`), with the honest statement of the residual degeneracy floor ($\sqrt{p_\varphi}/E_{\rm bg}$, Eq. `eq:joint_post_var`).
5. **Content.** Simulation with the existing FFT chain + posterior MAP; quantitative correlation table vs. prior strength.
6. **Impact.** +5–8%.

---

#### I6. The 2D anisotropy average $\langle A(\theta)\rangle=2\sqrt2/\pi\approx0.90$ is stated without derivation

(§2.1, "Anisotropy of the combined ITF".) The result is correct, but a theory paper should show the angular integral $\frac{1}{2\pi}\int_0^{2\pi}\sqrt{\cos^4\theta+\sin^4\theta}\,d\theta=2\sqrt2/\pi$, which is a standard elliptic integral; a one-line Lemma closes it. **Impact:** +2–3%; low-priority but zero-cost.

---

### OPTIONAL (would significantly increase impact)

- **O1. Exact two-dimensional Rayleigh solution (bi-sinusoidal grating).** The exact cross-validation (§4.4, supp S3) is 1D. A 2D Rayleigh-Floquet solution ($c_{nm}$ with $\sin(nk_xx)\sin(mk_yy)$) would extend the nonlinear Landweber/uniqueness results (`thm:landweber`, `prop:injective`) to two dimensions — currently the nonlinear theory is scalar ($A\in\mathbb R$).
- **O2. Finite-element validation.** The checklist's "finite-element validation" is absent; the FD panel (supp Fig. S11) covers a different discretisation family, but a FEM Laplace solve on the corrugated domain is the canonical third solver. FD+BEM+Rayleigh may suffice; add FEM only if a figure already exists.
- **O3. Image-charge density formula.** §2.1 states the equivalence to "a continuous distribution of image charges … proportional to the local curvature" without the formula. A Proposition giving the density $\varrho_{\rm im}(x,y)$ to $\mathcal O(z_s)$ would convert a heuristic into a theorem and is standard boundary-perturbation material.
- **O4. Micromotion-channel sensitivity and QFI.** $H_{\rm mm}\propto k_xe^{-kh}$ appears in `cor:electrostatic` and supp S3, but its noise-limited sensitivity $\delta A_{\min}^{\rm(mm)}(k;h)$ and its (two-channel, gradient+curvature) joint QFI are never computed. A short subsection would complete the multi-channel story quantitatively.
- **O5. Uncertainty budget for roughness parameters.** Table `tab:uncertainty` budgets the point height $z$; a propagated budget for $R_a$/$R_q$ (the ISO 25178 deliverable) would align the metrology section with the ISO discussion.
- **O6. Abstract harmonisation.** The abstract says $A_{\min}\sim1$ nm at the evaluation point while the text/table give 1.2 nm; and "$r=0.62\to0.76$" vs the text's more precise $0.617\to0.776$. Not errors, but referees do compare.
- **O7. Derivation of the adaptive-scanning submodularity guarantee.** The $(1-e^{-1})$ greedy bound is cited with a "Nemhauser–Wolsey-type" parenthetical; a one-paragraph statement (monotone submodularity of the GP information gain) with proof sketch would strengthen the Bayesian adaptive section.

---

## Part B — Prioritised roadmap

**CRITICAL (required before submission)**
1. C1 — null-space/injectivity statement (§2.1) — mathematical error.
2. C2 — "exactly once" claim in proof of `prop:spectrum`.
3. C3 — equilibrium-shift numbers ($10^{-10}$ m; $10^{-11}$ correction) recomputed ($\sim10^{-6}$ m; $\sim10^{-3}$).
4. C4 — Ramsey SQL formula missing $2\pi$ (printed $1/(2\tau\sqrt{N\langle n\rangle})$; values use $1/(4\pi\tau\sqrt{N\langle n\rangle})$).

**IMPORTANT (strongly recommended)**
5. I1 — analytic second-order frequency-shift kernel.
6. I2 — quantum (QFI) form of the degeneracy theorem.
7. I3 — quantitative spatial kernel of the heating channel (rank-2 identifiability).
8. I4 — correct logarithmic derivative in `lem:noise_cutoff` ($2\pi h/[x(x-2)]\approx1.65\times10^{-6}$ m).
9. I5 — synthetic two-parameter joint inversion demonstration.
10. I6 — derivation of $\langle A(\theta)\rangle=2\sqrt2/\pi$.

**OPTIONAL (would significantly increase impact)**
11–17. O1–O7 as above.

---

## Part C — Estimated acceptance probability (per cent)

| Journal | Before (current manuscript) | After CRITICAL+IMPORTANT | After ALL |
|---|---|---|---|
| Physical Review A | 40 | 55 | 60 |
| Physical Review Applied | 50 | 65 | 70 |
| New Journal of Physics | 55 | 70 | 75 |
| Measurement | 45 | 58 | 62 |
| Nature Communications | 8 | 12 | 15 |

*Basis for the estimates:* the manuscript is now a genuinely mature theory paper — compact-operator theory, Weyl asymptotics, minimax rates, QFI/Holevo analysis, metrology and roadmap are all present and mostly verified. The remaining CRITICAL items are localisable arithmetic/rigor errors (C1–C4) rather than missing pillars; fixing them removes the only reviewer-visible defects. The IMPORTANT items close the two remaining open ends of the framework (nonlinear electrostatics at the observable level, and the quantitative heating-channel story), which is what separates "very good" from "complete" at the theory journals. Nature Communications additionally requires a headline experimental advance or broader-audience synthesis; with I2/I3 in place the "transfer-function-limited quantum sensing" principle becomes quotable, which justifies the modest gain quoted.

---

## Part D — Verification note

All CRITICAL items were independently recomputed (see the arithmetic in each finding, labelled "referee's check"); none of the manuscript's other headline numbers were found inconsistent in this pass (the ITF constants $1.51h$, $3.395/h$, $Q=0.589$, $C=1.92\times10^3$ m/s; the four-regime cutoffs 18.7/16.0/14.4/13.7 µm; $h_{\rm opt}=11.8$ µm; joint rank 36 and slope 2.04; prolate eigenvalues; MC coverage 95.1%; roadmap SNRs 24/563 and $T_2$ values 10/3.4/0.77/0.06 ms — all reproduced independently today).

---

## Part E — Resolution record (30 August 2026)

An independent external review (framing/overclaiming audit) was processed
and its mandatory changes applied to main.tex and supplementary.tex:

1. "well-posed inversion" — verified **absent** from the manuscript (stale
   review item); the text already says "regularised inversion of a
   severely ill-posed inverse problem". No edit needed.
2. Sensitivity vs. reconstruction uncertainty — made explicit in the
   Abstract, in §3.1 ("Three distinct concepts": information sensitivity
   vs. identifiability), and in §3.14 (bold paragraph "Quantum measurement
   sensitivity is not topographic reconstruction uncertainty"; "the quantum
   sensor is not the fundamental limitation, electrostatic identifiability
   is").
3. Step-edge figure (supp. Fig. S14) — regenerated with the relative-error
   panel masked to $h/\sqrt3<|x|<3h$ (zero-crossing band shaded and
   excluded); relative RMS error $2.3\times10^{-3}$, maximum
   $3.9\times10^{-3}$; text in main and supplementary updated to match.
4. Resolution vocabulary — the undefined "${\sim}h/2$ lateral resolution"
   replaced in the limitations and the comparison table by the MTF cutoff
   $\lambda_{\min}^{\rm(MTF)}\approx1.51h$ (Rayleigh two-point criterion
   $\approx1.15h$), consistent with the step-edge FWHM $0.34h$.
5. Claim moderation — ISO section retitled "Potential classification
   within the ISO 25178 framework" with the explicit statement that no
   validation/calibration/inter-laboratory comparison has been performed;
   Introduction now frames the contribution as a "proposed theoretical
   framework ... experimental demonstration remains future work" and as a
   "candidate for a new measurement subclass"; the uniqueness claim
   "No other surface measurement method produces this spectral shape"
   removed ("a characteristic spectral signature predicted for the
   proposed ion probe"); "requires no fundamentally new technology"
   replaced by the integration-challenges sentence; NOON states presented
   as idealised (Heisenberg scaling only without decoherence and
   state-preparation overhead).
6. Geometry — new §2.1 paragraph "Model geometry versus experimental
   geometry" states that the half-space model describes the local planar
   surface region and that the planar-trap/laterally-translated-sample
   stage is arranged to approximate it.
7. Heating channel — hedged in three places: "potentially insensitive ...
   potentially providing", "requires an explicit microscopic model of the
   noise source and is not established here", and "The rank \emph{could}
   be restored only by a non-electrostatic channel ...".

Compile state after the changes: main 70 pp (0 errors / 0 undefined / 0
warnings / 0 textmu), supplementary 13 pp (clean).
