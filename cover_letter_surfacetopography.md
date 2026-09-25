# Cover Letter — Surface Topography: Metrology and Properties

Dawid Kucharski\
Poznan University of Technology\
Poznań, Poland\
dawid.kucharski@put.poznan.pl · 24 September 2026

---

The Editor\
Surface Topography: Metrology and Properties\
IOP Publishing

Dear Editor,

I am pleased to submit the manuscript entitled **"Quantum Limits of Surface Metrology: Spatial Information Extraction by a Single Trapped Ion"** for consideration as a research paper in *Surface Topography: Metrology and Properties*.  The submission comprises a 30-page main manuscript and a 35-page Supplementary Material.

**Fit to the journal's scope.** STMP's scope states that the journal looks at surfaces "from the fundamental, applied and natural sciences, at any and all length scales", covers "the modelling, design and characterisation of modified surfaces", and aims "to present the measurement of topography of surfaces and interfaces", with submissions required to use "thorough surface characterization and state-of-the-art surface metrology". This manuscript is a fundamental, modelling-led contribution to exactly that aim. It develops, from first principles, the complete measurement framework for a new non-contact, height-tunable (multiscale) surface-characterisation modality in which a single trapped ion — already present in surface-electrode ion traps — acts as an electrostatic probe of the surface beneath it. The paper delivers the metrological core the scope asks for: a defined instrument transfer function with quantified resolution limits and recoverable spatial bandwidth, a GUM-based measurement uncertainty budget, an SI traceability chain, a domain-of-validity diagram, and a pathway toward classification within the ISO 25178 framework. Although the paper is theoretical, its validation is not self-referential: the transfer function is cross-checked against exact Rayleigh–Floquet, boundary-element and finite-difference electrostatic solvers and exercised on published AFM topography data, and five quantitative falsifiable predictions together with a costed four-phase roadmap specify how the measurement will be realised on existing hardware. The manuscript therefore occupies the modelling-and-characterisation segment of the scope today and supplies the measurement framework on which the projected experimental work will build.

**What the paper does.** Boundary-perturbation theory on Laplace's equation yields the instrument transfer function H_x(k; h) ∝ k_x^2 e^(−kh), a distinctive height-tunable band-pass transfer function whose e^(−kh) factor imposes an irreversible spatial-information cutoff: no quantum enhancement can recover spatial frequencies the physical filter has destroyed. Singular-value analysis quantifies the resulting information ceiling and its dependence on ion height, and a theorem proves the topography–charge identifiability limit (rank-1 joint measurement matrix). The manuscript explicitly separates quantum measurement sensitivity from topographic reconstruction uncertainty: the realistic uncertainty budget (combined standard uncertainty 319 nm under the stated residual-charge baseline; 99.99 % of the variance from surface charge) is set by electrostatic identifiability rather than by noise. Analytical scaling laws, a five-row comparison of charge-limited and charge-free regimes, Bayesian adaptive scanning, five quantitative falsifiable experimental predictions, and a costed experimental roadmap are provided.

**Rigour and validation.** The analytical transfer function is independently cross-checked against exact Rayleigh–Floquet, boundary-element and finite-difference electrostatic solvers (relative error 0.01 % at A = 100 nm, growing as 0.14(kA)²), closed-form analytical benchmarks (step-edge response), Monte Carlo uncertainty propagation with coverage checks, and evaluated using published atomic-force-microscopy data of a stainless-steel surface, where the reconstruction correlation improves from r = 0.62 at h = 40 µm to r = 0.78 at h = 2 µm. This evaluation uses experimentally measured topography as input and is not presented as an experimental validation of the instrument. The experimental demonstration itself is transparently future work, and the paper's claims are confined to the analytically derived, independently cross-checked framework.

**Status.** Single-author paper; not under consideration elsewhere. A 35-page Supplementary Material (detailed derivations, operator-theoretic analysis, solver hierarchy, Monte Carlo uncertainty propagation, robustness studies, calibration and traceability chain, and the experimental roadmap) accompanies the submission. All code is openly available in a public GitHub repository (https://github.com/dawidkucharski/qespm), together with a synthetic benchmark dataset, archived in a permanent Zenodo deposit (DOI: 10.5281/zenodo.22721000, https://doi.org/10.5281/zenodo.22721000). The validation data are publicly available (Mendeley Data), as stated in the manuscript.

**Suggested referees** (independent experts, no personal connection):
1. Rafał Demkowicz-Dobrzański, University of Warsaw — quantum estimation theory and quantum Fisher information bounds.
2. Tracy Northup, Universität Innsbruck — trapped-ion interactions with surfaces and surface charge.
3. Per Christian Hansen, Technical University of Denmark — discrete inverse problems and regularisation.
4. Otmar Scherzer, University of Vienna — regularisation theory for ill-posed problems.
5. John Bollinger, NIST — quantum-enhanced sensing with trapped-ion crystals.
6. Richard Leach, University of Nottingham — areal surface metrology and ISO 25178 standardisation.
7. Han Haitjema, KU Leuven — surface texture measurement, calibration and ISO standardisation.

I would be glad to address any questions and thank you for considering this work for *Surface Topography: Metrology and Properties*.

Sincerely,

Dawid Kucharski\
Poznan University of Technology, Poznań, Poland
