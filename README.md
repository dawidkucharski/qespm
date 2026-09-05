# Quantum Limits of Surface Metrology: Spatial Information Extraction by a Single Trapped Ion

Repository accompanying the manuscript

> **Quantum Limits of Surface Metrology: Spatial Information Extraction by a
> Single Trapped Ion**
> D. Kucharski, Poznan University of Technology
> (submitted to *Surface Topography: Metrology and Properties*, IOP)

## What this project does

The work develops and validates a new non-contact surface-characterisation
modality, QESPM, in which a single trapped atomic ion acts as an electrostatic
probe of the conducting surface beneath it:

- derives the instrument transfer function H_x(k; h) ∝ k_x² e^(−kh) from
  boundary-perturbation theory, and analyses its band-pass structure;
- establishes the irreversible spatial-information cutoff imposed by e^(−kh),
  the effective-rank scaling r_eff ∝ 1/h, and the theorem that multi-height
  scanning alone cannot separate surface topography from surface charge;
- builds the complete metrological framework: GUM-compliant uncertainty
  budget, SI-traceable calibration chain, ISO 25178 classification pathway,
  Bayesian and adaptive scanning strategies, and five quantitative falsifiable
  experimental predictions;
- validates the framework against an exact Rayleigh–Floquet solver, independent
  boundary-element and finite-difference solvers, Monte Carlo uncertainty
  propagation, and published AFM data of a stainless-steel surface
  (reconstruction correlation r = 0.62 → 0.76 as the ion height decreases from
  40 µm to 5 µm).

The experimental demonstration is future work, specified in a costed
four-phase roadmap (Supplementary Material, Sec. S9).

## Repository contents

```
manuscript/          main.tex, supplementary.tex, references.bib (LaTeX sources)
*.py                 scripts generating all figures and numerical analyses
manuscript_audit_*.md    referee-style audit records (verification history)
cover_letter_surfacetopography.md   submission cover letter
```

## Data sources (not included in this repository)

- AFM images and processing routine, Camargo Jr. (2019), Mendeley Data,
  CC BY 4.0 — doi:10.17632/6dzmrjngcg.3
- Surface-Topography Challenge benchmark dataset — doi:10.5281/zenodo.15341939

## Environment

Python 3 with numpy, scipy and matplotlib. LaTeX (pdflatex + bibtex) for the
manuscript.

## License

MIT — see the code package (`qespm_code_v1.0.zip`) for LICENSE and a detailed
README.
