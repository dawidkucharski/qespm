# Single Trapped Ion Surface Topography Reconstruction

**A PhD-level research proposal in quantum metrology and surface science.**

---

## Project Structure

```
ion_surface_sensor/
├── README.md                       # This file
├── proposal.md                     # Comprehensive research proposal (all 6 tasks)
├── appendix_mathematical.tex       # Detailed mathematical derivations (LaTeX)
├── forward_model.py                # Numerical forward model & Tikhonov inversion
├── manuscript/                     # LaTeX manuscript directory (to be populated)
│   ├── main.tex                    # Main manuscript template
│   └── figures/                    # PDF vector graphics only
├── simulations/                    # Numerical experiment outputs
│   └── .gitkeep
└── literature/                     # Bibliography and notes
    └── .gitkeep
```

## Key Documents

| File | Contents |
|------|----------|
| `proposal.md` | Full 10-section proposal covering: literature review, mathematical model, inverse problem analysis, comparison with AFM/STM, experimental design, publication strategy, risk assessment, and timeline. |
| `appendix_mathematical.tex` | LaTeX document with complete derivations: Mathieu equations, boundary perturbation theory, Tikhonov regularisation, Bayesian hierarchical model, and quantum sensitivity limits. |
| `forward_model.py` | Python implementation of the forward model (topography → observables) using FFT-based Poisson kernel propagation, with Tikhonov inversion demonstration. |

## Quick Start: Forward Model Simulation

```bash
# Install dependencies
pip install numpy scipy matplotlib

# Run simulation
python forward_model.py
```

This generates:
- `forward_model_output.pdf` — surface topography, charge, and all motional observables
- `tikhonov_inversion_demo.pdf` — comparison of true vs. reconstructed surface potential

## Manuscript Preparation

All manuscripts use **LaTeX** with **PDF vector graphics only** (no raster formats for scientific figures). Figures should be generated as standalone `.pdf` files via matplotlib (`pdf` backend), TikZ/PGFPlots, or exported from vector-capable analysis tools.

## Key References

See `proposal.md` §10 for the core bibliography. The full `.bib` file will be maintained in `manuscript/references.bib`.

## Contact

Research proposal — internal working draft.
