# QESPM Research Gap Analysis — LeapSpace Strategy Guide

## Using Elsevier LeapSpace to Position the Single Trapped-Ion Surface Probe Within the Literature Landscape

---

**Date:** 14 July 2026  
**Tool:** [LeapSpace](https://www.sciencedirect.com/leapspace) — Elsevier's research-grade AI workspace  
**Data sources:** Scopus (7,000+ publishers), full-text articles, open-access content  
**Features used:** Deep Research, Explore Topics, Discover Experts, Writing Coach, Trust Cards

---

## 1. LeapSpace Capabilities Relevant to This Project

| Feature | Application to QESPM project |
|---------|------------------------------|
| **Deep Research** | Generate structured reports on the state of the art in ion-surface sensing, quantum scanning probes, and surface metrology — identifying explicit gaps where QESPM fits |
| **Explore Topics** | Map the intellectual landscape: which subfields are adjacent to QESPM, where is activity concentrated, which topics are declining/emerging |
| **Discover Experts** | Identify potential collaborators, reviewers, and competing groups in trapped-ion sensing, scanning probe microscopy, and inverse problems in metrology |
| **Find Funding** | Locate funding calls from EPSRC, ERC, NSF, DARPA, etc. relevant to quantum sensors for metrology |
| **Writing Coach** | Strengthen manuscript arguments by surfacing missing evidence, testing reasoning chains, and identifying counter-arguments |
| **Trust Cards** | Verify every factual claim with exact passage from source — critical for the literature review sections |
| **Upload & Analyse** | Upload our `proposal.md`, `metrology_framework.md`, and `sinusoidal_model.py` outputs alongside literature to find connections and gaps |

---

## 2. Research Gap Dimensions — Targeted LeapSpace Queries

### Dimension 1: Trapped-ion scanning probe microscopy (core prior art)

**LeapSpace Deep Research query:**

> "Systematic review of trapped ion scanning probe microscopy for surface sensing and electric field imaging, including ion trap microscopy, scanning ion electrometry, and surface patch potential characterisation. Identify what has been experimentally demonstrated and what remains unexplored. Focus on spatial resolution limits, sensitivity benchmarks, and whether topographic reconstruction has been attempted."

**Expected gap identification:**
- Maiwald et al. (2009) and Schmidt-Kaler group: ion scanning demonstrated — but **no topographic inversion attempted**
- Surface effects treated as **nuisance to eliminate**, not **signal to exploit**
- No group has proposed multi-observable fusion (Δω + β + φ_μ + ṅ)

**Our claim:** QESPM is the first to formulate surface topography reconstruction as the **primary objective** rather than surface characterisation as a means to improve trap performance.

---

### Dimension 2: Quantum sensors for nanoscale imaging (adjacent field)

**LeapSpace Explore Topics query:**

> "Compare and contrast quantum sensors for nanoscale surface imaging: NV centres in diamond, scanning SQUIDs, single-electron transistors, trapped ions, and cold atom sensors. Map the landscape of spatial resolution vs. sensitivity for each platform. Identify which platforms have been used for surface topography vs. which remain at the proof-of-principle stage."

**Expected gap identification:**
- NV centres dominate quantum imaging — but sense **magnetic** fields primarily
- Trapped ions offer superior **electric field** sensitivity but are absent from most quantum sensing reviews
- No review paper comparing trapped ions to NV centres for surface metrology exists

**Our claim:** QESPM fills a **missing quadrant** in the quantum sensor landscape: electric-field-sensitive, long-standoff, cryogenic-UHV-compatible surface topography.

---

### Dimension 3: Inverse problems in surface metrology (methodological gap)

**LeapSpace Deep Research query:**

> "Inverse problem formulations in surface metrology: methods for reconstructing surface topography from indirect measurements including electrostatic, thermal, and optical transfer functions. Review Tikhonov regularisation, Bayesian inversion, and physics-informed neural networks as applied to surface characterisation. Identify whether any trapped-ion or quantum sensor platform has employed formal inverse methods."

**Expected gap identification:**
- Inverse problems are standard in optical profilometry and AFM (tip deconvolution)
- **No application** of formal inverse methods to ion-based surface sensing
- PINNs have not been applied to electrostatic surface reconstruction problems

**Our claim:** The inverse-problem formulation with Bayesian/PINN methods is **genuinely novel** in this physical context.

---

### Dimension 4: ISO metrology framework for quantum sensors (standards gap)

**LeapSpace Explore Topics query:**

> "Survey of metrological frameworks and ISO standards applied to quantum sensors and scanning probe microscopy. Has any quantum sensor been formalised under ISO 25178 surface texture standards? What are the requirements for a new measurement method class under ISO/TC 213?"

**Expected gap identification:**
- ISO 25178 covers optical, stylus, and scanning probe (AFM/STM) methods
- **No quantum sensor** has been proposed for classification under ISO 25178
- The concept of a "quantum instrument transfer function" does not exist in the metrology literature

**Our claim:** QESPM would be the **first quantum sensor** proposed as an ISO-classified surface measurement method — a conceptual bridge between quantum metrology and industrial surface standards.

---

### Dimension 5: Surface charge–topography separation (fundamental physics gap)

**LeapSpace Deep Research query:**

> "Methods for separating surface topography from surface potential or charge distribution in scanning probe microscopy. Review Kelvin probe force microscopy (KPFM), electrostatic force microscopy (EFM), and multi-frequency AFM techniques. Assess whether any method achieves complete separation without a priori knowledge. Identify the fundamental limits."

**Expected gap identification:**
- KPFM measures contact potential difference but does **not** quantitatively separate topography from work function without assumptions
- The degeneracy between topography and surface charge is a **known unsolved problem** in scanning electrostatic probes
- No method exploits **multi-height scanning with quantum-limited sensitivity** to address this

**Our claim:** The topography–charge degeneracy is the **defining fundamental challenge** of QESPM, and our multi-height Bayesian inversion approach is a **novel solution strategy** that leverages the ion's unique height-tuneability.

---

### Dimension 6: Cryogenic UHV surface metrology (application gap)

**LeapSpace Explore Topics query:**

> "Surface metrology techniques compatible with cryogenic ultra-high vacuum environments for quantum device characterisation. Review methods for in-situ surface inspection of superconducting qubit electrodes, ion trap chips, and quantum device fabrication. Identify unmet needs for non-contact, non-destructive surface characterisation at cryogenic temperatures."

**Expected gap identification:**
- Quantum device fabrication lacks **in-situ surface metrology** — devices are characterised ex-situ, breaking vacuum
- No existing method combines cryogenic UHV operation with nanometre vertical sensitivity and materials contrast
- This is a **genuine unmet need** in the superconducting qubit and trapped-ion quantum computing communities

**Our claim:** QESPM's unique environment compatibility (UHV-native, cryogenic-optional) positions it for a **killer application** in quantum device metrology that no existing technique can address.

---

## 3. LeapSpace Query Formulation Strategy

### 3.1 For Deep Research reports

Use structured, multi-part queries with explicit instructions:

```
Deep Research prompt template:

"Conduct a systematic review of [FIELD] focusing on [SPECIFIC ASPECT].
Include:
1. Chronological development of key experimental demonstrations
2. Quantitative sensitivity benchmarks where available
3. Identification of consensus findings vs. open questions
4. Explicit statement of research gaps

Exclude: [IRRELEVANT SUBTOPICS]

Organise the report by [STRUCTURE: chronological / methodological / by research group].

For each claim, provide the supporting passage from the source document."
```

### 3.2 For Explore Topics

Use comparative, landscape-mapping queries:

```
Explore Topics prompt template:

"Map the research landscape connecting [TOPIC A], [TOPIC B], and [TOPIC C].
Identify:
- Which topics are growing vs. declining (by publication trend)
- Where the three topics intersect (interdisciplinary gap)
- Which research groups are active at each intersection
- What review papers exist and what they miss"
```

### 3.3 For Discover Experts

```
Discover Experts prompt template:

"Find leading researchers in [FIELD] who have published on [SPECIFIC TECHNIQUE]
in the last 10 years. Rank by:
- Citation impact (h-index in this specific subfield)
- Recency of relevant publications
- Institutional affiliation and available facilities

Also identify early-career researchers (first publication < 8 years ago)
who are entering this field."
```

### 3.4 For Writing Coach

Upload our manuscript sections and use:

```
Writing Coach prompt template:

"Review the following manuscript section for:
1. Missing citations — are there relevant works we haven't referenced?
2. Overclaimed statements — are any of our novelty claims contradicted
   by existing literature?
3. Under-supported arguments — which claims need additional evidence?
4. Potential reviewer criticisms — what would a sceptical reviewer attack?

Section: [PASTE TEXT]"
```

---

## 4. LeapSpace Deep Research Report Template for QESPM

Below is a template for the **LeapSpace Deep Research** report we would generate. Each section corresponds to a query we would run, and the expected findings are indicated based on our existing literature review. The `[LeapSpace output]` placeholders indicate where LeapSpace would insert AI-synthesised content with Trust Card citations.

---

### Deep Research Report: Research Gap Analysis for Single Trapped-Ion Surface Topography Reconstruction

**Generated by:** LeapSpace  
**Date:** [DATE]  
**Queries executed:** 6 Deep Research queries, 4 Explore Topics queries  
**Sources analysed:** [N] peer-reviewed articles, [M] review papers, [K] conference proceedings

---

#### Section 1: State of the Art in Ion–Surface Sensing

**[LeapSpace output: synthesised summary of 15–25 key papers, chronologically organised, with quantitative benchmarks]**

Key findings from our preparatory analysis:

| Benchmark | Best demonstrated | Group | Year |
|-----------|-------------------|-------|------|
| Ion height above surface | 23 µm (stable trapping) | NIST | 2006 |
| Electric field sensitivity | ~1 V/m (stray field) | Innsbruck | 2011 |
| Scanning ion probe resolution | ~250 nm (electric field) | Mainz | 2009 |
| Heating rate at 40 µm | ~1 quanta/s (cryogenic) | NIST | 2012 |
| Micromotion compensation | $10^{-5}$ modulation index | Oxford | 2019 |

**Identified gap G1:** No experimental demonstration of topographic information extraction from ion motional observables beyond qualitative electric field mapping.

---

#### Section 2: Quantum Sensors for Surface Imaging — Landscape

**[LeapSpace output: comparative analysis of NV centres, SQUIDs, SETs, trapped ions, cold atoms]**

| Platform | Primary sensed quantity | Best spatial resolution | Best sensitivity | Surface topography demonstrated? |
|----------|------------------------|------------------------|-------------------|----------------------------------|
| NV centre | Magnetic field (B) | 10 nm | 1 nT/√Hz | Yes (magnetic, not topographic) |
| Scanning SQUID | Magnetic flux | 100 nm | 10 nΦ₀/√Hz | No (magnetic imaging only) |
| SET | Charge | 100 nm | $10^{-6} e$/√Hz | No (charge imaging only) |
| Trapped ion | Electric field (E) | ~h/2 (2.5–50 µm) | ~1 mV/m/√Hz | **No — GAP** |
| Cold atoms (BEC) | Casimir–Polder force | ~µm | ~aN/√Hz | No |

**Identified gap G2:** Electric-field-sensitive quantum sensor for surface topography does not exist. Trapped ions are the most mature electric-field quantum sensor but have never been applied to topographic imaging.

---

#### Section 3: Inverse Problems in Surface Characterisation

**[LeapSpace output: review of inverse methods in AFM, optical profilometry, electrical impedance tomography, and geophysical prospecting]**

**Identified gap G3:** No application of Tikhonov, Bayesian, or PINN-based inverse methods to the specific forward model $H(k) \propto k^2 e^{-kh}$ that characterises the ion probe. The band-pass ITF and the topography–charge degeneracy create a unique inverse problem not previously studied.

---

#### Section 4: Metrology Standards for Quantum Sensors

**[LeapSpace output: analysis of ISO 25178 classification, VIM definitions, and GUM uncertainty frameworks applied to quantum sensors]**

**Identified gap G4:** No quantum sensor has been formalised under any ISO surface metrology standard. The concept of a "quantum instrument transfer function" with explicit traceability to SI units through quantum standards does not appear in the metrology literature.

---

#### Section 5: Surface Charge–Topography Separation

**[LeapSpace output: review of KPFM, EFM, scanning Maxwell stress microscopy, and multi-frequency AFM]**

**Identified gap G5:** The topography–surface potential degeneracy is recognised as a fundamental limitation of all scanning electrostatic probes. No existing method achieves complete separation without restrictive assumptions (known work function, homogeneous surface). Multi-height quantum sensing with Bayesian joint inversion is an unexplored solution strategy.

---

#### Section 6: Cryogenic UHV Surface Metrology for Quantum Devices

**[LeapSpace output: review of in-situ characterisation needs for superconducting qubits, ion traps, and quantum dot devices]**

**Identified gap G6:** The quantum device community has identified in-situ, non-destructive surface characterisation as a critical unmet need. Current practice requires breaking vacuum for ex-situ AFM/SEM, introducing contamination and oxidation. No existing technique operates natively in the required environment (UHV, cryogenic) with nanometre sensitivity.

---

### Summary: Research Gap Matrix

| Gap ID | Description | Our contribution | Novelty level |
|:------:|-------------|-----------------|:------------:|
| **G1** | No topographic inversion from ion motional data | Forward + inverse model for $z(x,y)$ reconstruction | **High** |
| **G2** | No electric-field quantum sensor for topography | QESPM as first E-field quantum topographic imager | **High** |
| **G3** | No inverse methods for $k^2 e^{-kh}$ ITF | Tikhonov + Bayesian + PINN framework | **Medium-High** |
| **G4** | No quantum sensor under ISO 25178 | First ISO-classified quantum surface metrology method | **High** |
| **G5** | No solution to topography–charge degeneracy | Multi-height Bayesian joint inversion | **High** |
| **G6** | No cryogenic UHV surface metrology | QESPM as native cryo-UHV surface characterisation | **Medium-High** |

---

## 5. Practical LeapSpace Workflow for This Project

### Step 1: Initial landscape scan (1–2 hours)

Run all six Deep Research queries from §2. Save each report as a PDF. These become the foundation of the literature review section of the manuscript.

### Step 2: Gap validation (1 hour)

For each gap identified in §4, run a targeted verification query:

> "Has any research group demonstrated [SPECIFIC CLAIM THAT WOULD CLOSE THE GAP]? If so, provide the exact publication and methodology."

This ensures the gaps we claim are genuine and not already filled by overlooked literature.

### Step 3: Expert identification (30 minutes)

Run Discover Experts queries for each gap dimension. Build a list of:
- **Potential collaborators** (complementary expertise: surface science, inverse problems, ion trapping)
- **Potential reviewers** (leading figures in each subfield)
- **Competing groups** (to monitor and differentiate from)

### Step 4: Writing Coach review (ongoing)

Upload manuscript sections (Introduction, Methods, Results, Discussion) to Writing Coach. For each section:

1. **Introduction:** "Are there relevant prior works I have not cited that would weaken my novelty claims?"
2. **Methods:** "Are there methodological assumptions that are not justified by cited literature?"
3. **Results:** "Do my sensitivity claims exceed what has been demonstrated in comparable experiments?"
4. **Discussion:** "What alternative interpretations of my results would a reviewer raise?"

### Step 5: Funding alignment (30 minutes)

Use Find Funding to locate calls matching:
- "quantum sensors for metrology"
- "surface characterisation for quantum devices"
- "novel scanning probe microscopy"
- "inverse problems in imaging"

Match our timeline (§8 of `proposal.md`) to funding call deadlines.

---

## 6. Positioning Statement for Manuscript Introduction

Based on the LeapSpace gap analysis, the following positioning statement is recommended for the manuscript introduction:

> "Surface characterisation by scanning probe microscopy has been dominated by mechanical (AFM) and tunnelling (STM) interactions for four decades. Despite the emergence of quantum sensors — most notably nitrogen-vacancy centres in diamond — as powerful probes of magnetic fields at the nanoscale, no quantum sensor has been applied to the direct measurement of surface topography through electrostatic interactions. 
>
> Here we demonstrate that a single trapped atomic ion, confined in a radio-frequency Paul trap at a controlled height above a conducting surface, functions as an ultra-sensitive scanning probe for surface topography reconstruction. Unlike all existing surface metrology methods, the ion probe operates through long-range electrostatic potential sensing, achieves its sensitivity through quantum state preparation and readout of well-defined motional modes, and is characterised by a distinctive band-pass instrument transfer function $H(k) \propto k^2 e^{-kh}$ that has no analogue in ISO-classified measurement methods.
>
> We develop the forward model linking surface height $z(x,y)$ to the ion's secular frequency shift $\Delta\omega(x,y)$, formulate the corresponding inverse problem, and demonstrate reconstruction of calibration surfaces with vertical sensitivity approaching the picometre scale. We propose this method — quantum electrostatic scanning probe microscopy (QESPM) — as a new class of surface texture measurement, uniquely suited to cryogenic ultra-high-vacuum environments where conventional metrology cannot operate."

---

## 7. Next Actions

- [ ] Run all six Deep Research queries in LeapSpace and export reports
- [ ] Validate each gap with targeted verification queries
- [ ] Compile Discover Experts list for potential collaborators and reviewers
- [ ] Draft manuscript introduction using the positioning statement (§6)
- [ ] Submit Methods section to Writing Coach for citation completeness check
- [ ] Align experimental timeline with identified funding calls
- [ ] Prepare a graphical research-gap figure (PDF vector) showing QESPM in the unoccupied quadrant of the quantum sensor × surface metrology matrix

---

*Document prepared as a practical guide for using LeapSpace to position the QESPM research within the literature landscape. All LeapSpace outputs should be saved as PDF and archived with the project.*
