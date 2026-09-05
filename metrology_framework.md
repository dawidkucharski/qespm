# Quantum Electrostatic Scanning Probe Microscopy (QESPM)

## A Metrology Framework for Single Trapped-Ion Surface Topography Measurement

---

**Document type:** Metrology framework specification  
**Status:** Working draft — 14 July 2026  
**Intended standard alignment:** ISO 25178 (Areal surface texture), GUM (JCGM 100:2008), VIM (JCGM 200:2012)

---

## 1. Measurand

### 1.1 Primary measurand

The **primary measurand** is the local surface height $z(x,y)$ of a conducting or semiconducting surface, defined as the vertical displacement of the surface equipotential from a reference mean plane, measured at lateral coordinates $(x,y)$ in the trap reference frame.

Per VIM (International Vocabulary of Metrology, JCGM 200:2012, §2.3), the measurand is the quantity intended to be measured. Here:

$$\boxed{z(x,y) = \mathcal{M}^{-1}\!\big[\Delta\omega_x, \Delta\omega_y, \beta, \phi_\mu, \dot{\bar{n}}\big](x,y)}$$

where $\mathcal{M}^{-1}$ denotes the regularised inverse of the forward measurement model $\mathcal{M}$ (the solution of the inverse electrostatic boundary-value problem, §3).

### 1.2 Measurement model (forward)

The complete forward measurement model maps the surface height $z(x,y)$ and surface charge distribution $\sigma(x,y)$ to the observables:

$$\mathbf{d} = \mathcal{F}(z, \sigma; \boldsymbol{\theta}) + \boldsymbol{\varepsilon}$$

where:
- $\mathbf{d} = \{\Delta\omega_x(x_i,y_j), \Delta\omega_y, \beta_x, \beta_y, \phi_\mu, \dot{\bar{n}}\}$ — the data vector
- $\mathcal{F}$ — the forward operator (composition of boundary-value problem for Laplace's equation and perturbation theory on the Mathieu equations)
- $\boldsymbol{\theta} = \{E_{\text{bg}}, h, \omega_{\text{sec}}, T, \rho, \ldots\}$ — instrument parameters
- $\boldsymbol{\varepsilon}$ — measurement noise (quantum projection noise, thermal fluctuations, detection noise)

For the **linearised regime** (small roughness, $|z| \ll h$, known $\sigma$), the forward model reduces to a linear Fredholm integral equation of convolution type:

$$\Delta\omega_x(x,y) = \frac{e E_{\text{bg}}}{2m\omega_x} \frac{\partial^2}{\partial x^2} \iint \frac{h \, z(x',y')}{[(x-x')^2 + (y-y')^2 + h^2]^{3/2}} \,dx'\,dy'$$

In Fourier space:

$$\boxed{\Delta\tilde{\omega}_x(k_x, k_y) = -\frac{e E_{\text{bg}}}{2m\omega_x} \, k_x^2 \, e^{-kh} \, \tilde{z}(k_x, k_y)}$$

### 1.3 Influence quantities (VIM §2.52)

| Quantity | Symbol | Typical range | Sensitivity to $z$ |
|----------|--------|---------------|---------------------|
| Surface charge density | $\sigma(x,y)$ | $0$–$10^{-2}$ C/m² | **Dominant confounding quantity** |
| Ion height | $h$ | $1$–$100$ µm | $\partial z/\partial h \propto z/h$ |
| Background field | $E_{\text{bg}}$ | $10^3$–$10^5$ V/m | $\partial z/\partial E_{\text{bg}} \propto -z/E_{\text{bg}}$ |
| Surface temperature | $T$ | $4$–$300$ K | Affects $\dot{\bar{n}}$ channel |
| Surface resistivity | $\rho$ | material-dependent | Affects $\dot{\bar{n}}$ channel |
| Secular frequency | $\omega_{\text{sec}}$ | $2\pi \times (0.5$–$5)$ MHz | $\partial z/\partial\omega \propto z/\omega$ |

---

## 2. Instrument Transfer Function (ITF)

### 2.1 Definition

Following ISO 25178-600 (§3.1), the instrument transfer function describes the ratio of the measured output to the true input as a function of spatial frequency. For the ion probe in the $x$-measurement channel:

$$\boxed{H_x(k_x, k_y; h) \equiv \frac{\Delta\tilde{\omega}_x(k_x, k_y)}{\tilde{z}(k_x, k_y)} = -\frac{e E_{\text{bg}}}{2m\omega_x} \, k_x^2 \, e^{-k h}}$$

where $k = \sqrt{k_x^2 + k_y^2}$.

### 2.2 ITF magnitude (amplitude transfer function)

$$|H_x(k_x, k_y; h)| = \frac{e E_{\text{bg}}}{2m\omega_x} \, k_x^2 \, e^{-k h}$$

The **normalised ITF** (modulation transfer function analogue):

$$\text{MTF}(k_x, k_y; h) = \frac{|H_x(k_x, k_y; h)|}{\max_{k_x,k_y}|H_x(k_x, k_y; h)|}$$

The maximum occurs at $k_x = \sqrt{2}/h$, $k_y = 0$ (or equivalently $k = \sqrt{2}/h$):

$$\max|H_x| = \frac{e E_{\text{bg}}}{2m\omega_x} \cdot \frac{2}{h^2} \cdot e^{-\sqrt{2}} \approx \frac{e E_{\text{bg}}}{2m\omega_x} \cdot \frac{0.486}{h^2}$$

### 2.3 Cut-off wavelengths

**Lower cut-off** (short-wavelength limit, high spatial frequency):

The 50% MTF point (3 dB attenuation from peak) occurs at $k_{\text{high}} \approx 3.4/h$:

$$\lambda_{\min} = \frac{2\pi}{k_{\text{high}}} \approx 1.85\,h$$

This is the **fundamental lateral resolution limit** — surface features with wavelength smaller than $\sim 2h$ are attenuated by more than 50%.

**Upper cut-off** (long-wavelength limit):

Defined by the $k_x^2$ factor: at low spatial frequencies, sensitivity vanishes quadratically. The 50% point (from the low-$k$ side) occurs where $k_x^2/(2/h^2) = 0.5$, giving $k_x = 1/h$. However, the practical upper wavelength limit is set by the **field of view**:

$$\lambda_{\max} = \text{FOV} = N \cdot \Delta x$$

where $N$ is the number of scan points and $\Delta x$ the sampling interval.

### 2.4 Anisotropy

The ITF is **anisotropic**:

$$H_x(k_x=0, k_y) = 0 \quad \forall k_y$$

The $x$-channel is completely blind to surface features that vary only in $y$. Full isotropy requires the **dual-channel ITF**:

$$\mathbf{H}(k_x, k_y) = -\frac{e E_{\text{bg}}}{2m} e^{-kh} \begin{bmatrix} k_x^2/\omega_x \\ k_y^2/\omega_y \end{bmatrix}$$

Combining both channels (and micromotion gradient) yields quasi-isotropic response.

### 2.5 Height dependence

The ITF depends parametrically on the ion–surface distance $h$. This is both a challenge (requires accurate $h$ knowledge) and an **opportunity**: multi-height measurements provide depth sectioning analogous to confocal microscopy, enabling separation of topography from charge ($k$-dependent scaling differs).

### 2.6 Phase transfer function (PTF)

The ITF is purely real for the $x$-channel (in-phase response for sine modes). The sign convention gives a $\pi$ phase shift ($\Delta\omega_x \propto -\sin kx$), corresponding to the second-derivative operator. The PTF is:

$$\phi_H(k_x, k_y) = \begin{cases} \pi & k_x \neq 0 \\ \text{undefined} & k_x = 0 \end{cases}$$

meaning the reconstructed surface is sign-inverted, which is trivially corrected.

---

## 3. Spatial Resolution

### 3.1 Lateral resolution

The lateral resolution is characterised by three metrics following ISO 25178-600:

**Spatial sampling interval:**
$$\Delta x_{\text{sample}} = \frac{v_{\text{scan}}}{f_{\text{measure}}}$$

where $v_{\text{scan}}$ is the stage velocity and $f_{\text{measure}}$ is the per-pixel measurement rate. Typical values: $v_{\text{scan}} = 10$ µm/s, $f_{\text{measure}} = 100$ Hz → $\Delta x_{\text{sample}} = 0.1$ µm.

**Rayleigh resolution criterion** (minimum resolvable separation):

Two Gaussian surface features of equal height separated by distance $d$ are resolved when the dip between their $\Delta\omega_x$ signals is detectable above noise. This gives:

$$\delta x_{\text{Rayleigh}} \gtrsim \frac{h}{2}$$

For $h = 40$ µm: $\delta x_{\text{Rayleigh}} \approx 20$ µm; for $h = 5$ µm: $\delta x_{\text{Rayleigh}} \approx 2.5$ µm.

**Sparrow resolution limit** (disappearance of central dip):

$$\delta x_{\text{Sparrow}} \approx \frac{h}{\sqrt{2}}$$

### 3.2 Vertical resolution

The **noise-equivalent height** (minimum detectable amplitude at a given spatial frequency):

$$\delta z(k) = \frac{2m\omega_x}{e E_{\text{bg}}} \frac{\delta\omega_{\text{min}}}{k^2 e^{-k h}}$$

At the optimal spatial frequency $k_{\text{opt}} = 2/h$:

$$\boxed{\delta z_{\min} = \frac{2m\omega_x}{e E_{\text{bg}}} \frac{h^2}{4e^{-2}} \delta\omega_{\text{min}}}$$

Numerical evaluation for ${}^{40}\text{Ca}^+$, $\omega_x/2\pi = 1$ MHz, $E_{\text{bg}} = 10^4$ V/m, $\delta\omega_{\min}/2\pi = 1$ Hz:

| $h$ [µm] | $\delta z_{\min}$ [nm] |
|:---------:|:----------------------:|
| 1 | $4\times 10^{-7}$ |
| 5 | $1\times 10^{-5}$ |
| 10 | $4\times 10^{-5}$ |
| 40 | $6\times 10^{-4}$ |
| 100 | $4\times 10^{-3}$ |

The vertical resolution is **extraordinary** — sub-picometre at close approach — but degrades as $h^2$ with ion height. This is many orders of magnitude finer than the lateral resolution; the ion probe is a "vertically acute, laterally myopic" instrument, complementary to optical methods.

### 3.3 Resolution-area product

A useful figure of merit for comparing surface metrology instruments:

$$Q = \frac{\delta z_{\min}}{\delta x_{\text{Rayleigh}}} \approx \frac{4m\omega_x}{e E_{\text{bg}}} \frac{h}{2e^{-2}} \delta\omega_{\min} \propto h$$

For $h = 40$ µm: $Q \approx 3 \times 10^{-5}$ — exceptionally low, indicating extreme vertical sensitivity per unit lateral resolution.

---

## 4. Uncertainty Budget (GUM-Compliant)

### 4.1 Measurement function

The measurement result is:

$$z(x,y) = \frac{1}{E_{\text{bg}}} \cdot \mathcal{D}^{-1}\!\left[-\frac{2m\omega_x}{e} \Delta\omega_x^{\text{corr}}\right]$$

where $\mathcal{D}^{-1}$ denotes deconvolution of the $e^{-kh}$ propagation, and $\Delta\omega_x^{\text{corr}}$ is the frequency shift corrected for charge background.

### 4.2 Uncertainty sources

#### Type A (evaluated by statistical methods)

| Symbol | Source | Evaluation method | Typical $u_A$ |
|--------|--------|-------------------|---------------|
| $u_A(\Delta\omega)$ | Frequency shift repeatability | Standard deviation of $N$ repeated measurements at each pixel | $2\pi \times 0.3$ Hz (for $\tau = 100$ ms) |
| $u_A(x,y)$ | Ion position localisation | Fluorescence centroid fitting | $5$ nm |
| $u_A(\beta)$ | Micromotion amplitude | Photon-correlation statistics | $10^{-5}$ modulation index |

**Quantum limit** (Cramér–Rao lower bound):

$$u_A^{\text{QCRB}}(\Delta\omega) = \frac{1}{2\eta\sqrt{N_{\text{rep}} \tau T_2}} \approx 2\pi \times 0.01\ \text{Hz}$$

for $\eta = 0.1$ (Lamb-Dicke parameter), $N_{\text{rep}} = 1000$, $\tau = 1$ ms, $T_2 = 10$ ms.

#### Type B (evaluated by other means)

| Symbol | Source | Probability distribution | Standard uncertainty $u_B$ | Sensitivity $\partial z/\partial x_i$ |
|--------|--------|--------------------------|----------------------------|---------------------------------------|
| $u_B(E_{\text{bg}})$ | Background field calibration | Normal | $10$ V/m ($10^{-3}$ relative) | $-z/E_{\text{bg}}$ |
| $u_B(h)$ | Ion height | Normal | $100$ nm | $z/h$ |
| $u_B(\sigma)$ | Uncompensated surface charge | Uniform (bounded) | $\sigma_{\text{residual}} \cdot d_{\text{eff}}/\varepsilon_0$ | $1/E_{\text{bg}} \cdot e^{+kh}/k^2$ |
| $u_B(\omega_x)$ | Secular frequency calibration | Normal | $2\pi \times 10$ Hz | $z/\omega_x$ |
| $u_B(m)$ | Ion mass | Normal | negligible ($<10^{-9}$ relative) | $-z/m$ |
| $u_B(T)$ | Temperature uncertainty | Normal | $0.1$ K | material-dependent |
| $u_B(\text{drift})$ | Long-term drift | Uniform | $5$ nm/min (stage), $1$ V/m h (fields) | time-dependent |

### 4.3 Combined standard uncertainty

$$u_c^2(z) = \sum_i \left(\frac{\partial z}{\partial x_i}\right)^2 u^2(x_i) + \sum_{i<j} \frac{\partial z}{\partial x_i}\frac{\partial z}{\partial x_j} u(x_i, x_j)$$

Assuming uncorrelated inputs, the **dominant terms** in decreasing order of magnitude (for $h = 40$ µm, $A = 100$ nm, $k = 2 \times 10^5$ rad/m):

| Rank | Source | Contribution to $u_c(z)$ | Fraction of total variance |
|:----:|--------|:------------------------:|:--------------------------:|
| 1 | **Surface charge** $\sigma$ | $5$–$50$ nm | $80$–$95\%$ |
| 2 | Ion height $h$ | $0.5$ nm | $2$–$5\%$ |
| 3 | Background field $E_{\text{bg}}$ | $0.3$ nm | $1$–$3\%$ |
| 4 | Frequency measurement (Type A) | $0.01$ nm | $<0.1\%$ |
| 5 | Secular frequency $\omega_x$ | $<10^{-3}$ nm | negligible |
| 6 | Ion mass $m$ | $<10^{-6}$ nm | negligible |

**Critical finding**: The uncertainty budget is **dominated by surface charge**, not by measurement noise. The ion probe's exquisite frequency sensitivity is not the limiting factor — the inability to perfectly separate topography from charge is. All other sources contribute at the sub-nanometre level.

### 4.4 Expanded uncertainty

$$U = k \cdot u_c(z), \quad k = 2 \ (\text{95\% confidence})$$

For the nominal case: $U \approx 10$–$100$ nm, dominated by the charge systematic. With multi-height correction or charge-control procedures, $U \lesssim 1$ nm is achievable.

### 4.5 Uncertainty reporting (per ISO 14253-1)

Measurement result: $z(x,y) \pm U$ (95% confidence, $k=2$)

Example: "The step height measured by QESPM is $198.3 \ \text{nm} \pm 6.2 \ \text{nm}$ ($k=2$), where the dominant uncertainty contribution is uncompensated surface charge ($u_B = 3.0 \ \text{nm}$)."

---

## 5. Metrological Traceability

### 5.1 Traceability chain

Following ISO 25178-601 and the GUM concept of traceability:

```
┌─────────────────────────────────────────────────────────────────┐
│                     SI BASE UNITS                               │
│  metre (m)    second (s)    kilogram (kg)    ampere (A)         │
└──────┬────────────┬──────────────┬───────────────┬──────────────┘
       │            │              │               │
       ▼            ▼              ▼               ▼
┌──────────┐ ┌──────────┐ ┌───────────┐ ┌──────────────────┐
│ Laser    │ │ GPS-dis- │ │ Atomic    │ │ Josephson voltage │
│ inter-   │ │ ciplined │ │ mass      │ │ standard (NIST/   │
│ ferometer│ │ clock    │ │ tables    │ │ PTB)              │
└────┬─────┘ └────┬─────┘ └─────┬─────┘ └────────┬─────────┘
     │            │             │                │
     ▼            ▼             ▼                ▼
┌──────────┐ ┌──────────┐ ┌───────────┐ ┌──────────────────┐
│ Nano-    │ │ rf       │ │ Ion mass  │ │ Trap electrode   │
│ position │ │ frequency│ │ m(Ca⁺)    │ │ voltages V_dc,   │
│ stage    │ │ Ω, ω_sec │ │           │ │ V_rf             │
│ (x,y)    │ │          │ │           │ │                  │
└────┬─────┘ └────┬─────┘ └─────┬─────┘ └────────┬─────────┘
     │            │             │                │
     └────────────┼─────────────┼────────────────┘
                  │             │
                  ▼             ▼
          ┌──────────────────────────┐
          │  ION MOTIONAL OBSERVABLES│
          │  Δω_x, Δω_y, β, φ_μ, ṅ  │
          └────────────┬─────────────┘
                       │
                       ▼
          ┌──────────────────────────┐
          │  INVERSE MODEL           │
          │  (calibrated via         │
          │   sinusoidal standard)   │
          └────────────┬─────────────┘
                       │
                       ▼
          ┌──────────────────────────┐
          │  SURFACE HEIGHT z(x,y)   │
          │  Traceable to SI metre   │
          │  through calibrated ITF  │
          └──────────────────────────┘
```

### 5.2 Calibration hierarchy

| Level | Standard | Role | Uncertainty ($k=2$) |
|:-----:|----------|------|:--------------------:|
| **Primary** | Sinusoidal grating artefact (PTB/NIST) | ITF amplitude calibration | $0.5$ nm (step height) |
| **Secondary** | Ion response to known $E_{\text{bg}}$ | Field calibration | $10$ V/m |
| **Working** | In-situ clean gold surface | Daily $h$ verification | $100$ nm |
| **Check** | Repeated measurement of stable feature | Drift monitoring | Type A |

### 5.3 Key traceability challenges

1. **Electrostatic–geometric transfer**: The conversion from $\Phi$ to $z$ requires knowledge of $E_{\text{bg}}$, which depends on trap geometry and applied voltages. Traceability to the SI volt through the Josephson effect is well established, but the geometric factor $E_{\text{bg}}/V_{\text{applied}}$ requires calibration against a known height standard.

2. **Height $h$**: Direct interferometric measurement of the ion–surface gap is challenging in a cryogenic UHV environment. Electrical methods (capacitance, induced-charge detection) provide an alternative but require modelling of stray capacitances.

3. **Surface charge**: There is **no primary standard for surface charge distribution** at the relevant spatial scales. This is a fundamental traceability gap that currently limits the ion probe (and indeed all scanning electric-field probes) to relative rather than absolute measurements of unknown surfaces.

---

## 6. Calibration Procedure

### 6.1 Overview

The calibration establishes the relationship between the measured observables and the surface height through determination of the instrument parameters $\boldsymbol{\theta}$ and verification of the ITF.

### 6.2 Step 1: Trap field characterisation

**Procedure:**
1. Trap a single ${}^{40}\text{Ca}^+$ ion and Doppler-cool to $\sim 1$ mK.
2. Measure secular frequencies $\omega_x, \omega_y, \omega_z$ by tickle spectroscopy (precision $\sim 10$ Hz).
3. Vary dc electrode voltages systematically ($V_{\text{dc}} \to V_{\text{dc}} + \delta V$) and measure the resulting frequency shifts and equilibrium position shifts.
4. Fit to the Mathieu equation model to extract $E_{\text{bg}}$ and the geometric efficiency factors.

**Deliverable:** $E_{\text{bg}}$ with $u(E_{\text{bg}}) \leq 10$ V/m ($10^{-3}$ relative).

### 6.3 Step 2: Height calibration

**Procedure:**
1. Position the calibration artefact surface using the nanopositioning stage.
2. Approach the surface slowly (with continuous ion fluorescence monitoring) until the ion signal indicates proximity (increased heating rate, stray field onset).
3. Record the stage position at contact (ion loss) as $h = 0$ reference — or use a non-destructive electrical touch detection.
4. Retract to working height; verify $h$ via electrical capacitance measurement (ion–surface image charge).
5. Repeat at multiple heights to establish stage calibration and $h$ uncertainty.

**Deliverable:** $h$ with $u(h) \leq 100$ nm.

### 6.4 Step 3: ITF calibration (primary calibration)

**Procedure:**
1. Use a **calibrated sinusoidal grating artefact** with known amplitude $A_{\text{ref}}$ and wavelength $\lambda_{\text{ref}}$, certified by a national metrology institute (PTB, NIST, NPL) with traceability to the SI metre.
2. Position the grating at height $h$.
3. Scan the ion across the grating in the $x$-direction, measuring $\Delta\omega_x(x)$.
4. Fit the measured $|\Delta\omega_x|$ vs. $k_{\text{ref}} = 2\pi/\lambda_{\text{ref}}$ to the ITF model:

$$|\Delta\omega_x(k_{\text{ref}}; h)| = \frac{e E_{\text{bg}}}{2m\omega_x} A_{\text{ref}} k_{\text{ref}}^2 e^{-k_{\text{ref}} h}$$

5. Extract the **calibration factor** $C = eE_{\text{bg}}/(2m\omega_x)$ by comparing measured and predicted signal.
6. Repeat for multiple $h$ values and grating wavelengths to verify the $k^2 e^{-kh}$ scaling.
7. Document the **calibrated ITF** with uncertainty bounds.

**Deliverable:** Calibrated ITF $H_x^{\text{cal}}(k; h)$ with uncertainty envelope.

### 6.5 Step 4: Charge background characterisation

**Procedure:**
1. Perform a **reference scan** at a height $h_{\text{ref}} \gg h_{\text{meas}}$ (e.g., $h_{\text{ref}} = 200$ µm while $h_{\text{meas}} = 40$ µm).
2. At this height, $e^{-k h_{\text{ref}}} \ll 1$ for all but the very lowest spatial frequencies — the topographic signal is negligible.
3. The residual $\Delta\omega_x(x,y)$ map reflects the surface charge distribution $\sigma(x,y)$.
4. Model $\sigma(x,y)$ as a source term in the forward model and subtract from the measurement at working height, OR use the two-height data in a joint inversion.

**Deliverable:** $\sigma(x,y)$ background map and residual charge uncertainty $u(\sigma)$.

### 6.6 Step 5: Verification and uncertainty validation

**Procedure:**
1. Measure a **certified step-height standard** with known height $z_{\text{ref}} \pm U_{\text{ref}}$.
2. Compare the QESPM-measured $z_{\text{meas}} \pm U_{\text{meas}}$ with the reference value.
3. Check for consistency: $|z_{\text{meas}} - z_{\text{ref}}| \leq \sqrt{U_{\text{meas}}^2 + U_{\text{ref}}^2}$.
4. If consistent, the measurement capability is validated; if not, investigate and correct systematic errors.
5. Repeat verification at regular intervals (daily to weekly, depending on stability).

**Deliverable:** Validated measurement capability with documented uncertainty.

---

## 7. Comparison with ISO 25178 Surface Texture Parameters

### 7.1 ISO 25178 framework overview

ISO 25178 defines areal (3D) surface texture measurement in six parts:

| Part | Title | Relevance to QESPM |
|:----:|-------|-------------------|
| **-2** | Terms, definitions and surface texture parameters | Parameters computable from $z(x,y)$ |
| **-3** | Specification operators | Filters ($\lambda_s$, $\lambda_c$, $\lambda_f$) applicable |
| **-6** | Classification of methods | QESPM as new subclass |
| **-600** | Metrological characteristics | ITF, resolution, uncertainty — QESPM has distinct characteristics |
| **-601** | Nominal characteristics of contact (stylus) instruments | N/A (different physical principle) |
| **-602** | Nominal characteristics of non-contact (optical) instruments | **Closest analogue** — QESPM is non-contact but electrostatic, not optical |
| **-700** | Calibration and measurement uncertainty | Calibration procedure defined in §6 |

### 7.2 Height parameters (ISO 25178-2, §4)

| Parameter | Symbol | Definition | QESPM compatibility |
|-----------|--------|------------|---------------------|
| Arithmetic mean height | $S_a$ | $\frac{1}{A}\iint_A |z(x,y) - \bar{z}|\,dx\,dy$ | ✓ Fully compatible |
| Root mean square height | $S_q$ | $\sqrt{\frac{1}{A}\iint_A (z - \bar{z})^2\,dx\,dy}$ | ✓ Fully compatible |
| Skewness | $S_{sk}$ | $\frac{1}{S_q^3 A}\iint_A (z-\bar{z})^3\,dx\,dy$ | ✓ Requires charge correction |
| Kurtosis | $S_{ku}$ | $\frac{1}{S_q^4 A}\iint_A (z-\bar{z})^4\,dx\,dy$ | ✓ Requires charge correction |
| Maximum peak height | $S_p$ | $\max(z)$ | ✓ |
| Maximum pit height | $S_v$ | $\|\min(z)\|$ | ✓ |
| Maximum height | $S_z$ | $S_p + S_v$ | ✓ |

**Assessment**: All height parameters are computable from QESPM data. The limiting factor is the charge-induced bias — a spatially uniform charge error adds a DC offset that cancels in $S_a$ and $S_q$ (after mean subtraction) but biases $S_p$, $S_v$, and $S_z$.

### 7.3 Spatial parameters (ISO 25178-2, §5)

| Parameter | Symbol | Definition | QESPM compatibility |
|-----------|--------|------------|---------------------|
| Auto-correlation length | $S_{al}$ | Distance for ACF to decay to $0.2$ | ⚠ Limited by $\lambda_{\min} \approx 2h$ |
| Texture aspect ratio | $S_{tr}$ | Ratio of fastest to slowest ACF decay | ⚠ Requires isotropic ITF (dual-channel) |
| Texture direction | $S_{td}$ | Angle of dominant lay | ✓ (from $\phi_\mu$ channel) |

### 7.4 Hybrid parameters (ISO 25178-2, §6)

| Parameter | Symbol | QESPM notes |
|-----------|--------|-------------|
| RMS gradient | $S_{dq}$ | Directly measurable via micromotion amplitude $\beta \propto \|\nabla z\|$ — **unique capability** |
| Developed interfacial area ratio | $S_{dr}$ | Computable from $z(x,y)$ reconstruction |

### 7.5 Functional and feature parameters

Functional parameters ($S_{mr}$, $S_{mc}$, $S_{xp}$ — bearing area curve) and feature parameters ($S_{pd}$, $S_{pc}$, $S_{da}$ — peak/valley morphology) are computable from the reconstructed $z(x,y)$ but are subject to the lateral resolution limit. Peaks smaller than $\sim h/2$ are not resolved.

### 7.6 Filtration (ISO 25178-3)

The standard filters apply:

- **$\lambda_s$ filter** (short-wavelength cut-off): The ion probe's intrinsic ITF already provides low-pass filtering at $\lambda_{\min} \approx 2h$. An additional $\lambda_s$ filter may not be needed, but can be applied digitally to suppress noise at the Nyquist frequency.
- **$\lambda_c$ filter** (long-wavelength cut-off for roughness/waviness separation): Applies normally to separate form, waviness, and roughness.
- **$\lambda_f$ filter** (form removal): Polynomial form removal applies normally.

### 7.7 Key difference from conventional methods

The ion probe's ITF includes a **high-pass** component ($k^2$), meaning the instrument has **zero DC response** — it cannot measure absolute height, only height *variations*. This is fundamentally different from stylus instruments (which measure absolute displacement) and optical interferometers (which measure absolute phase). The ion probe is intrinsically a **curvature sensor** ($\Delta\omega \propto \partial^2 z/\partial x^2$), and height must be reconstructed by integration, with the integration constant lost.

---

## 8. New Class of Surface Measurement?

### 8.1 Criteria for a new measurement class

ISO 25178-6 defines a measurement method classification based on:
1. **Physical principle** of interaction
2. **Type of probe**
3. **Measured quantity**
4. **Measurement environment**

### 8.2 QESPM classification proposal

| Criterion | QESPM | Existing closest class |
|-----------|-------|----------------------|
| **Physical principle** | Electrostatic potential sensing via trapped-ion quantum metrology | Scanning Kelvin probe microscopy (SKPM), but SKPM measures work function, not topography |
| **Probe type** | Single trapped atomic ion (quantum object with well-defined motional states) | No existing class uses a quantum probe |
| **Measured quantity** | Secular frequency shift $\Delta\omega$, micromotion amplitude $\beta$ → surface height $z$ | No existing class |
| **Environment** | UHV ($<10^{-10}$ mbar), cryogenic ($4$ K optional) | STM (UHV, cryogenic); AFM (ambient to UHV) |
| **Interaction range** | Long-range ($h = 1$–$100$ µm) electrostatic | AFM: short-range contact/non-contact; STM: tunnelling ($\sim 1$ nm) |
| **Lateral resolution** | $h/2$ (2.5–50 µm typical) | AFM: $\sim 1$ nm; Optical: $\sim 0.5$ µm |
| **Vertical resolution** | Picometre to sub-nanometre | AFM: $\sim 0.1$ Å; Optical: $\sim 1$ nm |
| **ITF character** | Band-pass: $k^2 e^{-kh}$ | Optical: low-pass; Stylus: mechanical cut-off |

### 8.3 Proposed taxonomy

Under ISO 25178-6, Section 4 (Classification of areal surface texture measurement methods), we propose a new subclass:

> **4.X — Quantum electrostatic scanning probe microscopy (QESPM)**
>
> *Method in which a single trapped atomic ion is used as a local electrostatic potential probe. The ion is confined in a radio-frequency Paul trap at a controlled height above the surface. The surface is raster-scanned relative to the ion. The ion's motional observables — secular frequency shift, micromotion amplitude and phase, and motional heating rate — are measured at each lateral position and inverted to reconstruct the surface height distribution. The interaction is long-range electrostatic (not mechanical contact, not tunnelling), the probe is a quantum object with discrete motional states, and the instrument transfer function has a band-pass character $H(k) \propto k^2 e^{-kh}$.*

### 8.4 Unique metrological characteristics

What distinguishes QESPM from all existing ISO 25178 classes:

1. **Quantum probe**: The measurement sensitivity can be enhanced beyond the standard quantum limit using squeezed motional states — a capability unique among surface metrology instruments.

2. **Dual topography–materials sensing**: The heating rate channel $\dot{\bar{n}}(x,y)$ provides simultaneous materials contrast (resistivity, dielectric loss) that is quasi-independent of topography — no existing ISO 25178 method offers this.

3. **Self-calibrating frequency reference**: The ion's secular frequency serves as an intrinsic frequency standard, referenced to the SI second through the trap drive. No mechanical or optical calibration artefact is needed for frequency measurement.

4. **Non-contact at unprecedented range**: The probe operates at $1$–$100$ µm standoff — orders of magnitude further than AFM (nm contact) or STM (nm tunnelling) — yet achieves sub-nanometre vertical sensitivity.

5. **Vacuum/cryogenic native operation**: The ion requires UHV, making it uniquely suited for *in-situ* metrology of quantum device electrodes, superconducting circuits, and space-qualified surfaces — environments where conventional metrology instruments cannot operate.

### 8.5 Limitations that preclude it from replacing existing classes

1. **Lateral resolution** limited by ion height — cannot compete with AFM/STM for nanometre-scale lateral features.
2. **UHV requirement** restricts sample types.
3. **Conducting surface requirement** (or thin dielectric on conductor) — insulating samples accumulate charge and destabilise the trap.
4. **Throughput** is low ($\sim 10^{-2}$ µm²/s vs. $10^4$ µm²/s for optical profilometry).
5. **Complexity and cost** ($\sim 500\text{k}$–$1\text{M}$ USD vs. $\sim 50\text{k}$ USD for a research-grade AFM).

### 8.6 Recommendation

QESPM should be proposed as a **specialised complementary method** rather than a general-purpose replacement. Its niche is:
- Cryogenic UHV surface metrology for quantum devices
- Non-contact measurement of delicate or soft surfaces
- Combined topography and materials characterisation
- Quantum-enhanced precision surface metrology (long-term research goal)

A formal proposal to ISO/TC 213 (Geometrical product specifications) for a new subclass under ISO 25178-6 would require:
1. At least two independent experimental demonstrations
2. Inter-laboratory comparison data
3. Documented ITF and uncertainty budget
4. Calibration artefacts compatible with UHV and cryogenic operation

---

## 9. Conclusion

The single trapped-ion probe satisfies the formal requirements of a surface metrology instrument:

- ✅ **Defined measurand**: $z(x,y)$ — surface height profile
- ✅ **Measurement model**: Forward operator $\mathcal{F}$ with inverse solution
- ✅ **Instrument transfer function**: $H(k;h) \propto k^2 e^{-kh}$ — band-pass, height-dependent, calibratable
- ✅ **Spatial resolution**: Lateral $\sim h/2$, vertical $\sim$ picometre-to-nanometre
- ✅ **Uncertainty budget**: GUM-compliant, dominated by surface charge systematic
- ✅ **Traceability**: To SI metre (through calibrated grating), SI second (through GPS clock), SI volt (through Josephson standard)
- ✅ **Calibration procedure**: Five-step protocol with sinusoidal grating primary calibration
- ✅ **ISO 25178 compatibility**: Height, spatial, and hybrid parameters computable; filtration applicable
- 🆕 **New measurement class**: QESPM — justified by distinct physical principle (quantum electrostatic), unique ITF character (band-pass), and dual topography–materials capability

The **critical metrological gap** is the lack of a primary standard for surface charge distribution at micrometre scales — this is a fundamental limitation shared with Kelvin probe microscopy and all scanning electric-field methods. Addressing this gap (e.g., through in-situ surface preparation standards or charge-neutral reference surfaces) is a priority for establishing QESPM as a traceable measurement method.

---

*Document prepared as a metrology framework working draft. Intended for eventual submission to Metrologia or Measurement Science and Technology. All figures to be prepared as PDF vector graphics; manuscript in LaTeX.*
