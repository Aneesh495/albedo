# ALBEDO: Mathematical Foundations & Algorithmic Derivations

This document provides rigorous mathematical derivations for the algorithms implemented across the ALBEDO engine.

---

## 1. Atmospheric Radiative Transfer & Solar Normalization

### Top-of-Atmosphere (TOA) Spectral Reflectance
Satellite sensors measure calibrated spectral radiance $L_{at-sensor}(\lambda)$ in units of $\mathrm{W} / (\mathrm{m}^2 \cdot \mathrm{sr} \cdot \mathrm{nm})$. The apparent Top-of-Atmosphere (TOA) planetary reflectance $\rho_{TOA}(\lambda)$ is given by:

$$\rho_{TOA}(\lambda) = \frac{\pi \cdot L_{at-sensor}(\lambda) \cdot d^2(t)}{E_0(\lambda) \cdot \cos(\theta_z)}$$

Where:
- $E_0(\lambda)$ is the extraterrestrial solar spectral irradiance at mean Earth-Sun distance (Thuillier 2003 solar reference model).
- $\theta_z$ is the solar zenith angle, related to solar elevation $\beta$ by $\theta_z = 90^\circ - \beta$.
- $d(t)$ is the Earth-Sun astronomical distance correction factor as a function of Day of Year (DOY) $t \in [1, 365]$:

$$d(t) = 1.0 - 0.01672 \cdot \cos\left(\frac{2\pi (t - 4)}{365.25}\right)$$

$$d^2(t) \approx 1.000110 + 0.034221 \cos(\Gamma) + 0.001280 \sin(\Gamma) + 0.000719 \cos(2\Gamma) + 0.000077 \sin(2\Gamma)$$

where $\Gamma = \frac{2\pi(t - 1)}{365}$ is the fractional day angle in radians.

### Dark Object Subtraction (DOS-1 & COST Models)
Surface reflectance $\rho_{surf}(\lambda)$ requires eliminating path radiance $L_p(\lambda)$ caused by atmospheric Rayleigh scattering and aerosol backscatter:

$$L_{surf}(\lambda) = L_{at-sensor}(\lambda) - L_p(\lambda)$$

$$\rho_{surf}(\lambda) = \frac{\pi \cdot (L_{at-sensor}(\lambda) - L_p(\lambda)) \cdot d^2}{E_0(\lambda) \cdot \cos(\theta_z) \cdot T_v(\lambda) \cdot T_z(\lambda)}$$

In the Chavez (1996) COST model, atmospheric transmittance along the solar path $T_z(\lambda)$ and view path $T_v(\lambda)$ is approximated by $T_z(\lambda) \approx \cos(\theta_z)$ and $T_v(\lambda) \approx 1.0$, yielding:

$$\rho_{surf}(\lambda) = \frac{\pi \cdot (L_{at-sensor}(\lambda) - L_p(\lambda)) \cdot d^2}{E_0(\lambda) \cdot \cos^2(\theta_z)}$$

---

## 2. Continuum Removal & Absorption Feature Analysis

Continuum removal normalizes spectral reflectance by fitting a convex hull over the local spectrum, isolating diagnostic absorption bands:

$$\rho_{cr}(\lambda) = \frac{\rho(\lambda)}{C(\lambda)}$$

Where $C(\lambda)$ is the upper convex hull formed by linear tie lines connecting local reflectance maxima.

### Feature Metrics
1. **Absorption Band Depth**:
   $$D(\lambda) = 1.0 - \rho_{cr}(\lambda)$$
   Peak absorption depth: $D_{max} = \max_\lambda D(\lambda) = 1.0 - \min_\lambda \rho_{cr}(\lambda)$.

2. **Absorption Center**:
   $$\lambda_{center} = \arg\max_\lambda D(\lambda)$$

3. **Full Width at Half Maximum (FWHM)**:
   $$\mathrm{FWHM} = \lambda_{right} - \lambda_{left}$$
   where $D(\lambda_{left}) = D(\lambda_{right}) = \frac{1}{2} D_{max}$.

4. **Integrated Absorption Area**:
   $$\mathrm{Area} = \int_{\lambda_{start}}^{\lambda_{stop}} (1.0 - \rho_{cr}(\lambda)) \, d\lambda \approx \sum_{i} \frac{D(\lambda_i) + D(\lambda_{i+1})}{2} (\lambda_{i+1} - \lambda_i)$$

---

## 3. Subspace Projection & Noise Whitening

### Green's Maximum Noise Fraction (MNF)
The MNF transform (Green et al., 1988) seeks linear combinations $\mathbf{y} = \mathbf{w}^T \mathbf{x}$ that maximize the Signal-to-Noise Ratio (SNR):

$$\mathrm{SNR} = \frac{\mathbf{w}^T \mathbf{\Sigma}_s \mathbf{w}}{\mathbf{w}^T \mathbf{\Sigma}_n \mathbf{w}} = \frac{\mathbf{w}^T \mathbf{\Sigma}_x \mathbf{w}}{\mathbf{w}^T \mathbf{\Sigma}_n \mathbf{w}} - 1$$

Where:
- $\mathbf{\Sigma}_x$ is the total scene covariance matrix.
- $\mathbf{\Sigma}_n$ is the noise covariance matrix, estimated via spatial difference operators between adjacent horizontal and vertical pixels:

$$\mathbf{n}(h, w) = \frac{1}{2} \left[ (\mathbf{x}(h, w) - \mathbf{x}(h, w+1)) + (\mathbf{x}(h, w) - \mathbf{x}(h+1, w)) \right]$$

The generalized eigenvalue problem $\mathbf{\Sigma}_x \mathbf{w} = \lambda \mathbf{\Sigma}_n \mathbf{w}$ is solved via a two-step symmetric eigen-decomposition:
1. Whiten noise covariance: $\mathbf{E}_n^T \mathbf{\Sigma}_n \mathbf{E}_n = \mathbf{\Lambda}_n$, giving whitening transform $\mathbf{F} = \mathbf{E}_n \mathbf{\Lambda}_n^{-1/2}$.
2. Transform total covariance: $\mathbf{\Sigma}_{adj} = \mathbf{F}^T \mathbf{\Sigma}_x \mathbf{F}$.
3. Diagonalize adjusted covariance: $\mathbf{G}^T \mathbf{\Sigma}_{adj} \mathbf{G} = \mathbf{\Lambda}_{mnf}$.
4. Combined MNF projector: $\mathbf{W}_{MNF} = \mathbf{F} \mathbf{G}$.

### HySime: Hyperspectral Subspace Identification
HySime (Bioucas-Dias & Nascimento, 2008) estimates the intrinsic subspace dimensionality $k$ by minimizing the sum of projection error and noise power:

$$k^* = \arg\min_k \left[ \sum_{i=1}^k (\sigma_i^2 - 2 \hat{\sigma}_{n, i}^2) \right]$$

---

## 4. Spectral Unmixing: Endmember Extraction & Abundance Inversion

The Linear Mixing Model (LMM) assumes observed pixel spectrum $\mathbf{r} \in \mathbb{R}^B$ is a convex combination of $p$ pure spectral endmembers $\mathbf{M} \in \mathbb{R}^{B \times p}$:

$$\mathbf{r} = \mathbf{M} \mathbf{a} + \mathbf{n} = \sum_{j=1}^p a_j \mathbf{m}_j + \mathbf{n}$$

Subject to physical constraints:
1. **Abundance Non-negativity Constraint (ANC)**: $a_j \ge 0, \quad \forall j \in \{1, \dots, p\}$.
2. **Abundance Sum-to-One Constraint (ASC)**: $\sum_{j=1}^p a_j = 1$.

### Vertex Component Analysis (VCA)
VCA (Nascimento & Bioucas-Dias, 2005) extracts the vertices of the endmember simplex by iteratively projecting data onto random orthogonal directions:
1. Project data onto affine subspace of dimension $p-1$.
2. For each iteration $i = 1, \dots, p$:
   - Generate random unit vector $\mathbf{w} \sim \mathcal{N}(0, \mathbf{I})$.
   - Orthogonalize $\mathbf{w}$ against previously identified endmembers:
     $$\mathbf{f} = \left(\mathbf{I} - \mathbf{V}_{i-1} \mathbf{V}_{i-1}^\dagger\right) \mathbf{w}$$
   - Identify vertex maximizing absolute inner product:
     $$j^* = \arg\max_{j} \left| \mathbf{f}^T \mathbf{x}_j \right|$$
   - Append endmember $\mathbf{m}_i = \mathbf{r}_{j^*}$.

### Fully Constrained Least Squares (FCLS)
FCLS solves the quadratic programming problem for each pixel:

$$\min_{\mathbf{a}} \frac{1}{2} \|\mathbf{r} - \mathbf{M} \mathbf{a}\|_2^2 \quad \text{s.t.} \quad \mathbf{a} \ge \mathbf{0}, \quad \mathbf{1}^T \mathbf{a} = 1$$

We formulate the augmented Lagrangian with active-set constraint enforcement:
$$\mathbf{M}_{aug} = \begin{bmatrix} \delta \mathbf{M} \\ \mathbf{1}^T \end{bmatrix}, \quad \mathbf{r}_{aug} = \begin{bmatrix} \delta \mathbf{r} \\ 1 \end{bmatrix}$$
where $\delta > 0$ balances the residual fit with the sum-to-one constraint. The Lawson-Hanson active set method iterates between passive (free) set $\mathcal{P}$ and active (clamped to zero) set $\mathcal{Z}$:
$$\mathbf{a}_{\mathcal{P}} = (\mathbf{M}_{\mathcal{P}}^T \mathbf{M}_{\mathcal{P}})^{-1} \mathbf{M}_{\mathcal{P}}^T \mathbf{r}$$
Negative components in $\mathbf{a}_{\mathcal{P}}$ are projected back to zero along the active search direction.

### Sparse Unmixing (SUNSAL) via ADMM
SUNSAL (Bioucas-Dias & Figueiredo, 2010) finds sparse endmember representations over large spectral libraries ($\mathbf{M} \in \mathbb{R}^{B \times L}$ with $L \gg B$):

$$\min_{\mathbf{A}} \frac{1}{2} \|\mathbf{Y} - \mathbf{M} \mathbf{A}\|_F^2 + \lambda \|\mathbf{A}\|_1 + \mathbb{I}_{\mathbb{R}_+}(\mathbf{A})$$

Using the Alternating Direction of Method of Multipliers (ADMM) with variable splitting $\mathbf{V} = \mathbf{A}$:
$$\mathbf{A}^{(k+1)} = (\mathbf{M}^T \mathbf{M} + \mu \mathbf{I})^{-1} (\mathbf{M}^T \mathbf{Y} + \mu (\mathbf{V}^{(k)} - \mathbf{D}^{(k)}))$$
$$\mathbf{V}^{(k+1)} = \max\left(0, \mathcal{S}_{\lambda/\mu}(\mathbf{A}^{(k+1)} + \mathbf{D}^{(k)})\right)$$
$$\mathbf{D}^{(k+1)} = \mathbf{D}^{(k)} + \mathbf{A}^{(k+1)} - \mathbf{V}^{(k+1)}$$
where $\mathcal{S}_\tau(x) = \mathrm{sign}(x) \max(0, |x| - \tau)$ is the soft-thresholding operator.

---

## 5. Target and Anomaly Detection

Hyperspectral detection is formulated as binary hypothesis testing for pixel $\mathbf{r} \in \mathbb{R}^B$:
- $H_0$: Target absent (background only), $\mathbf{r} \sim \mathcal{N}(\boldsymbol{\mu}_b, \mathbf{\Sigma}_b)$
- $H_1$: Target present, $\mathbf{r} \sim \mathcal{N}(\alpha \mathbf{d} + \boldsymbol{\mu}_b, \mathbf{\Sigma}_b)$

### Reed-Xiaoli (RX) Anomaly Detector
The Generalized Likelihood Ratio Test (GLRT) for unspecified anomalous targets with unknown covariance yields the quadratic Mahalanobis metric:

$$D_{RX}(\mathbf{r}) = (\mathbf{r} - \boldsymbol{\mu}_b)^T \mathbf{\Sigma}_b^{-1} (\mathbf{r} - \boldsymbol{\mu}_b)$$

In Dual-Window Local RX, $\boldsymbol{\mu}_b$ and $\mathbf{\Sigma}_b$ are computed dynamically over a local sliding annular background window, excluding a central guard window:
- Inner Guard Window $W_{in} \times W_{in}$: Isolates target spatial signature.
- Outer Background Window $W_{out} \times W_{out}$: Samples surrounding local clutter.

### Constrained Energy Minimization (CEM)
CEM designs a finite impulse response (FIR) filter $\mathbf{w} \in \mathbb{R}^B$ that minimizes output energy while passing desired target signature $\mathbf{d}$ with unit gain:

$$\min_{\mathbf{w}} \mathbf{w}^T \mathbf{R} \mathbf{w} \quad \text{s.t.} \quad \mathbf{w}^T \mathbf{d} = 1$$

Where $\mathbf{R} = \frac{1}{N} \mathbf{X} \mathbf{X}^T$ is the sample correlation matrix. Solving via Lagrange multipliers:

$$\mathbf{w}_{CEM} = \frac{\mathbf{R}^{-1} \mathbf{d}}{\mathbf{d}^T \mathbf{R}^{-1} \mathbf{d}}, \quad y_{CEM}(\mathbf{r}) = \mathbf{w}_{CEM}^T \mathbf{r} = \frac{\mathbf{d}^T \mathbf{R}^{-1} \mathbf{r}}{\mathbf{d}^T \mathbf{R}^{-1} \mathbf{d}}$$

### Adaptive Cosine Estimator (ACE)
ACE provides scale and illumination invariance by normalizing with respect to both target and pixel Mahalanobis lengths:

$$D_{ACE}(\mathbf{r}) = \frac{(\mathbf{d}^T \mathbf{\Sigma}^{-1} (\mathbf{r} - \boldsymbol{\mu}))^2}{(\mathbf{d}^T \mathbf{\Sigma}^{-1} \mathbf{d}) \cdot ((\mathbf{r} - \boldsymbol{\mu})^T \mathbf{\Sigma}^{-1} (\mathbf{r} - \boldsymbol{\mu}))} \in [0, 1]$$

---

## 6. Multi-Way Tensor Decompositions

A 3D hyperspectral cube is modeled as a 3rd-order tensor $\mathcal{X} \in \mathbb{R}^{I_1 \times I_2 \times I_3}$.

### Tucker3 / Higher-Order SVD (HOSVD)
Tucker decomposition decomposes $\mathcal{X}$ into a core tensor $\mathcal{G} \in \mathbb{R}^{R_1 \times R_2 \times R_3}$ multiplied by factor matrices along each mode:

$$\mathcal{X} \approx \mathcal{G} \times_1 \mathbf{U}_1 \times_2 \mathbf{U}_2 \times_3 \mathbf{U}_3$$

HOSVD computes each orthogonal factor matrix $\mathbf{U}_n$ as the leading $R_n$ left singular vectors of mode-$n$ unfolding $\mathbf{X}_{(n)}$:
$$\mathbf{X}_{(n)} = \mathbf{U}_n \mathbf{\Sigma}_n \mathbf{V}_n^T$$
The core tensor is retrieved via multi-mode contraction:
$$\mathcal{G} = \mathcal{X} \times_1 \mathbf{U}_1^T \times_2 \mathbf{U}_2^T \times_3 \mathbf{U}_3^T$$

### Canonical Polyadic (CP) Decomposition via ALS
CP-ALS approximates $\mathcal{X}$ as the sum of $R$ rank-1 tensors:

$$\mathcal{X} \approx \sum_{r=1}^R \lambda_r (\mathbf{a}_r \circ \mathbf{b}_r \circ \mathbf{c}_r)$$

The Alternating Least Squares (ALS) update fixes two factor matrices and solves the unconstrained least-squares subproblem for the third:
$$\mathbf{A} \leftarrow \mathbf{X}_{(0)} (\mathbf{C} \odot \mathbf{B}) \left[(\mathbf{C}^T \mathbf{C}) * (\mathbf{B}^T \mathbf{B})\right]^{-1}$$
$$\mathbf{B} \leftarrow \mathbf{X}_{(1)} (\mathbf{C} \odot \mathbf{A}) \left[(\mathbf{C}^T \mathbf{C}) * (\mathbf{A}^T \mathbf{A})\right]^{-1}$$
$$\mathbf{C} \leftarrow \mathbf{X}_{(2)} (\mathbf{B} \odot \mathbf{A}) \left[(\mathbf{B}^T \mathbf{B}) * (\mathbf{A}^T \mathbf{A})\right]^{-1}$$
where $\odot$ denotes the Khatri-Rao product and $*$ is the Hadamard (element-wise) product.

---

## 7. Grünwald-Letnikov Fractional Calculus Differencing

Fractional-order differentiation enhances weak absorption features obscured by broad continuum slopes. The Grünwald-Letnikov fractional derivative of order $\alpha \in \mathbb{R}^+$ for discrete spectrum $f[n]$ is defined as:

$$\Delta^\alpha f[n] = \sum_{j=0}^{K} (-1)^j \binom{\alpha}{j} f[n - j] = \sum_{j=0}^K w_j^{(\alpha)} f[n - j]$$

Where the filter weights $w_j^{(\alpha)}$ are computed via the Gamma recurrence:

$$w_0^{(\alpha)} = 1, \quad w_j^{(\alpha)} = \left(1 - \frac{\alpha + 1}{j}\right) w_{j-1}^{(\alpha)}, \quad j = 1, 2, \dots, K$$

For intermediate orders such as $\alpha = 0.5$ (half-order derivative), the transform simultaneously removes linear baseline trends while preserving peak locations and asymmetric shoulder geometries.
