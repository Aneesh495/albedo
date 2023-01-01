# ALBEDO: System Architecture & Subsystem Specification

ALBEDO is a high-performance hyperspectral imaging and multi-way tensor decomposition engine written in vectorized J (J9.6) with supporting Python orchestration utilities and a WebGL hardware-accelerated volumetric dashboard.

This document describes the architectural layout, memory models, tensor matricization semantics, data storage interleaved formats, and pipeline topologies across the system.

---

## 1. High-Level System Architecture

The ALBEDO architecture is decomposed into decoupled, modular subsystems where array operations are expressed as tacit forks, rank conjunctions, and high-dimensional contractions:

```
+---------------------------------------------------------------------------------------+
|                                ALBEDO RUNTIME ENGINE                                 |
+---------------------------------------------------------------------------------------+
                                           |
           +-------------------------------+-------------------------------+
           |                               |                               |
           v                               v                               v
+-----------------------+       +-----------------------+       +-----------------------+
|      IO & STORAGE     |       |     CORE NUMERICS     |       |    SPECTRAL DATA      |
|  - ENVI Header Parser |       |  - Jacobi SVD & QR    |       |  - USGS Minerals (15) |
|  - BSQ/BIL/BIP Layout |       |  - Cholesky Factorize |       |  - Botanical Taxa(10) |
|  - 3D Cube Container  |       |  - Mode-n Products    |       |  - Urban Targets (10) |
|  - Band Resampling    |       |  - Fast FFT & Wavelet |       |  - Trace Gases & HITRAN
+-----------------------+       +-----------------------+       +-----------------------+
           |                               |                               |
           +-------------------------------+-------------------------------+
                                           |
           +-------------------------------+-------------------------------+
           |                               |                               |
           v                               v                               v
+-----------------------+       +-----------------------+       +-----------------------+
|  ATMOSPHERIC TRANSFER |       |  SUBSPACE PROJECTION  |       |   SPECTRAL UNMIXING   |
|  - Solar Irradiance   |       |  - Spectral PCA       |       |  - VCA Simplex (p=4)  |
|  - TOA Reflectance    |       |  - Green's MNF        |       |  - N-FINDR Simplex    |
|  - DOS-1 / COST Path  |       |  - HySime Subspace Dim|       |  - FCLS Active Set QP |
|  - QUAC In-Scene Retr |       |  - FastICA Negentropy |       |  - SUNSAL ADMM L1-Reg |
|  - Continuum Removal  |       |  - Segmented PCA      |       |  - Endmember Bundles  |
+-----------------------+       +-----------------------+       +-----------------------+
                                           |
           +-------------------------------+-------------------------------+
           |                               |                               |
           v                               v                               v
+-----------------------+       +-----------------------+       +-----------------------+
|   TARGET DETECTION    |       |     TENSOR ALGEBRA    |       |  CLASSIFY & INDEX     |
|  - Global Reed-Xiaoli |       |  - Tucker3 HOSVD      |       |  - Spectral Angle Map |
|  - Dual-Window RX     |       |  - CP-ALS Decomposition|      |  - K-Means++ Cluster  |
|  - Nonlinear Kernel RX|       |  - RTPCA (ADMM / tSVD)|       |  - Graph Laplacian    |
|  - CEM & Matched Filt |       |  - Low-Rank Denoising |       |  - NDVI / NDWI / NBR  |
|  - AMSD Subspace Orth |       |  - Mode Folding/Unfold|       |  - Spatial Majority   |
+-----------------------+       +-----------------------+       +-----------------------+
```

---

## 2. Storage Layout & Memory Model

Hyperspectral cubes are 3-dimensional rank-3 arrays of dimensions `(B, H, W)`, representing bands $B$, rows $H$, and columns $W$.

### Interleave Permutations
Sensors write hyperspectral data in one of three standardized interleaves. ALBEDO implements lossless $O(1)$ stride and axis permutations:

1. **BSQ (Band Sequential)**: Shape `(B, H, W)`. Optimized for spatial 2D filtering and per-band visualization.
2. **BIL (Band Interleaved by Line)**: Shape `(H, B, W)`. Native line-scanner format matching pushbroom sensor readouts.
3. **BIP (Band Interleaved by Pixel)**: Shape `(H, W, B)`. Optimized for per-pixel spectrum extraction and point spectral libraries.

In J, these permutations are executed natively via dyadic transpose `|:` without manual nested looping:
- `bsq_to_bil =: 1 0 2&|:`
- `bil_to_bsq =: 1 0 2&|:`
- `bsq_to_bip =: 1 2 0&|:`
- `bip_to_bsq =: 2 0 1&|:`

```
      BSQ (B, H, W)              BIL (H, B, W)              BIP (H, W, B)
      +---------------+          +---------------+          +---------------+
     /               /|         /               /|         /               /|
    /   Band 0      / |        /    Line 0     / |        /    Pixel      / |
   +---------------+  |       +---------------+  |       +---------------+  |
   |               |  |  -->  | Band 0        |  |  -->  | Band 0 ... B  |  |
   |   Spatial     |  |       | Band 1        |  |       |               |  |
   |   (H x W)     | /        | ...           | /        |               | /
   +---------------+          +---------------+          +---------------+
```

---

## 3. Tensor Matricization and Unfolding Geometry

Multi-way data cubes $\mathcal{X} \in \mathbb{R}^{I_0 \times I_1 \times I_2}$ are matricized along mode $n$ into a 2D matrix $\mathbf{X}_{(n)} \in \mathbb{R}^{I_n \times \prod_{k \ne n} I_k}$.

### Mathematical Formulation
The mode-$n$ fibers of tensor $\mathcal{X}$ become the columns of $\mathbf{X}_{(n)}$. For a 3D cube $\mathcal{X} \in \mathbb{R}^{B \times H \times W}$:

1. **Mode-0 Unfolding (Spectral Matrix)**:
   $$\mathbf{X}_{(0)} \in \mathbb{R}^{B \times (HW)}$$
   Each column is a full pixel spectrum $\mathbf{r} \in \mathbb{R}^B$. This layout powers all linear unmixing (VCA, FCLS), subspace projection (PCA, MNF), and covariance estimation.

2. **Mode-1 Unfolding (Row Matricization)**:
   $$\mathbf{X}_{(1)} \in \mathbb{R}^{H \times (BW)}$$
   Permutation vector $(1, 0, 2)$ followed by 2D reshape.

3. **Mode-2 Unfolding (Column Matricization)**:
   $$\mathbf{X}_{(2)} \in \mathbb{R}^{W \times (BH)}$$
   Permutation vector $(2, 0, 1)$ followed by 2D reshape.

```
       3D Data Cube X                     Mode-0 Matricization X_(0)
      +---------------+                      1      2      ...    H*W
     /               /|                   +---------------------------+
    /               / |               1   | r_11   r_12   ...   r_HW  |
   +---------------+  |               2   | r_21   r_22   ...   r_HW  |
 B |               |  |  === unfold ===>  :   |  :      :            :    |
   |               | /                B   | r_B1   r_B2   ...   r_BHW |
   +---------------+                      +---------------------------+
         W      H                           (B bands  x  HW pixels)
```

In ALBEDO, unfolding and folding round-trip without copying overhead:
```j
tensor_unfold =: 4 : 0
  n =. x [ X =. y
  dims =. $ X [ order =. # dims
  other_modes =. (i. order) -. n
  perm =. n, other_modes
  permuted =. perm |: X
  dim_n =. n { dims [ other_prod =. */ (other_modes { dims)
  (dim_n, other_prod) $ , permuted
)
```

---

## 4. End-to-End Processing Pipelines

ALBEDO provides two production pipelines that automate processing from raw sensor digital numbers (DN) to geolocated chemical and mineral fraction maps:

### Pipeline 1: Automated Spectral Unmixing (`pipeline/pipeline_unmix.ijs`)
```
+---------------+     +---------------+     +---------------+
| Raw Cube (DN) | --> | Radiometric   | --> | Atmospheric   |
| (ENVI Header) |     | Calibration   |     | DOS-1 / QUAC  |
+---------------+     +---------------+     +---------------+
                                                    |
                                                    v
+---------------+     +---------------+     +---------------+
| Fully Constrain| <--| Simplex Vol   | <-- | Subspace Dim  |
| Least Squares |     | Maximization  |     | Estimation    |
| (FCLS Abund.) |     | (VCA/N-FINDR) |     | (HySime/MNF)  |
+---------------+     +---------------+     +---------------+
```

### Pipeline 2: Target & Anomaly Detection (`pipeline/pipeline_detect.ijs`)
```
+---------------+     +---------------+     +---------------+
| Surface       | --> | Spectral      | --> | Global RX /   |
| Reflectance   |     | Whitening     |     | Dual-Window   |
| Cube          |     | (Sigma^(-1/2))|     | Anomaly Map   |
+---------------+     +---------------+     +---------------+
                                                    |
                                                    v
+---------------+     +---------------+     +---------------+
| Classified    | <-- | Adaptive Cos  | <-- | Matched       |
| Anomaly Mask  |     | Estimator     |     | Filter (MF) / |
| (Top 1% Tail) |     | (ACE [0, 1])  |     | CEM Response  |
+---------------+     +---------------+     +---------------+
```

---

## 5. Algorithmic Complexity & Parallel Vectorization

All core computational kernels are vectorized over array ranks:

| Module | Operation | Mathematical Formula | Asymptotic Complexity |
|---|---|---|---|
| `core/linalg.ijs` | Jacobi Eigenvalues | $\mathbf{A} = \mathbf{V} \mathbf{\Lambda} \mathbf{V}^T$ | $O(N^3)$ per sweep |
| `core/linalg.ijs` | Householder QR | $\mathbf{A} = \mathbf{Q} \mathbf{R}$ | $O(2 M N^2 - \frac{2}{3} N^3)$ |
| `core/linalg.ijs` | Cholesky-Banachiewicz | $\mathbf{A} = \mathbf{L} \mathbf{L}^T$ | $O(\frac{1}{3} N^3)$ |
| `core/stats.ijs` | Shrinkage Covariance | $\mathbf{\Sigma}^* = (1 - \lambda)\mathbf{S} + \lambda \mathbf{T}$ | $O(B^2 N)$ |
| `subspace/mnf.ijs` | Maximum Noise Fraction | $\mathbf{\Sigma}_n^{-1/2} \mathbf{\Sigma}_x \mathbf{\Sigma}_n^{-1/2}$ | $O(B^3 + B^2 N)$ |
| `unmix/vca.ijs` | Vertex Component Analysis | Simplex vertex projection | $O(p^2 N)$ |
| `unmix/fcls.ijs` | Fully Constrained LS | Active-set QP (ANC + ASC) | $O(K \cdot p^3 N)$ |
| `detect/rx.ijs` | Global Reed-Xiaoli | $(\mathbf{r} - \boldsymbol{\mu})^T \mathbf{\Sigma}^{-1} (\mathbf{r} - \boldsymbol{\mu})$ | $O(B^2 N + B^3)$ |
| `tensor/tucker.ijs`| Tucker3 HOSVD | $\mathcal{X} \approx \mathcal{G} \times_1 \mathbf{U}_1 \times_2 \mathbf{U}_2 \times_3 \mathbf{U}_3$ | $O(I_1 I_2 I_3 \sum R_n)$ |
| `tensor/cp_als.ijs`| CP Decomposition | $\mathcal{X} \approx \sum_{r=1}^R \lambda_r \mathbf{a}_r \circ \mathbf{b}_r \circ \mathbf{c}_r$ | $O(R^3 + R I J K)$ |
