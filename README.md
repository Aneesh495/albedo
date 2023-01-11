# ALBEDO

[![Build Status](https://img.shields.io/badge/build-passing-brightgreen.svg)]()
[![Tests](https://img.shields.io/badge/tests-35%20passed%2C%200%20failed-success.svg)]()
[![J Version](https://img.shields.io/badge/language-J9.6%20vectorized-blue.svg)]()
[![License](https://img.shields.io/badge/license-MIT-green.svg)]()

**ALBEDO** is a high-dimensional hyperspectral imaging, satellite spectroscopy, and multi-way tensor decomposition engine written in vectorized J (J9.6) with Python orchestration glue and an interactive WebGL hardware-accelerated volumetric dashboard.

Built from first principles, ALBEDO solves real-world Earth Observation challenges: atmospheric radiative transfer normalization, continuum removal absorption spectroscopy, unsupervised subspace dimensionality estimation, simplex volume maximization, active-set convex quadratic programming for physical abundance fractions, subpixel gas plume detection (e.g. greenhouse methane CH4 emissions), multi-way Tucker3/CP tensor factorizations, and multiscale wavelet denoising.

```
                  _   _              _         
        /\       | | | |            | |        
       /  \      | | | |__     ___  | |   ___  
      / /\ \     | | | '_ \   / _ \ | |  / _ \ 
     / ____ \    | | | |_) | |  __/ | | | (_) |
    /_/    \_\   |_| |_.__/   \___| |_|  \___/ 
                                               
    Hyperspectral Spectroscopy & Tensor Decomposition Engine
```

---

## Architecture Overview

```
+---------------------------------------------------------------------------------------+
|                                  ALBEDO ENGINE CORE                                   |
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
           |                               |                               |
           +-------------------------------+-------------------------------+
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

## Key Modules & Mathematical Capabilities

### 1. Atmospheric Radiative Transfer & Normalization (`atmos/`)
- **Thuillier Solar Irradiance Model**: $E_0(\lambda)$ calibrated extraterrestrial solar spectrum.
- **Top-of-Atmosphere (TOA) Reflectance**: Earth-Sun distance factor $d^2(t)$, solar declination, and solar zenith angle corrections.
- **Dark Object Subtraction (DOS-1 & COST)**: Path radiance $L_p(\lambda)$ estimation and optical thickness Rayleigh scattering removal.
- **Quick Atmospheric Correction (QUAC)**: Endmember-guided scene statistics reflectance retrieval without requiring aerosol optical depth Look-Up Tables.
- **Continuum Removal**: Upper convex hull normalization, diagnostic absorption depth, band center, integrated absorption area, and Full Width at Half Maximum (FWHM).

### 2. Subspace Dimensionality Reduction (`subspace/`)
- **Spectral Principal Component Analysis (PCA)**: Covariance eigen-decomposition with bidirectional forward projection and inverse subspace reconstruction.
- **Green's Maximum Noise Fraction (MNF)**: Signal-to-Noise Ratio maximization using spatial difference noise covariance estimation and dual-stage whitening.
- **Hyperspectral Subspace Identification by Minimum Error (HySime)**: Unsupervised, parameter-free estimation of intrinsic signal subspace dimension.
- **FastICA**: Fixed-point negentropy maximization with symmetric deflationary Gram-Schmidt orthogonalization.
- **Segmented PCA**: Window-partitioned subspace projections across VNIR, SWIR-1, and SWIR-2 atmospheric transmission windows.

### 3. Spectral Unmixing & Optimization (`unmix/`)
- **Vertex Component Analysis (VCA)**: Simplex vertex identification via successive random orthogonal subspace projections.
- **N-FINDR**: Simplex volume maximization via determinant evaluations.
- **Pixel Purity Index (PPI)**: Random skewer projection histograms on extreme convex vertices.
- **Non-Negative Least Squares (NNLS)**: Lawson-Hanson active-set quadratic programming for physical Abundance Non-negativity Constraints (ANC).
- **Fully Constrained Least Squares (FCLS)**: Active-set dual Lagrangian solver satisfying both ANC ($a_j \ge 0$) and Abundance Sum-to-One Constraints (ASC, $\sum a_j = 1$).
- **Sparse Unmixing via Variable Splitting (SUNSAL)**: Alternating Direction Method of Multipliers (ADMM) $\ell_1$-regularized unmixing over extensive spectral libraries.
- **Endmember Bundles & Extended Linear Mixing Model (ELMM)**: Per-pixel spectral variability modeling accommodating natural illumination and topographic shading.
- **Non-Negative Matrix Factorization (NMF)**: Lee-Seung multiplicative updates for joint blind endmember and abundance matrix recovery.

### 4. Target & Anomaly Detection Kernels (`detect/`)
- **Global Reed-Xiaoli (RX)**: Multivariate Gaussian benchmark anomaly detector evaluating quadratic Mahalanobis distances against scene background covariance.
- **Dual-Window Local RX**: Annular background estimation sliding kernel with a central guard window isolating target spatial signatures.
- **Nonlinear Kernel RX**: Gaussian RBF reproducing kernel Hilbert space (RKHS) projection for nonlinearly distributed spectral clutter.
- **Constrained Energy Minimization (CEM)**: Finite impulse response filter minimizing output background variance while maintaining unit gain on target signatures.
- **Spectral Matched Filter (MF)**: Mean-subtracted covariance-whitened target detector.
- **Adaptive Cosine Estimator (ACE)**: Scale-invariant and illumination-invariant angle cosine detector normalized in $[0, 1]$.
- **Adaptive Matched Subspace Detector (AMSD)**: Target and background subspace orthogonal complement projections.

### 5. Multi-Way Tensor Algebra & Decomposition (`core/tensor.ijs`, `tensor/`)
- **Tensor Matricization & Folding**: Efficient mode-$n$ unfolding $\mathbf{X}_{(n)}$ and rank-preserving tensor folding.
- **Mode-$n$ Matrix-Tensor Products**: Multi-mode contractions ($\mathcal{Y} = \mathcal{X} \times_1 \mathbf{U}_1 \times_2 \mathbf{U}_2 \times_3 \mathbf{U}_3$).
- **Kronecker & Khatri-Rao Products**: Vectorized column-wise Kronecker products powering multi-way decompositions.
- **Tucker3 / Higher-Order SVD (HOSVD)**: Low-rank core tensor truncation and orthonormal factor matrix estimation.
- **Canonical Polyadic (CP) Decomposition**: Alternating Least Squares (CP-ALS) rank-$R$ outer-product factor matrix recovery.
- **Robust Tensor PCA (RTPCA)**: ADMM decomposition into low-rank background tensor $\mathcal{L}$ and sparse target anomaly tensor $\mathcal{S}$ via tensor singular value thresholding (t-SVT).

### 6. Signal Transforms & Spatial Operators (`signal/`)
- **Discrete Wavelet Transform (DWT)**: Multilevel 1D and 2D Mallat pyramidal filterbanks with Haar and Daubechies (db4) orthogonal wavelets.
- **VisuShrink Wavelet Denoising**: Universal soft-thresholding for spectral signal recovery.
- **Grünwald-Letnikov Fractional Calculus Differentiation**: Arbitrary non-integer order differentiation ($\alpha \in \mathbb{R}^+$) via Gamma recurrence coefficients, enhancing weak absorption shoulders.
- **Spatial Filters**: 2D convolution, Gaussian blur, Sobel edge gradients, Laplacian sharpening, and bilateral edge-preserving smoothing.

### 7. Calibrated Reference Libraries (`data/`)
- **USGS Minerals**: 15 calibrated laboratory mineral spectra (Kaolinite, Alunite, Calcite, Dolomite, Hematite, Goethite, Montmorillonite, Illite, Muscovite, Chlorite, Epidote, Gypsum, Quartz, Jarosite).
- **Botanical Taxa**: 10 vegetation species with red-edge dynamics and chlorophyll fluorescence (broadleaf, conifer, maize, soybean, vineyard, sedge, moss, algae, lichen, turf).
- **Urban & Anthropogenic**: 10 infrastructural materials (fresh/aged asphalt, concrete, galvanized steel, terracotta, EPDM, PVC, red brick, sand).
- **ECOSTRESS Lithology & Soils**: 30 rock types (igneous, metamorphic, sedimentary) and 25 soil orders.
- **Polymers & Marine Debris**: 20 polymer spectral profiles with diagnostic C-H, N-H, O-H, and C=O overtone doublets (HDPE, LDPE, PP, PET, PVC, PS, PMMA, PTFE).
- **HITRAN Trace Gases**: 8 atmospheric trace gas absorption cross-sections (CH4, CO2, CO, SO2, NH3, NO2, N2O, H2O) with Beer-Lambert plume attenuation.

---

## Performance Benchmarks

ALBEDO achieves execution latencies measured on an Apple M-series workstation:

```
================================================================================
                    ALBEDO MICROBENCHMARK AUDIT TABLE
================================================================================
| #   | Benchmark Kernel                       | Complexity         | Latency (ms) | Implementation               |
|-----|----------------------------------------|--------------------|--------------|------------------------------|
| 1   | Thin SVD (50x20)                       | O(m n^2)           | 10           | Dense Jacobi SVD             |
| 2   | Cholesky Factorization (100x100)       | O(n^3 / 3)         | 3            | Cholesky-Banachiewicz        |
| 3   | Covariance Estimation (50x500)         | O(p^2 N)           | < 1          | Two-pass shrinkage           |
| 4   | VCA Simplex Extraction (p=4, N=1000)   | O(p^2 N)           | 38           | Vertex Component Analysis    |
| 5   | Reed-Xiaoli Detection (20x25x20)       | O(p^2 N)           | < 1          | Quadratic Mahalanobis metric |
| 6   | Tensor Mode-0 Contraction (20x20x20)   | O(I J K J')        | < 1          | Matricized tensor product    |
| 7   | 4-Level Wavelet DWT (4096 samples)     | O(N)               | < 1          | Mallat pyramidal filterbank  |
| 8   | Fractional Derivative (alpha=0.5, 224) | O(B^2)             | < 1          | Grunwald-Letnikov difference |
| 9   | Tucker3 HOSVD (15x15x15, rank 4)       | O(I1 I2 I3 R)      | 12           | Higher-Order SVD             |
| 10  | CP-ALS Factorization (10x10x10, rank 3) | O(R^3 + R I J K)   | < 1          | Alternating Least Squares    |
================================================================================
Total benchmarks executed: 10
Status: 100% BENCHMARKS PASSED
```

---

## Verification & Test Suite

The engine is backed by 35 automated test suites with 100% pass rate:

```bash
$ make test
==> Running full ALBEDO test suite (35 suites)...
================================================================================
                     ALBEDO AUTOMATED TEST SUITE RUNNER                         
                     J Interpreter: /Users/aneeshkrishna/.local/bin/jconsole
                     Discovered 35 Test Suites
================================================================================
--> Running tests/test_amsd.ijs                           ... [PASS] (0.021s)
--> Running tests/test_aquatic.ijs                        ... [PASS] (0.020s)
--> Running tests/test_atmos.ijs                          ... [PASS] (0.021s)
--> Running tests/test_atmos_lut.ijs                      ... [PASS] (0.021s)
--> Running tests/test_bundles.ijs                        ... [PASS] (0.024s)
--> Running tests/test_detect.ijs                         ... [PASS] (0.026s)
--> Running tests/test_envi_io.ijs                        ... [PASS] (0.021s)
--> Running tests/test_fault_reproduction.ijs             ... [PASS] (0.024s)
--> Running tests/test_fft.ijs                            ... [PASS] (0.020s)
--> Running tests/test_fractional.ijs                     ... [PASS] (0.020s)
--> Running tests/test_gases.ijs                          ... [PASS] (0.020s)
--> Running tests/test_ica.ijs                            ... [PASS] (0.021s)
--> Running tests/test_indices_classify.ijs               ... [PASS] (0.023s)
--> Running tests/test_krx.ijs                            ... [PASS] (0.025s)
--> Running tests/test_linalg.ijs                         ... [PASS] (0.021s)
--> Running tests/test_linalg_ext.ijs                     ... [PASS] (0.021s)
--> Running tests/test_morphology.ijs                     ... [PASS] (0.020s)
--> Running tests/test_multiway_tensor.ijs                ... [PASS] (0.023s)
--> Running tests/test_nmf.ijs                            ... [PASS] (0.021s)
--> Running tests/test_pipelines.ijs                      ... [PASS] (0.040s)
--> Running tests/test_plastics.ijs                       ... [PASS] (0.021s)
--> Running tests/test_quac.ijs                           ... [PASS] (0.035s)
--> Running tests/test_quadrature.ijs                     ... [PASS] (0.020s)
--> Running tests/test_rocks.ijs                          ... [PASS] (0.020s)
--> Running tests/test_rtpca.ijs                          ... [PASS] (0.021s)
--> Running tests/test_soils.ijs                          ... [PASS] (0.022s)
--> Running tests/test_spatial_filter.ijs                 ... [PASS] (0.022s)
--> Running tests/test_spectral_graph.ijs                 ... [PASS] (0.032s)
--> Running tests/test_spectral_libraries.ijs             ... [PASS] (0.021s)
--> Running tests/test_spline.ijs                         ... [PASS] (0.021s)
--> Running tests/test_stats.ijs                          ... [PASS] (0.021s)
--> Running tests/test_subspace.ijs                       ... [PASS] (0.023s)
--> Running tests/test_tensor.ijs                         ... [PASS] (0.022s)
--> Running tests/test_unmix.ijs                          ... [PASS] (0.029s)
--> Running tests/test_wavelet.ijs                        ... [PASS] (0.020s)
================================================================================
RESULTS: 35 PASSED, 0 FAILED out of 35 suites in 0.80s
================================================================================
```

---

## Quickstart

### Prerequisites
- J Interpreter: [Jsoftware](https://code.jsoftware.com/wiki/System/Installation) (`jconsole`, J9.5 or J9.6).
- Python 3.8+ (for CLI driver and WebGL dashboard).

### Clone and Verify
```bash
git clone https://github.com/aneesh495/albedo.git
cd albedo

# Run all 35 test suites
make test

# Run microbenchmarks
make bench

# Verify zero em-dash compliance
make check
```

### Command-Line Interface (`bin/albedo`)
```bash
# 1. Simulate a synthetic 32x32 hyperspectral scene (224 bands)
./bin/albedo sim --rows 32 --cols 32

# 2. Extract 4 pure endmembers via VCA and compute FCLS abundance fractions
./bin/albedo unmix --endmembers 4 --method vca

# 3. Detect anomalies via Global Reed-Xiaoli (RX)
./bin/albedo detect --method rx

# 4. Detect target signatures via Constrained Energy Minimization (CEM)
./bin/albedo detect --method cem

# 5. Launch interactive WebGL hyperspectral visualizer
./bin/albedo viewer --port 8000
```

### Interactive J Session
```j
$ jconsole albedo.ijs
========================================================================
        ALBEDO: Hyperspectral Imaging & Tensor Decomposition Engine    
        Version: 1.0.0 | Pure J Vectorized Numerics      
========================================================================
Modules initialized: core, io, atmos, subspace, unmix, detect, tensor,
                     indices, classify, signal, data, sim, pipeline.

   NB. Synthesize 24x24 scene (224 bands, SNR = 30 dB)
   scene =: albedo_sim_generate_scene 24 24 30.0
   cube =: 0 {:: scene

   NB. Run Global Reed-Xiaoli Anomaly Detector
   scores =: albedo_rx_detect_global_rx cube
   echo 'Anomaly map shape: ', ": $ scores
Anomaly map shape: 24 24

   NB. Extract 4 endmembers via VCA
   E =: 4 albedo_vca_extract_vca cube
   echo 'Endmember matrix shape: ', ": $ E
Endmember matrix shape: 224 4
```

---

## Documentation Links

- [System Architecture & Memory Models](docs/architecture.md)
- [Mathematical Foundations & Derivations](docs/mathematics.md)
- [Tacit Programming Guide in J](docs/tacit_programming_guide.md)

---

## License

MIT License. Developed with passion by Aneesh Krishna Parthasarathy.
