NB. ==============================================================================
NB. ALBEDO - High-Dimensional Hyperspectral Imaging Engine
NB. scripts/benchmarks.ijs : Throughput and Latency Microbenchmarks
NB.
NB. Evaluates performance across dense linear algebra, tensor decompositions,
NB. spectral unmixing optimization, signal transforms, and anomaly detection.
NB. ==============================================================================

cocurrent 'base'

load 'core/linalg.ijs'
load 'core/linalg_ext.ijs'
load 'core/stats.ijs'
load 'core/tensor.ijs'
load 'core/fft.ijs'
load 'subspace/pca.ijs'
load 'subspace/mnf.ijs'
load 'unmix/vca.ijs'
load 'unmix/fcls.ijs'
load 'detect/rx.ijs'
load 'detect/matched_filter.ijs'
load 'tensor/cp_als.ijs'
load 'tensor/tucker.ijs'
load 'signal/wavelet.ijs'
load 'signal/fractional_derivative.ijs'
load 'sim/synthetic_scene.ijs'

format_ms=: 3 : 0
  ms=. y * 1000.0
  ": <. (ms * 100.0) % 100.0
)

echo '========================================================================'
echo '                     ALBEDO BENCHMARK SUITE                             '
echo '========================================================================'

NB. Benchmark 1: Thin SVD (50x20 matrix)
A_svd=: ? 50 20 $ 1000
t_svd=: 6!:2 'albedo_linalg_thin_svd A_svd'
echo '1. Thin SVD (50x20 matrix):                  ' , (format_ms t_svd) , ' ms'

NB. Benchmark 2: Cholesky decomposition of 100x100 positive definite matrix
A_chol_raw=: ? 100 100 $ 100
C_chol=: (A_chol_raw albedo_linalg_matmul albedo_linalg_tr A_chol_raw) + 100.0 * albedo_linalg_eye 100
t_chol=: 6!:2 'albedo_linalg_cholesky C_chol'
echo '2. Cholesky Factorization (100x100):         ' , (format_ms t_chol) , ' ms'

NB. Benchmark 3: Spectral Covariance calculation for 500 pixels x 50 bands
X_cov=: ? 50 500 $ 1000
t_cov=: 6!:2 'albedo_stats_covariance_matrix X_cov'
echo '3. Covariance Estimation (50 bands x 500 px):' , (format_ms t_cov) , ' ms'

NB. Benchmark 4: VCA Endmember Extraction (4 endmembers, 1000 pixels x 30 bands)
X_vca=: ? 30 1000 $ 1000
t_vca=: 6!:2 '4 albedo_vca_extract_vca X_vca'
echo '4. VCA Simplex Extraction (p=4, N=1000):     ' , (format_ms t_vca) , ' ms'

NB. Benchmark 5: Reed-Xiaoli Global Anomaly Detector (20 bands x 25 x 20 cube)
X_rx=: ? 20 25 20 $ 1000
t_rx=: 6!:2 'albedo_rx_detect_global_rx X_rx'
echo '5. Reed-Xiaoli Anomaly Detection (20x25x20): ' , (format_ms t_rx) , ' ms'

NB. Benchmark 6: Multi-way Tensor Mode-n Product (20x20x20 tensor)
T_box=: ? 20 20 20 $ 1000
U_mat=: ? 10 20 $ 1000
t_tensor=: 6!:2 '(0; U_mat) albedo_tensor_mode_n_prod T_box'
echo '6. Tensor Mode-0 Contraction (20x20x20):     ' , (format_ms t_tensor) , ' ms'

NB. Benchmark 7: 1D Wavelet DWT decomposition (4096-sample signal, 4 levels)
sig_dwt=: ? 4096 $ 1000
t_dwt=: 6!:2 '4 albedo_wavelet_dwt1_multilevel sig_dwt'
echo '7. 4-Level Discrete Wavelet Transform:       ' , (format_ms t_dwt) , ' ms'

NB. Benchmark 8: Grünwald-Letnikov Fractional Derivative (order 0.5, 224 bands)
spec_frac=: ? 224 $ 1000
t_frac=: 6!:2 '0.5 albedo_frac_fractional_derivative_1d spec_frac'
echo '8. Fractional Differentiation (alpha=0.5):   ' , (format_ms t_frac) , ' ms'

NB. Benchmark 9: Tucker3 HOSVD Tensor Decomposition (15x15x15 cube, rank 4 4 4)
T_tuck=: ? 15 15 15 $ 1000
t_tucker=: 6!:2 '(4 4 4) albedo_tucker_fit_hosvd T_tuck'
echo '9. Tucker3 HOSVD Decomposition (15x15x15):    ' , (format_ms t_tucker) , ' ms'

NB. Benchmark 10: CP-ALS Tensor Factorization (10x10x10 cube, rank 3, 5 iters)
T_cp=: ? 10 10 10 $ 1000
t_cp=: 6!:2 '(3; 5; 1e_4) albedo_cp_als_fit_cp_als T_cp'
echo '10. CP-ALS Factorization (10x10x10, rank 3):  ' , (format_ms t_cp) , ' ms'

echo '========================================================================'
echo 'BENCHMARK RUN COMPLETED'
exit 0
