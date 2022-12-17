NB. ==============================================================================
NB. ALBEDO - High-Dimensional Hyperspectral Imaging Engine
NB. scripts/benchmarks.ijs : Throughput and Latency Microbenchmarks
NB. Measures execution time across linear algebra, tensor decompositions,
NB. unmixing optimization, and anomaly detection kernels.
NB. ==============================================================================

load 'core/linalg.ijs'
load 'core/linalg_ext.ijs'
load 'core/stats.ijs'
load 'core/tensor.ijs'
load 'subspace/pca.ijs'
load 'subspace/mnf.ijs'
load 'unmix/vca.ijs'
load 'unmix/fcls.ijs'
load 'detect/rx.ijs'
load 'detect/matched_filter.ijs'
load 'tensor/cp_als.ijs'
load 'tensor/tucker.ijs'
load 'signal/wavelet.ijs'

echo '========================================================================'
echo '                     ALBEDO BENCHMARK SUITE                             '
echo '========================================================================'

NB. Benchmark 1: SVD of 500x100 matrix
n =: 500
p =: 100
A =: (n, p) $ (i. n * p) % (n * p)
t_start =: (6!:1 '')
svd_res =: thin_svd_albedo_linalg_ A
t_svd =: (6!:1 '') - t_start
echo '1. Thin SVD (500x100 matrix):             ' , (": 1000 * t_svd) , ' ms'

NB. Benchmark 2: Cholesky decomposition of 150x150 positive definite matrix
K =: 150
C =: ((K, K) $ 0.5) + =/ ~ i. K
t_start =: (6!:1 '')
L =: cholesky_albedo_linalg_ C
t_chol =: (6!:1 '') - t_start
echo '2. Cholesky Factorization (150x150):      ' , (": 1000 * t_chol) , ' ms'

NB. Benchmark 3: Spectral Covariance calculation for 5000 pixels x 50 bands
X =: (5000, 50) $ (i. 5000 * 50) % (5000 * 50)
t_start =: (6!:1 '')
cov =: sample_covariance_albedo_stats_ X
t_cov =: (6!:1 '') - t_start
echo '3. Covariance Estimation (5000x50):       ' , (": 1000 * t_cov) , ' ms'

NB. Benchmark 4: VCA Endmember Extraction (4 endmembers, 1000 pixels x 30 bands)
Y =: (30, 1000) $ (i. 30 * 1000) % (30 * 1000)
t_start =: (6!:1 '')
vca_res =: 4 vca_albedo_vca_ Y
t_vca =: (6!:1 '') - t_start
echo '4. VCA Endmember Extraction (p=4, N=1000):' , (": 1000 * t_vca) , ' ms'

NB. Benchmark 5: Reed-Xiaoli Global Anomaly Detector (2000 pixels x 20 bands)
X_rx =: (2000, 20) $ (i. 2000 * 20) % (2000 * 20)
t_start =: (6!:1 '')
rx_scores =: rx_detector_albedo_rx_ X_rx
t_rx =: (6!:1 '') - t_start
echo '5. Reed-Xiaoli Anomaly Detection (2000x20): ' , (": 1000 * t_rx) , ' ms'

NB. Benchmark 6: Multi-way Tensor Matricization & Mode-n Product (30x30x20 tensor)
T =: (30 30 20) $ (i. 30 * 30 * 20) % (30 * 30 * 20)
U =: (15, 30) $ 0.1
t_start =: (6!:1 '')
T_prod =: 0 mode_n_product_albedo_tensor_ (T ; U)
t_tensor =: (6!:1 '') - t_start
echo '6. Tensor Mode-n Contraction (30x30x20):   ' , (": 1000 * t_tensor) , ' ms'

NB. Benchmark 7: 1D Wavelet DWT decomposition (4096-sample signal, 5 levels)
sig =: ? 4096 # 0
t_start =: (6!:1 '')
dwt_res =: 5 dwt1_multilevel_albedo_wavelet_ sig
t_dwt =: (6!:1 '') - t_start
echo '7. 5-Level Discrete Wavelet Transform:     ' , (": 1000 * t_dwt) , ' ms'

echo '========================================================================'
echo 'BENCHMARK RUN COMPLETED'
exit 0
