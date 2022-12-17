NB. ==============================================================================
NB. ALBEDO - High-Dimensional Hyperspectral Imaging Engine
NB. detect/kernel_rx.ijs : Nonlinear Kernel Reed-Xiaoli (KRX) Anomaly Detector
NB. Maps spectral pixels into reproducing kernel Hilbert space (RKHS) via Gaussian RBF
NB. kernels to detect nonlinear manifolds and subtle non-Gaussian subpixel anomalies.
NB. ==============================================================================

cocurrent 'base'

load 'core/linalg.ijs'

NB. Gaussian Radial Basis Function (RBF) kernel between two spectral vectors:
NB. k(x, y) = exp( - gamma * ||x - y||^2 )
rbf_kernel =: 4 : 0
  gamma =. x
  d2 =. +/ *: y - x NB. dyad handled if called correctly
  ^. - gamma * d2
)

NB. Compute full Gram matrix K (N x N) for matrix X (N pixels x B bands)
NB. gamma: kernel parameter
rbf_gram_matrix =: 4 : 0
  gamma =. x
  X =. y
  n =. # X
  NB. Vectorized squared Euclidean distances
  NB. ||xi - xj||^2 = ||xi||^2 + ||xj||^2 - 2 xi . xj
  norms2 =. +/"1 *: X
  dot_prod =. X +/ .* |: X
  dist2 =. (norms2 +/ norms2) - 2.0 * dot_prod
  NB. Zero negative values due to floating-point rounding
  dist2 =. 0.0 >. dist2
  ^ - gamma * dist2
)

NB. Center Gram matrix in Hilbert space:
NB. K_tilde = K - 1_N K - K 1_N + 1_N K 1_N
NB. where 1_N is (1/N) * J_N
center_gram_matrix =: 3 : 0
  n =. # y
  H =. (=/ ~ i. n) - (% n)
  H +/ .* y +/ .* H
)

NB. Run Kernel RX Anomaly Detection on spectral dataset:
NB. x: (gamma ; reg_lambda)
NB. y: X (N pixels x B bands)
NB. Returns: anomaly score vector of length N
kernel_rx =: 4 : 0
  'gamma lambda' =. x
  X =. y
  n =. # X

  NB. 1. Compute RBF Gram matrix K
  K =. gamma rbf_gram_matrix X

  NB. 2. Center K in RKHS feature space
  K_tilde =. center_gram_matrix K

  NB. 3. Regularized pseudo-inverse: (K_tilde + lambda * I)^-1
  eye_n =. =/ ~ i. n
  K_reg =. K_tilde + lambda * eye_n
  K_inv =. albedo_linalg_pinv_svd K_reg

  NB. 4. Kernel RX score for pixel i is the i-th diagonal element of K_tilde * K_inv * K_tilde
  NB. score = diag(K_tilde * K_inv * K_tilde)
  quad =. (K_tilde +/ .* K_inv) +/ .* K_tilde
  (<"1 (i. n) ,. (i. n)) { quad
)


NB. Exported Aliases
albedo_krx_rbf_gram_matrix =: rbf_gram_matrix
albedo_krx_center_gram_matrix =: center_gram_matrix
albedo_krx_kernel_rx =: kernel_rx
