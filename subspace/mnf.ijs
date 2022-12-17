NB. ==============================================================================
NB. Albedo - Minimum Noise Fraction (MNF) Transform
NB.
NB. Implements Green's Maximum/Minimum Noise Fraction transform:
NB. 1. Spatial difference noise covariance estimation Sigma_n.
NB. 2. Noise whitening transformation F such that F^T Sigma_n F = I.
NB. 3. Principal component rotation of noise-whitened data.
NB. 4. Component ranking strictly by decreasing Signal-to-Noise Ratio (SNR).
NB. ==============================================================================

cocurrent 'base'

load 'core/linalg.ijs'
load 'core/stats.ijs'
load 'core/tensor.ijs'

NB. ------------------------------------------------------------------------------
NB. Spatial Difference Noise Covariance Estimation
NB. ------------------------------------------------------------------------------

NB. Estimate noise covariance matrix Sigma_n (B x B) from 3D cube (B x H x W)
NB. Uses horizontal spatial shift difference: Delta = X[h, w+1] - X[h, w]
NB. Sigma_n = 0.5 * Cov(Delta)
estimate_noise_covariance=: 3 : 0
  Cube=. y
  'B H W'=. $ Cube
  assert. W > 1 [ 'Width must be > 1 to compute spatial differences'

  NB. Compute horizontal differences across lines
  diff_cube=. (1 {"2 Cube) - (_1 {"2 Cube)
  NB. diff_cube has shape (B, H, W-1)
  N_diff=. H * (W - 1)
  diff_mat=. (B, N_diff) $ diff_cube

  NB. Compute sample covariance of differences
  diff_cov=. albedo_stats_covariance_matrix diff_mat
  0.5 * diff_cov
)

NB. ------------------------------------------------------------------------------
NB. MNF Transform Model Training
NB. ------------------------------------------------------------------------------

NB. Train MNF model on 3D cube (B x H x W)
NB. Returns: (mu; F_white; G_rot; T_mnf; snr_evals)
train_mnf=: 3 : 0
  Cube=. y
  'B H W'=. $ Cube
  N=. H * W
  X_mat=. (B, N) $ Cube
  mu=. albedo_stats_mean_vector X_mat
  Sigma_total=. albedo_stats_covariance_matrix X_mat

  NB. Step 1: Estimate noise covariance Sigma_n
  Sigma_n=. estimate_noise_covariance Cube
  NB. Regularize noise covariance for invertibility
  reg_Sigma_n=. Sigma_n + (1e_6 * albedo_linalg_eye B)

  NB. Step 2: Noise Whitening matrix F
  NB. Sigma_n = V_n Lambda_n V_n^T -> F = V_n Lambda_n^(-1/2)
  'evals_n V_n'=. albedo_linalg_eigen_sym reg_Sigma_n
  inv_sqrt_evals=. % %: 1e_8 >. evals_n
  F_white=. V_n albedo_linalg_matmul (albedo_linalg_diag_mat inv_sqrt_evals)

  NB. Step 3: Transform total covariance to whitened space
  NB. Sigma_adj = F^T Sigma_total F
  Sigma_adj=. (albedo_linalg_tr F_white) albedo_linalg_matmul Sigma_total albedo_linalg_matmul F_white

  NB. Step 4: Eigen-decomposition of whitened covariance
  NB. Sigma_adj = G Lambda_adj G^T
  'evals_adj G_rot'=. albedo_linalg_eigen_sym Sigma_adj

  NB. Step 5: Full forward MNF projection matrix T = G^T F^T (shape B x B)
  T_mnf=. (albedo_linalg_tr G_rot) albedo_linalg_matmul (albedo_linalg_tr F_white)

  mu; F_white; G_rot; T_mnf; evals_adj
)

NB. ------------------------------------------------------------------------------
NB. Forward & Inverse MNF Transforms
NB. ------------------------------------------------------------------------------

NB. Forward MNF transform to k components
NB. Dyadic: (k; mnf_model) project_mnf X_mat
project_mnf=: 4 : 0
  'k model'=. x
  X_mat=. y
  'mu F_white G_rot T_mnf evals_adj'=. model
  assert. k <: # evals_adj [ 'k exceeds total MNF bands'

  Tk=. k {. T_mnf
  Xc=. X_mat - mu
  Tk albedo_linalg_matmul Xc
)

NB. Inverse MNF transform from k components
NB. Dyadic: (k; mnf_model) reconstruct_mnf Y_mnf
reconstruct_mnf=: 4 : 0
  'k model'=. x
  Y_mnf=. y
  'mu F_white G_rot T_mnf evals_adj'=. model
  Tk=. k {. T_mnf

  NB. Pseudoinverse of Tk: Tk^+ = Tk^T (Tk Tk^T)^(-1)
  Tk_pinv=. albedo_linalg_pinv_svd Tk
  (Tk_pinv albedo_linalg_matmul Y_mnf) + mu
)


albedo_mnf_train_mnf =: train_mnf
albedo_mnf_project_mnf =: project_mnf
