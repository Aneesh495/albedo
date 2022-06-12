NB. ==============================================================================
NB. Albedo - Spectral Principal Component Analysis (PCA)
NB.
NB. Implements spectral dimension reduction, covariance eigen-decomposition,
NB. explained variance ratios, forward subspace projection, and inverse
NB. low-rank reconstruction.
NB. ==============================================================================

cocurrent 'albedo_pca'

load 'core/linalg.ijs'
load 'core/stats.ijs'

NB. ------------------------------------------------------------------------------
NB. PCA Model Training
NB. ------------------------------------------------------------------------------

NB. Train PCA model from spectral data matrix X of shape (B, N)
NB. Returns: (mu; evals; evecs; explained_var_ratio; cum_var_ratio)
train_pca=: 3 : 0
  X=. y
  'B N'=. $ X
  mu=. albedo_stats_mean_vector X
  Cov=. albedo_stats_covariance_matrix X

  'evals evecs'=. albedo_linalg_eigen_sym Cov
  NB. Protect negative numerical eigenvalues
  evals=. 0.0 >. evals
  tot_var=. (+/ evals) >. 1e_12
  var_ratio=. evals % tot_var
  cum_var=. +/\ var_ratio

  mu; evals; evecs; var_ratio; cum_var
)

NB. ------------------------------------------------------------------------------
NB. Forward Projection & Inverse Reconstruction
NB. ------------------------------------------------------------------------------

NB. Forward PCA projection to k components
NB. Dyadic: (k; pca_model) project_pca X
project_pca=: 4 : 0
  'k model'=. x
  X=. y
  'mu evals evecs var_ratio cum_var'=. model
  assert. k <: # evals [ 'k exceeds total number of components'

  NB. Select top k eigenvectors (columns)
  Vk=. k {"1 evecs

  NB. Center data: Xc = X - mu
  Xc=. X - mu

  NB. Project: Y = Vk^T Xc (shape k x N)
  (albedo_linalg_tr Vk) albedo_linalg_matmul Xc
)

NB. Inverse PCA reconstruction from k components
NB. Dyadic: (k; pca_model) reconstruct_pca Y
reconstruct_pca=: 4 : 0
  'k model'=. x
  Y=. y
  'mu evals evecs var_ratio cum_var'=. model
  Vk=. k {"1 evecs

  NB. Reconstruct: X_hat = Vk Y + mu
  (Vk albedo_linalg_matmul Y) + mu
)

NB. Compute Reconstruction Mean Squared Error (MSE)
NB. Dyadic: X_orig reconstruction_mse X_hat
reconstruction_mse=: 4 : 0
  X_orig=. x
  X_hat=. y
  'B N'=. $ X_orig
  diff=. X_orig - X_hat
  (+/ , *: diff) % (B * N)
)

NB. Find minimum number of components k to retain target variance (e.g. 0.99 for 99%)
NB. Dyadic: target_var_ratio select_k_components pca_model
select_k_components=: 4 : 0
  target_var=. x
  model=. y
  cum_var=. 4 {:: model
  1 + (+/ cum_var < target_var) <. (<: # cum_var)
)
