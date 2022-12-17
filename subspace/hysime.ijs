NB. ==============================================================================
NB. Albedo - Hyperspectral Signal Identification by Minimum Error (HySime)
NB.
NB. Unsupervised identification of the intrinsic signal subspace dimension
NB. (number of spectrally distinct endmembers k) without parameter tuning.
NB. Formulates minimum mean squared error (MSE) between projection error
NB. and noise power.
NB. Reference: Bioucas-Dias & Nascimento, IEEE TGRS 2008.
NB. ==============================================================================

cocurrent 'base'

load 'core/linalg.ijs'
load 'core/stats.ijs'

NB. ------------------------------------------------------------------------------
NB. Multiple Linear Regression Noise Estimation
NB. ------------------------------------------------------------------------------

NB. Estimates additive noise vector for each pixel using multiple linear regression (MLR)
NB. where each band is regressed against all other bands.
estimate_noise_mlr=: 3 : 0
  X=. y
  'B N'=. $ X
  Noise=. (B, N) $ 0.0

  for_b. i. B do.
    y_b=. b { X
    NB. Design matrix Z containing all bands except b
    other_bands=. (i. B) -. b
    Z=. other_bands { X

    NB. Regularized regression: beta = (Z Z^T + lambda I)^(-1) Z y_b
    Zt=. albedo_linalg_tr Z
    ZZt=. Z albedo_linalg_matmul Zt
    reg_ZZt=. ZZt + (1e_5 * albedo_linalg_eye (B - 1))
    beta=. reg_ZZt %. (Z albedo_linalg_matmul y_b)

    NB. Residual noise: e = y_b - Z^T beta
    e_b=. y_b - (Zt albedo_linalg_matmul beta)
    Noise=. e_b (<b)} Noise
  end.
  Noise
)

NB. ------------------------------------------------------------------------------
NB. HySime Intrinsic Subspace Dimension Estimation
NB. ------------------------------------------------------------------------------

NB. Estimate intrinsic subspace dimension k from spectral data matrix X (B, N)
NB. Returns: (estimated_k; mse_curve; eigenvectors_U)
hysime_subspace_dim=: 3 : 0
  X=. y
  'B N'=. $ X

  NB. Step 1: Estimate noise and noise correlation matrix Rn
  Noise=. estimate_noise_mlr X
  Rn=. (Noise albedo_linalg_matmul (albedo_linalg_tr Noise)) % N

  NB. Step 2: Signal correlation matrix Rx
  Rx=. (X albedo_linalg_matmul (albedo_linalg_tr X)) % N

  NB. Step 3: Eigen-decomposition of Rx
  'evals_x U'=. albedo_linalg_eigen_sym Rx

  NB. Step 4: Projection error and noise power for each subset of eigenvectors
  NB. p_i = - (u_i^T Rx u_i) + 2 * (u_i^T Rn u_i)
  NB. If p_i < 0, adding eigenvector u_i reduces total MSE.
  costs=. B $ 0.0
  for_i. i. B do.
    ui=. i {"1 U
    NB. Projection variance
    sig_pow=. +/ ui * (Rx albedo_linalg_matmul ui)
    noise_pow=. +/ ui * (Rn albedo_linalg_matmul ui)
    cost=. (_1.0 * sig_pow) + (2.0 * noise_pow)
    costs=. cost (<i)} costs
  end.

  NB. Count negative cost values (where signal dominates noise)
  dim_k=. +/ costs < 0.0
  NB. Clamp k between 1 and B-1
  dim_k=. 1 >. ((<: B) <. dim_k)

  dim_k; costs; U
)


albedo_hysime_hysime_subspace_dim =: hysime_subspace_dim
