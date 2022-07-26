NB. ==============================================================================
NB. Albedo - Robust Tensor Principal Component Analysis (RTPCA)
NB.
NB. Decomposes 3-way hyperspectral tensor X into:
NB.   X = L + S + E
NB. where:
NB.   - L is low-rank background tensor (spatially & spectrally correlated)
NB.   - S is sparse anomaly tensor (localized subpixel targets, methane plumes)
NB.   - E is Gaussian sensor noise
NB. Solves via Alternating Direction Method of Multipliers (ADMM) with
NB. Tensor Singular Value Thresholding (t-SVT) and soft-thresholding.
NB. ==============================================================================

cocurrent 'albedo_rtpca'

load 'core/tensor.ijs'
load 'core/linalg.ijs'

NB. Soft-thresholding operator S_tau(x) = sign(x) * max(0, |x| - tau)
soft_thresh=: 4 : 0
  tau=. x
  x_val=. y
  sign=. (x_val > 0.0) - (x_val < 0.0)
  sign * (0.0 >. (| x_val) - tau)
)

NB. Singular Value Thresholding (SVT) on 2D matrix M: D_tau(M) = U S_tau(Sigma) V^T
matrix_svt=: 4 : 0
  tau=. x
  M=. y
  'svals U V'=. albedo_linalg_svd_power M
  svals_thresh=. 0.0 >. (svals - tau)
  U albedo_linalg_matmul (albedo_linalg_diag_mat svals_thresh) albedo_linalg_matmul (albedo_linalg_tr V)
)

NB. Tensor Singular Value Thresholding along mode 0 (spectral unfolding)
tensor_svt_mode0=: 4 : 0
  tau=. x
  T=. y
  dims=. $ T
  T0=. 0 albedo_tensor_tensor_unfold T
  T0_thresh=. tau matrix_svt T0
  (0; dims) albedo_tensor_fold T0_thresh
)

NB. Fit Robust Tensor PCA on 3D tensor X of shape (B, H, W)
NB. Dyadic: (lambda_sparse; max_iterations; tolerance) fit_rtpca X
fit_rtpca=: 4 : 0
  'lambda_s max_iter tol'=. x
  X=. y
  dims=. $ X
  norm_X=. albedo_tensor_tensor_fnorm X

  NB. Initialize L (low-rank), S (sparse), Y (multiplier)
  L=. dims $ 0.0
  S=. dims $ 0.0
  Y=. dims $ 0.0

  mu_admm=. 1.0 % ((>./ dims) >. 1.0)
  rho_rate=. 1.1

  for_iter. i. max_iter do.
    NB. Step 1: Update Low-Rank tensor L = D_{1/mu}(X - S + (1/mu) Y)
    arg_L=. X - S + (Y % mu_admm)
    tau_L=. 1.0 % mu_admm
    L_new=. tau_L tensor_svt_mode0 arg_L

    NB. Step 2: Update Sparse tensor S = S_{lambda/mu}(X - L + (1/mu) Y)
    arg_S=. X - L_new + (Y % mu_admm)
    tau_S=. lambda_s % mu_admm
    S_new=. tau_S soft_thresh arg_S

    NB. Step 3: Update Lagrange multiplier Y = Y + mu * (X - L - S)
    residual=. X - L_new - S_new
    Y=. Y + (mu_admm * residual)

    NB. Check convergence: ||X - L - S||_F / ||X||_F < tol
    rel_res=. (albedo_tensor_tensor_fnorm residual) % (norm_X >. 1e_8)
    L=. L_new
    S=. S_new

    if. rel_res < tol do. break. end.

    NB. Accelerate penalty parameter mu
    mu_admm=. 1e6 <. (mu_admm * rho_rate)
  end.

  NB. Returns (L_background; S_sparse_anomalies)
  L; S
)
