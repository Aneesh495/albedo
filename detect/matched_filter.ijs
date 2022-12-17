NB. ==============================================================================
NB. Albedo - Spectral Target Detectors: CEM, Matched Filter & ACE
NB.
NB. Implements subpixel spectral target detectors for identifying specific
NB. materials (e.g. methane leaks, rare earth minerals, manmade objects)
NB. against background clutter:
NB.   1. Constrained Energy Minimization (CEM - Farrand & Harsanyi)
NB.   2. Spectral Matched Filter (MF)
NB.   3. Adaptive Cosine / Coherence Estimator (ACE)
NB. ==============================================================================

cocurrent 'base'

load 'core/linalg.ijs'
load 'core/stats.ijs'

NB. ------------------------------------------------------------------------------
NB. Constrained Energy Minimization (CEM)
NB. Filter w = (R^(-1) d) / (d^T R^(-1) d)
NB. Output y_CEM = w^T r
NB. ------------------------------------------------------------------------------

NB. Compute CEM filter weights w for target signature d given correlation matrix R
NB. Dyadic: target_d compute_cem_filter R_corr
compute_cem_filter=: 4 : 0
  d=. x
  R=. y
  B=. # d
  reg_R=. R + (1e_6 * albedo_linalg_eye B)
  inv_R=. %. reg_R

  invR_d=. inv_R albedo_linalg_matmul d
  denom=. +/ d * invR_d
  assert. denom > 1e_12 [ 'Singular target projection in CEM'
  invR_d % denom
)

NB. Run CEM detector on 3D cube (B x H x W)
NB. Dyadic: target_d detect_cem cube
detect_cem=: 4 : 0
  d=. x
  cube=. y
  'B H W'=. $ cube
  X_mat=. (B, H * W) $ cube

  NB. Sample correlation matrix R = (1/N) * X X^T
  R=. (X_mat albedo_linalg_matmul (albedo_linalg_tr X_mat)) % (H * W)
  w_cem=. d compute_cem_filter R

  NB. Filter all pixels: y = w^T X
  scores_1d=. w_cem albedo_linalg_matmul X_mat
  (H, W) $ scores_1d
)

NB. ------------------------------------------------------------------------------
NB. Spectral Matched Filter (MF)
NB. Filter w = (Sigma^(-1) (d - mu)) / ((d - mu)^T Sigma^(-1) (d - mu))
NB. Output y_MF = w^T (r - mu)
NB. ------------------------------------------------------------------------------

NB. Run Spectral Matched Filter on 3D cube
NB. Dyadic: target_d detect_matched_filter cube
detect_matched_filter=: 4 : 0
  d=. x
  cube=. y
  'B H W'=. $ cube
  X_mat=. (B, H * W) $ cube

  mu=. albedo_stats_mean_vector X_mat
  Sigma=. albedo_stats_covariance_matrix X_mat
  reg_Sigma=. Sigma + (1e_6 * albedo_linalg_eye B)
  inv_Sigma=. %. reg_Sigma

  d_cent=. d - mu
  invSig_d=. inv_Sigma albedo_linalg_matmul d_cent
  denom=. (+/ d_cent * invSig_d) >. 1e_12
  w_mf=. invSig_d % denom

  Xc=. X_mat - mu
  scores_1d=. w_mf albedo_linalg_matmul Xc
  (H, W) $ scores_1d
)

NB. ------------------------------------------------------------------------------
NB. Adaptive Cosine / Coherence Estimator (ACE)
NB. Scale-invariant detector:
NB. ACE(r) = (d^T Sigma^(-1) (r - mu))^2 / ((d^T Sigma^(-1) d) * ((r - mu)^T Sigma^(-1) (r - mu)))
NB. Output in [0, 1]
NB. ------------------------------------------------------------------------------

NB. Run ACE detector on 3D cube
NB. Dyadic: target_d detect_ace cube
detect_ace=: 4 : 0
  d=. x
  cube=. y
  'B H W'=. $ cube
  X_mat=. (B, H * W) $ cube

  mu=. albedo_stats_mean_vector X_mat
  Sigma=. albedo_stats_covariance_matrix X_mat
  reg_Sigma=. Sigma + (1e_6 * albedo_linalg_eye B)
  inv_Sigma=. %. reg_Sigma

  d_cent=. d - mu
  invSig_d=. inv_Sigma albedo_linalg_matmul d_cent
  denom_d=. (+/ d_cent * invSig_d) >. 1e_12

  Xc=. X_mat - mu
  NB. Numerator: (d^T Sigma^(-1) Xc)^2
  num=. *: invSig_d albedo_linalg_matmul Xc

  NB. Denominator: denom_d * mahalanobis_sq(x)
  mahal_x=. (mu; Sigma) albedo_stats_batch_mahalanobis_sq X_mat
  denom_tot=. denom_d * (mahal_x >. 1e_12)

  ace_scores=. 0.0 >. 1.0 <. (num % denom_tot)
  (H, W) $ ace_scores
)


albedo_target_detect_detect_cem =: detect_cem
albedo_target_detect_detect_matched_filter =: detect_matched_filter
albedo_target_detect_detect_ace =: detect_ace
