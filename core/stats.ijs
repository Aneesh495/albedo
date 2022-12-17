NB. ==============================================================================
NB. Albedo - High-Dimensional Spectral Statistics
NB.
NB. Implements multi-band mean vectors, sample covariance, Ledoit-Wolf style
NB. shrinkage covariance, spectral correlation matrices, batch Mahalanobis
NB. distances, SNR channel estimation, and robust outlier-resistant statistics.
NB. ==============================================================================

cocurrent 'base'

NB. Load linear algebra dependencies
load 'core/linalg.ijs'

NB. ------------------------------------------------------------------------------
NB. Mean Vectors & Centering
NB. ------------------------------------------------------------------------------

NB. Multi-band mean vector from spectral data matrix X of shape (B, N)
NB. B is number of spectral bands, N is number of spatial pixels.
NB. Returns column vector / vector of length B
mean_vector=: 3 : 0
  X=. y
  'B N'=. $ X
  (+/ "1 X) % N
)

NB. Center data matrix by subtracting mean vector from each pixel column
NB. Returns centered matrix (B, N)
center_matrix=: 3 : 0
  X=. y
  mu=. mean_vector X
  X - mu
)

NB. ------------------------------------------------------------------------------
NB. Covariance & Shrinkage Estimation
NB. ------------------------------------------------------------------------------

NB. Sample covariance matrix: Sigma = (1 / (N - 1)) * (X - mu) (X - mu)^T
NB. X has shape (B, N) -> Sigma has shape (B, B)
covariance_matrix=: 3 : 0
  X=. y
  'B N'=. $ X
  assert. N > 1 [ 'Number of samples must be greater than 1 for sample covariance'
  Xc=. center_matrix X
  (Xc albedo_linalg_matmul (albedo_linalg_tr Xc)) % (N - 1)
)

NB. Shrinkage Covariance: (1 - alpha) * Sigma + alpha * (tr(Sigma) / B) * I
NB. Essential for hyperspectral data where B (bands) can be large and matrix ill-conditioned.
NB. Dyadic: alpha shrinkage_cov X
shrinkage_cov=: 4 : 0
  alpha=. x
  X=. y
  B=. {. $ X
  Sigma=. covariance_matrix X
  mu_var=. (albedo_linalg_trace Sigma) % B
  ((1.0 - alpha) * Sigma) + (alpha * mu_var * albedo_linalg_eye B)
)

NB. Spectral Correlation Matrix: R_ij = Sigma_ij / (std_i * std_j)
correlation_matrix=: 3 : 0
  Sigma=. y
  B=. {. $ Sigma
  vars=. albedo_linalg_diag_vec Sigma
  stds=. %: 0.0 >. vars
  NB. Outer product of stds
  denom=. stds */ stds
  NB. Protect against zero division
  denom=. denom + (denom = 0.0)
  Sigma % denom
)

NB. ------------------------------------------------------------------------------
NB. Mahalanobis Distance
NB. ------------------------------------------------------------------------------

NB. Single vector Mahalanobis distance squared: (x - mu)^T Sigma^(-1) (x - mu)
NB. Dyadic: (mu; inv_Sigma) mahalanobis_sq x
mahalanobis_sq=: 4 : 0
  'mu inv_Sigma'=. x
  x_vec=. y
  diff=. x_vec - mu
  +/ diff * (inv_Sigma albedo_linalg_matmul diff)
)

NB. Batch Mahalanobis distance squared for all pixels in X (B, N)
NB. Uses Cholesky whitening: Sigma = L L^T -> L^(-1) (X - mu) -> sum of squares
NB. Dyadic: (mu; Sigma) batch_mahalanobis_sq X
batch_mahalanobis_sq=: 4 : 0
  'mu Sigma'=. x
  X=. y
  'B N'=. $ X
  Xc=. X - mu

  NB. Add small diagonal jitter for numerical stability
  reg_Sigma=. Sigma + (1e_7 * albedo_linalg_eye B)
  inv_Sigma=. %. reg_Sigma

  NB. Compute (Xc^T inv_Sigma Xc) diagonal efficiently
  NB. Y = inv_Sigma matmul Xc (shape B, N)
  Y=. inv_Sigma albedo_linalg_matmul Xc
  NB. Column-wise inner product
  +/ Xc * Y
)

NB. ------------------------------------------------------------------------------
NB. Signal-to-Noise Ratio (SNR) Estimation
NB. ------------------------------------------------------------------------------

NB. Estimates per-band SNR (dB) using spatial first-difference noise proxy
NB. For 3D cube (B, H, W):
NB. Signal is mean radiance across scene.
NB. Noise variance is estimated from horizontal spatial differences: Delta = X[h, w+1] - X[h, w]
NB. var(noise) = var(Delta) / 2
snr_per_band=: 3 : 0
  Cube=. y
  'B H W'=. $ Cube
  snr_db=. B $ 0.0

  for_b. i. B do.
    band=. b { Cube
    mean_sig=. (+/ , band) % (H * W)
    NB. Horizontal spatial differences
    diff_h=. (1 {"1 band) - (_1 {"1 band)
    diff_var=. (+/ *: , diff_h) % ((# , diff_h) - 1)
    noise_var=. 0.5 * diff_var
    if. noise_var > 1e_12 do.
      snr_linear=. (*: mean_sig) % noise_var
      snr_db=. (10.0 * 10 ^. (1.0 >. snr_linear)) (<b)} snr_db
    else.
      snr_db=. 60.0 (<b)} snr_db
    end.
  end.
  snr_db
)

NB. ------------------------------------------------------------------------------
NB. Robust Statistics
NB. ------------------------------------------------------------------------------

NB. Median Absolute Deviation (MAD) of a 1D vector
NB. MAD = median(|x - median(x)|) * 1.4826
median=: 3 : 0
  s=. /:~ y
  n=. # s
  half=. <. n % 2
  if. 2 | n do.
    half { s
  else.
    -:(half { s) + (<: half) { s
  end.
)

mad=: 3 : 0
  v=. y
  med=. median v
  devs=. | v - med
  1.4826 * median devs
)

NB. Robust z-score using median and MAD
robust_zscore=: 3 : 0
  v=. y
  med=. median v
  scale=. mad v
  if. scale = 0.0 do. v * 0.0 else. (v - med) % scale end.
)


albedo_stats_mean_vector =: mean_vector
albedo_stats_center_matrix =: center_matrix
albedo_stats_covariance_matrix =: covariance_matrix
albedo_stats_shrinkage_cov =: shrinkage_cov
albedo_stats_correlation_matrix =: correlation_matrix
albedo_stats_mahalanobis_sq =: mahalanobis_sq
albedo_stats_batch_mahalanobis_sq =: batch_mahalanobis_sq
albedo_stats_channel_snr =: channel_snr
albedo_stats_cube_snr_profile =: cube_snr_profile
albedo_stats_median =: median
albedo_stats_mad =: mad
albedo_stats_winsorize =: winsorize
