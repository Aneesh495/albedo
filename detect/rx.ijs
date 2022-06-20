NB. ==============================================================================
NB. Albedo - Global Reed-Xiaoli (RX) Multivariate Anomaly Detector
NB.
NB. Implements the benchmark Reed-Xiaoli (RX) algorithm (Reed & Yu, IEEE T-AES 1990)
NB. for unsupervised detection of spectral anomalies and rare materials.
NB. Evaluates quadratic Mahalanobis metric:
NB. D_RX(r) = (r - mu_b)^T Sigma_b^(-1) (r - mu_b)
NB. ==============================================================================

cocurrent 'albedo_rx'

load 'core/linalg.ijs'
load 'core/stats.ijs'

NB. ------------------------------------------------------------------------------
NB. Global RX Detection on 3D Hyperspectral Cube
NB. ------------------------------------------------------------------------------

NB. Compute global RX anomaly score map for 3D cube (B x H x W)
NB. Returns 2D matrix of anomaly scores (H x W)
detect_global_rx=: 3 : 0
  cube=. y
  'B H W'=. $ cube
  N=. H * W
  X_mat=. (B, N) $ cube

  NB. Background statistics across entire scene
  mu_b=. albedo_stats_mean_vector X_mat
  Sigma_b=. albedo_stats_covariance_matrix X_mat

  NB. Batch Mahalanobis distance squared across all pixels
  scores_1d=. (mu_b; Sigma_b) albedo_stats_batch_mahalanobis_sq X_mat

  NB. Reshape 1D scores back to 2D spatial map (H x W)
  (H, W) $ scores_1d
)

NB. ------------------------------------------------------------------------------
NB. Anomaly Mask Thresholding
NB. ------------------------------------------------------------------------------

NB. Threshold anomaly map at specified percentile (e.g. 0.99 for top 1% anomalies)
NB. Dyadic: percentile threshold_anomalies score_map
threshold_anomalies=: 4 : 0
  pct=. x
  scores=. y
  'H W'=. $ scores
  vals=. /:~ , scores
  n=. # vals
  thresh=. (<. pct * n) { vals
  scores >: thresh
)
