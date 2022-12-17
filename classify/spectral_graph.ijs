NB. ==============================================================================
NB. Albedo - Graph-Based Spectral Clustering & Manifold Learning
NB.
NB. Implements Laplacian Eigenmaps and spectral graph clustering:
NB.   1. Affinity matrix W via Gaussian RBF kernel.
NB.   2. Degree matrix D and unnormalized/normalized graph Laplacians.
NB.   3. Spectral embedding via bottom eigenvectors of L_sym = D^(-1/2) L D^(-1/2).
NB.   4. Manifold clustering on nonlinear spectral geometries.
NB. ==============================================================================

cocurrent 'base'

load 'core/linalg.ijs'
load 'classify/kmeans.ijs'

NB. ------------------------------------------------------------------------------
NB. Affinity Matrix Construction
NB. ------------------------------------------------------------------------------

NB. Compute Gaussian RBF affinity matrix W (N x N) for data matrix X (B x N)
NB. W_ij = exp(- ||x_i - x_j||^2 / (2 sigma^2))
NB. Dyadic: sigma compute_rbf_affinity X
compute_rbf_affinity=: 4 : 0
  sigma=. x
  X=. y
  'B N'=. $ X
  two_sigma_sq=. 2.0 * *: sigma

  NB. Compute pairwise squared Euclidean distance matrix
  W=. (N, N) $ 0.0
  for_i. i. N do.
    xi=. i {"1 X
    diffs=. X - xi
    dists_sq=. +/ *: diffs
    weights=. ^ - (dists_sq % two_sigma_sq)
    W=. weights (<i)} W
  end.
  NB. Zero out self-affinity diagonal
  (1.0 - albedo_linalg_eye N) * W
)

NB. ------------------------------------------------------------------------------
NB. Normalized Symmetric Graph Laplacian
NB. L_sym = I - D^(-1/2) W D^(-1/2)
NB. ------------------------------------------------------------------------------

NB. Compute normalized symmetric Graph Laplacian L_sym
NB. Returns (L_sym; D_sqrt_inv)
normalized_laplacian=: 3 : 0
  W=. y
  N=. {. $ W
  degrees=. +/ "1 W
  inv_sqrt_d=. % %: degrees >. 1e_8
  D_inv_sqrt=. albedo_linalg_diag_mat inv_sqrt_d

  NB. Normalized affinity: W_norm = D^(-1/2) W D^(-1/2)
  W_norm=. D_inv_sqrt albedo_linalg_matmul W albedo_linalg_matmul D_inv_sqrt

  NB. L_sym = I - W_norm
  L_sym=. (albedo_linalg_eye N) - W_norm
  L_sym; inv_sqrt_d
)

NB. ------------------------------------------------------------------------------
NB. Spectral Embedding & Clustering
NB. ------------------------------------------------------------------------------

NB. Compute k-dimensional spectral embedding using the smallest eigenvectors of L_sym
NB. Dyadic: k_dim spectral_embedding X
spectral_embedding=: 4 : 0
  k=. x
  X=. y
  'B N'=. $ X
  sigma=. 10.0

  W=. sigma compute_rbf_affinity X
  'L_sym inv_sqrt_d'=. normalized_laplacian W

  NB. Eigen-decomposition of L_sym (eigen_sym returns descending order)
  'evals V'=. albedo_linalg_eigen_sym L_sym

  NB. Smallest k eigenvectors correspond to the tail of sorted evals
  bottom_indices=. (<: N) - i. k
  U=. bottom_indices {"1 V

  NB. Normalize rows of embedding to unit length
  row_norms=. (albedo_linalg_norm2 "1 U) >. 1e_12
  U_norm=. U % row_norms
  albedo_linalg_tr U_norm
)

NB. Full Spectral Graph Clustering
NB. Dyadic: (k_clusters; sigma) cluster_spectral X
cluster_spectral=: 4 : 0
  'K sigma'=. x
  X=. y
  'B N'=. $ X

  W=. sigma compute_rbf_affinity X
  'L_sym inv_sqrt_d'=. normalized_laplacian W
  'evals V'=. albedo_linalg_eigen_sym L_sym

  NB. Extract K smallest eigenvectors
  bottom_indices=. (<: N) - i. K
  U=. bottom_indices {"1 V
  row_norms=. (albedo_linalg_norm2 "1 U) >. 1e_12
  U_norm=. U % row_norms
  Y_norm=. albedo_linalg_tr U_norm

  NB. Run K-Means on normalized spectral coordinates (K x N)
  'labels centroids'=. (K; 25) albedo_kmeans_fit_kmeans Y_norm
  labels; Y_norm
)


albedo_spectral_graph_compute_rbf_affinity =: compute_rbf_affinity
albedo_spectral_graph_normalized_laplacian =: normalized_laplacian
albedo_spectral_graph_spectral_embedding =: spectral_embedding
albedo_spectral_graph_cluster_spectral =: cluster_spectral
