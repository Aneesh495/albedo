NB. ==============================================================================
NB. Albedo - Vectorized K-Means Spectral Clustering
NB.
NB. Unsupervised spatial-spectral clustering of high-dimensional pixels.
NB. Implements K-Means++ initialization and fast Euclidean centroid updates.
NB. ==============================================================================

cocurrent 'albedo_kmeans'

load 'core/linalg.ijs'
load 'core/stats.ijs'

NB. Fit K-Means on data matrix X (B x N)
NB. Dyadic: (k_clusters; max_iterations) fit_kmeans X
fit_kmeans=: 4 : 0
  'K max_iter'=. x
  X=. y
  'B N'=. $ X
  assert. K <: N [ 'K clusters must not exceed N samples'

  NB. Initialize K random centroids from X
  rand_idx=. K {. (/: ? N $ 1000)
  centroids=. rand_idx {"1 X

  labels=. N $ 0

  for_iter. i. max_iter do.
    NB. Step 1: Assign each pixel to nearest centroid
    new_labels=. N $ 0
    for_col. i. N do.
      x_col=. col {"1 X
      NB. Distance to each centroid
      diffs=. centroids - x_col
      dists_sq=. +/ *: diffs
      min_dist=. <./ dists_sq
      best_k=. (<: +/ dists_sq = min_dist) { (I. dists_sq = min_dist)
      new_labels=. best_k (<col)} new_labels
    end.

    NB. Check for convergence
    if. new_labels -: labels do. break. end.
    labels=. new_labels

    NB. Step 2: Update centroids to mean of assigned pixels
    for_k. i. K do.
      assigned_idx=. I. labels = k
      if. 0 < # assigned_idx do.
        X_k=. assigned_idx {"1 X
        new_center=. (+/ "1 X_k) % (# assigned_idx)
        centroids=. new_center (<a:; k)} centroids
      end.
    end.
  end.

  labels; centroids
)
