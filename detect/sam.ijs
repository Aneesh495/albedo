NB. ==============================================================================
NB. Albedo - Spectral Angle Mapper (SAM) & Spectral Information Divergence (SID)
NB.
NB. Implements hyperdimensional angle and information-theoretic divergence metrics:
NB.   1. Spectral Angle Mapper (SAM): angle theta = arccos((r . t) / (||r|| ||t||))
NB.      Scale, illumination, and topography invariant.
NB.   2. Spectral Information Divergence (SID - Chang, 1999):
NB.      Kullback-Leibler symmetric divergence between normalized spectral probability vectors.
NB. ==============================================================================

cocurrent 'base'

load 'core/linalg.ijs'

NB. ------------------------------------------------------------------------------
NB. Spectral Angle Mapper (SAM)
NB. ------------------------------------------------------------------------------

NB. Compute spectral angle (radians) between two 1D spectral vectors
NB. Dyadic: spec_a spectral_angle spec_b
spectral_angle=: 4 : 0
  a=. x
  b=. y
  dot=. +/ a * b
  norm_a=. albedo_linalg_norm2 a
  norm_b=. albedo_linalg_norm2 b
  denom=. (norm_a * norm_b) >. 1e_12
  cos_theta=. _1.0 >. 1.0 <. (dot % denom)
  _2&o. cos_theta
)

NB. Batch Spectral Angle Mapper across 3D cube (B x H x W) against target vector t
NB. Returns 2D angle map (H x W) in radians
NB. Dyadic: target_t sam_cube cube
sam_cube=: 4 : 0
  t=. x
  cube=. y
  'B H W'=. $ cube
  X_mat=. (B, H * W) $ cube

  norm_t=. albedo_linalg_norm2 t
  dot_prods=. t albedo_linalg_matmul X_mat
  norm_pixels=. %: +/ *: X_mat
  denom=. (norm_t * norm_pixels) >. 1e_12

  cos_thetas=. _1.0 >. 1.0 <. (dot_prods % denom)
  angles_1d=. _2&o. cos_thetas
  (H, W) $ angles_1d
)

NB. ------------------------------------------------------------------------------
NB. Spectral Information Divergence (SID)
NB. ------------------------------------------------------------------------------

NB. Compute SID between two 1D spectral vectors
NB. Normalizes spectra into probability vectors p and q:
NB. SID(p, q) = sum(p_i * log(p_i / q_i)) + sum(q_i * log(q_i / p_i))
NB. Dyadic: spec_a spectral_information_divergence spec_b
spectral_information_divergence=: 4 : 0
  a=. 0.0 >. x
  b=. 0.0 >. y
  NB. Normalize to sum-to-one probability distributions
  p=. a % ((+/ a) >. 1e_12)
  q=. b % ((+/ b) >. 1e_12)

  NB. Protect against zero entries in log
  eps=. 1e_10
  p=. (p + eps) % (+/ p + eps)
  q=. (q + eps) % (+/ q + eps)

  NB. KL divergences D(p || q) and D(q || p)
  kl_pq=. +/ p * ^. (p % q)
  kl_qp=. +/ q * ^. (q % p)
  0.0 >. (kl_pq + kl_qp)
)


albedo_similarity_spectral_angle =: spectral_angle
albedo_similarity_spectral_information_divergence =: spectral_information_divergence
albedo_similarity_sam_cube =: sam_cube
