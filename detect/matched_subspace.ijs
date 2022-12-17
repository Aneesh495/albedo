NB. ==============================================================================
NB. ALBEDO - High-Dimensional Hyperspectral Imaging Engine
NB. detect/matched_subspace.ijs : Adaptive Matched Subspace Detector (AMSD)
NB. Implements GLRT subspace projection detector against composite hypotheses
NB. where both target and background clutter occupy low-dimensional linear subspaces.
NB. ==============================================================================

cocurrent 'base'

load 'core/linalg.ijs'

NB. Construct orthogonal projection matrix onto orthogonal complement of column space S:
NB. P_S_perp = I - S * (S^T * S)^-1 * S^T
proj_perp =: 3 : 0
  S =. y
  'b p' =. $ S
  eye =. =/ ~ i. b
  NB. S^T * S
  StS =. (|: S) +/ .* S
  NB. Invert or pseudoinvert
  inv_StS =. pinv_albedo_linalg_ StS
  proj =. S +/ .* (inv_StS +/ .* |: S)
  eye - proj
)

NB. Compute AMSD detection statistic for a single pixel vector x (length B):
NB. x: (P_B_perp ; P_Z_perp) where Z = [T , B]
NB. y: pixel vector x
amsd_pixel =: 4 : 0
  'PB_perp PZ_perp' =. x
  vec =. y
  NB. Quadratic forms
  num =. (+/ vec * PB_perp +/ .* vec) - (+/ vec * PZ_perp +/ .* vec)
  denom =. (+/ vec * PZ_perp +/ .* vec) + 1e_12
  num % denom
)

NB. Batch AMSD detection across image matrix X (N pixels x B bands):
NB. x: (T_subspace ; B_subspace) where T is (B x p_t) and B is (B x p_b)
NB. y: X (N pixels x B bands)
NB. Returns: anomaly score vector of length N
amsd_detector =: 4 : 0
  'T B' =. x
  X =. y
  NB. Composite subspace Z = [T , B]
  Z =. T ,. B
  PB_perp =. proj_perp B
  PZ_perp =. proj_perp Z

  NB. Map over rows of X
  (PB_perp ; PZ_perp) & amsd_pixel "1 X
)


NB. Exported Aliases
albedo_amsd_amsd_detector =: amsd_detector
