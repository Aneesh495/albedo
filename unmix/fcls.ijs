NB. ==============================================================================
NB. Albedo - Fully Constrained Least Squares (FCLS) Spectral Unmixing
NB.
NB. Solves linear spectral unmixing: y = M s + e
NB. subject to:
NB.   1. Non-negativity Constraint (ANC): s_i >= 0, for all i = 1..p
NB.   2. Abundance Sum-to-One Constraint (ASC): sum(s_i) = 1
NB. Implements Heinz & Chang (IEEE TGRS, 2001) augmented active-set solver.
NB. ==============================================================================

cocurrent 'albedo_fcls'

load 'unmix/nnls.ijs'
load 'core/linalg.ijs'

NB. ------------------------------------------------------------------------------
NB. Single Pixel FCLS Solver
NB. ------------------------------------------------------------------------------

NB. Augmented system formulation:
NB. Appends a row of delta * ones to M and delta to y
NB. where delta is a scaling parameter weighting the sum-to-one constraint.
NB. Dyadic: M solve_fcls y_pixel
solve_fcls=: 4 : 0
  M=. x
  y_pixel=. y
  'B p'=. $ M
  delta=. 1e2

  NB. Augmented endmember matrix M_aug of shape (B + 1, p)
  aug_row=. p $ delta
  M_aug=. M , aug_row

  NB. Augmented pixel observation y_aug of length (B + 1)
  y_aug=. y_pixel , delta

  NB. Solve NNLS on augmented system: ensures ANC (s >= 0) and forces ASC (sum s ~ 1)
  s_raw=. M_aug albedo_nnls_solve_nnls y_aug

  NB. Exact normalization to guarantee sum-to-one to machine precision
  sum_s=. (+/ s_raw) >. 1e_12
  s_raw % sum_s
)

NB. ------------------------------------------------------------------------------
NB. Batch Hyperspectral Cube Unmixing
NB. ------------------------------------------------------------------------------

NB. Unmix full 3D data cube (B x H x W) against endmember matrix M (B x p)
NB. Returns abundance fraction cube of shape (p x H x W)
NB. Dyadic: M unmix_cube cube_3d
unmix_cube=: 4 : 0
  M=. x
  cube=. y
  'B H W'=. $ cube
  'BM p'=. $ M
  assert. B = BM [ 'Spectral band count mismatch between cube and endmembers'

  abundances=. (p, H, W) $ 0.0

  for_h. i. H do.
    for_w. i. W do.
      pixel=. cube <a:; h; w
      s_pixel=. M solve_fcls pixel
      abundances=. s_pixel (<a:; h; w)} abundances
    end.
  end.
  abundances
)

NB. Compute Reconstruction Root Mean Square Error (RMSE) of unmixing
NB. Dyadic: (M; abundances) unmix_rmse cube_3d
unmix_rmse=: 4 : 0
  'M A'=. x
  cube=. y
  'B H W'=. $ cube
  'p HA WA'=. $ A

  NB. Reconstruct cube: X_hat = M A
  X_mat=. (B, H * W) $ cube
  A_mat=. (p, H * W) $ A
  X_hat=. M albedo_linalg_matmul A_mat

  diff=. X_mat - X_hat
  %:(+/ , *: diff) % (B * H * W)
)
