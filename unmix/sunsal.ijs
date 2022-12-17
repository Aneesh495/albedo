NB. ==============================================================================
NB. Albedo - Sparse Unmixing via Augmented Lagrangian (SUNSAL)
NB.
NB. Implements Sparse Unmixing by Variable Splitting and Augmented Lagrangian
NB. (Bioucas-Dias & Figueiredo, IEEE TGRS 2010). Solves:
NB. min (1/2) ||M s - y||_2^2 + lambda ||s||_1  subject to  s >= 0
NB. using Alternating Direction Method of Multipliers (ADMM).
NB. ==============================================================================

cocurrent 'base'

load 'core/linalg.ijs'

NB. Soft thresholding operator: S_tau(v) = sign(v) * max(0, |v| - tau)
soft_thresh=: 4 : 0
  tau=. x
  v=. y
  sign=. (v > 0.0) - (v < 0.0)
  sign * (0.0 >. (| v) - tau)
)

NB. Solve sparse unmixing for single pixel
NB. Dyadic: (lambda; M) solve_sunsal y_pixel
solve_sunsal=: 4 : 0
  'lambda M'=. x
  y_pixel=. y
  'B p'=. $ M

  mu_admm=. 1e_2
  max_iter=. 50
  tol=. 1e_4

  MtM=. (albedo_linalg_tr M) albedo_linalg_matmul M
  Mty=. (albedo_linalg_tr M) albedo_linalg_matmul y_pixel
  inv_matrix=. %. (MtM + (mu_admm * albedo_linalg_eye p))

  NB. Initialize variables: s (primal), v (auxiliary), d (Lagrange multiplier)
  s_vec=. p $ 0.0
  v_vec=. p $ 0.0
  d_vec=. p $ 0.0

  for_iter. i. max_iter do.
    NB. Step 1: Update s = (M^T M + mu I)^(-1) (M^T y + mu (v - d))
    rhs=. Mty + (mu_admm * (v_vec - d_vec))
    s_new=. inv_matrix albedo_linalg_matmul rhs

    NB. Step 2: Update v with non-negativity and L1 thresholding
    NB. v = max(0, soft_thresh(s_new + d, lambda / mu))
    tau=. lambda % mu_admm
    arg=. s_new + d_vec
    v_new=. 0.0 >. (tau soft_thresh arg)

    NB. Step 3: Update Lagrange multiplier d = d + (s - v)
    d_vec=. d_vec + (s_new - v_new)

    NB. Check primal residual convergence: ||s - v||_2 < tol
    prim_res=. albedo_linalg_norm2 (s_new - v_new)
    s_vec=. s_new
    v_vec=. v_new
    if. prim_res < tol do. break. end.
  end.

  0.0 >. v_vec
)


albedo_sunsal_solve_sunsal =: solve_sunsal
