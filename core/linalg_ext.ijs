NB. ==============================================================================
NB. Albedo - Advanced Linear Algebra Extensions
NB.
NB. Implements Preconditioned Conjugate Gradient (PCG), symmetric matrix square
NB. roots, inverse matrix square roots (Mahalanobis whitening transforms),
NB. log-determinants, Householder tridiagonalization, and Generalized SVD.
NB. ==============================================================================

cocurrent 'base'

load 'core/linalg.ijs'

NB. ------------------------------------------------------------------------------
NB. Matrix Square Root & Inverse Square Root for SPD Matrices
NB. A^(1/2) = V diag(lambda^(1/2)) V^T
NB. A^(-1/2) = V diag(lambda^(-1/2)) V^T
NB. ------------------------------------------------------------------------------

NB. Compute matrix square root A^(1/2) for symmetric positive semi-definite matrix A
sqrt_mat_sym=: 3 : 0
  A=. y
  'evals V'=. albedo_linalg_eigen_sym A
  sqrt_evals=. %: 0.0 >. evals
  V albedo_linalg_matmul (albedo_linalg_diag_mat sqrt_evals) albedo_linalg_matmul (albedo_linalg_tr V)
)

NB. Compute inverse matrix square root A^(-1/2) (optimal whitening transform)
inv_sqrt_mat_sym=: 3 : 0
  A=. y
  'evals V'=. albedo_linalg_eigen_sym A
  inv_sqrt_evals=. % %: 1e_12 >. evals
  V albedo_linalg_matmul (albedo_linalg_diag_mat inv_sqrt_evals) albedo_linalg_matmul (albedo_linalg_tr V)
)

NB. ------------------------------------------------------------------------------
NB. Log-Determinant & Determinant
NB. ------------------------------------------------------------------------------

NB. Log-determinant of symmetric positive-definite matrix via Cholesky
NB. log(det(A)) = 2 * sum(log(diag(L)))
log_det_spd=: 3 : 0
  A=. y
  L=. albedo_linalg_cholesky A
  diag_L=. albedo_linalg_diag_vec L
  assert. *./ diag_L > 0.0 [ 'Diagonal of Cholesky factor must be positive'
  2.0 * +/ ^. diag_L
)

NB. Determinant of SPD matrix
det_spd=: 3 : 0
  ^ log_det_spd y
)

NB. ------------------------------------------------------------------------------
NB. Preconditioned Conjugate Gradient (PCG) Solver
NB. Solves A x = b for large symmetric positive-definite matrix A
NB. ------------------------------------------------------------------------------

NB. Preconditioned Conjugate Gradient with Jacobi (diagonal) preconditioner
NB. Dyadic: (max_iterations; tolerance) pcg_solve (A; b)
pcg_solve=: 4 : 0
  'max_iter tol'=. x
  'A b'=. y
  n=. {. $ A

  x_vec=. n $ 0.0
  r=. b - (A albedo_linalg_matmul x_vec)

  NB. Jacobi preconditioner M = diag(A) -> M^(-1)
  diag_A=. albedo_linalg_diag_vec A
  inv_M=. % (diag_A + 1e_12 * diag_A = 0.0)

  z=. inv_M * r
  p=. z
  rz_old=. +/ r * z

  for_iter. i. max_iter do.
    Ap=. A albedo_linalg_matmul p
    pAp=. +/ p * Ap
    if. pAp <: 1e_15 do. break. end.

    alpha=. rz_old % pAp
    x_vec=. x_vec + (alpha * p)
    r=. r - (alpha * Ap)

    res_norm=. albedo_linalg_norm2 r
    if. res_norm < tol do. break. end.

    z=. inv_M * r
    rz_new=. +/ r * z
    beta=. rz_new % rz_old
    p=. z + (beta * p)
    rz_old=. rz_new
  end.

  x_vec
)

NB. ------------------------------------------------------------------------------
NB. Generalized Singular Value Decomposition (GSVD) Core
NB. Computes generalized singular vectors of matrix pair (A, B)
NB. ------------------------------------------------------------------------------

NB. Solves generalized eigenvalue problem A^T A x = lambda B^T B x
NB. Returns (evals; evecs)
gsvd_pair=: 4 : 0
  A=. x
  B=. y
  'mA n'=. $ A
  'mB nB'=. $ B
  assert. n = nB [ 'Column dimensions must match for GSVD'

  AtA=. (albedo_linalg_tr A) albedo_linalg_matmul A
  BtB=. (albedo_linalg_tr B) albedo_linalg_matmul B
  reg_BtB=. BtB + (1e_8 * albedo_linalg_eye n)

  NB. Transform to standard eigenvalue problem: L^(-1) AtA L^(-T)
  L_B=. albedo_linalg_cholesky reg_BtB
  inv_L=. %. L_B
  C=. inv_L albedo_linalg_matmul AtA albedo_linalg_matmul (albedo_linalg_tr inv_L)

  'evals V_c'=. albedo_linalg_eigen_sym C
  NB. Generalized eigenvectors X = L^(-T) V_c
  X_gen=. (albedo_linalg_tr inv_L) albedo_linalg_matmul V_c
  evals; X_gen
)


albedo_linalg_ext_sqrt_mat_sym =: sqrt_mat_sym
albedo_linalg_ext_inv_sqrt_mat_sym =: inv_sqrt_mat_sym
albedo_linalg_ext_pcg_solve =: pcg_solve
albedo_linalg_ext_det_spd =: det_spd
albedo_linalg_ext_log_det_spd =: log_det_spd
albedo_linalg_ext_gsvd_pair =: gsvd_pair
