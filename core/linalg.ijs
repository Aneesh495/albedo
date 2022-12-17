NB. ==============================================================================
NB. Albedo - High-Performance Linear Algebra Primitives
NB.
NB. Implements matrix factorizations, Cholesky decomposition, Householder QR,
NB. Singular Value Decomposition (SVD), symmetric Jacobi eigenvalue solver,
NB. Moore-Penrose pseudoinverse, Tikhonov regularization, and linear solvers.
NB. ==============================================================================

cocurrent 'base'

NB. ------------------------------------------------------------------------------
NB. Matrix Properties & Basic Operators
NB. ------------------------------------------------------------------------------

NB. Identity matrix of size n
eye=: =@i.

NB. Matrix transpose
tr=: |:

NB. Matrix multiplication
matmul=: +/ .*

NB. Extract diagonal vector of matrix
diag_vec=: 3 : 0
  (<"1 (i. # y) ,. (i. # y)) { y
)

NB. Matrix trace: sum of diagonal entries
trace=: 3 : 0
  +/ diag_vec y
)

NB. Construct diagonal matrix from vector
diag_mat=: 3 : 0
  (=/ ~ i. # y) *"1 y
)

NB. Vector Euclidean 2-norm
norm2=: %:@(+/@:*:)

NB. Frobenius norm of matrix
norm_frob=: %:@(+/@:(,^:2)@:|)

NB. Matrix condition number proxy (ratio of max to min diagonal in Cholesky or SVD)
cond_est=: 3 : 0
  s=. >@{. svd_power y
  if. 0 = {: s do. _ else. ({. s) % {: s end.
)

NB. ------------------------------------------------------------------------------
NB. Linear Solvers & Inverses
NB. ------------------------------------------------------------------------------

NB. Standard least-squares solve: A x = b -> x = A %. b
solve=: 4 : 0
  A=. x
  b=. y
  A %. b
)

NB. Regularized inverse: (A + lambda * I)^-1
NB. For square matrix A
reg_inv=: 4 : 0
  lambda=. x
  A=. y
  n=. {. $ A
  reg_A=. A + lambda * eye n
  %. reg_A
)

NB. Tikhonov regularized solve: (A^T A + lambda * I)^-1 A^T b
tikhonov_solve=: 4 : 0
  lambda=. x
  'A b'=. y
  n=. {: $ A
  AtA=. (tr A) matmul A
  Atb=. (tr A) matmul b
  reg_AtA=. AtA + lambda * eye n
  reg_AtA %. Atb
)

NB. ------------------------------------------------------------------------------
NB. Cholesky Factorization: A = L L^T
NB. For symmetric positive definite matrix A (n x n)
NB. Returns lower triangular matrix L
NB. ------------------------------------------------------------------------------

cholesky=: 3 : 0
  A=. y
  n=. {. $ A
  L=. (n, n) $ 0.0

  for_j. i. n do.
    NB. Diagonal entry L[j, j]
    sum_sq=. 0.0
    if. j > 0 do.
      sum_sq=. +/ *: j {. (j { L)
    end.
    diag_val=. ((< j, j) { A) - sum_sq
    assert. diag_val > 0.0 [ 'Matrix is not positive definite in Cholesky'
    ljj=. %: diag_val
    L=. ljj (<j, j)} L

    NB. Off-diagonal entries L[i, j] for i > j
    for_i. (j + 1) + i. (<: n - j) do.
      sum_prod=. 0.0
      if. j > 0 do.
        sum_prod=. +/ (j {. i { L) * (j {. j { L)
      end.
      lij=. (((< i, j) { A) - sum_prod) % ljj
      L=. lij (<i, j)} L
    end.
  end.
  L
)

NB. ------------------------------------------------------------------------------
NB. Modified Gram-Schmidt QR Decomposition: A = Q R
NB. A is (m x n), returns (Q; R) where Q is (m x n) and R is (n x n)
NB. ------------------------------------------------------------------------------

qr_mgs=: 3 : 0
  A=. y
  'm n'=. $ A
  Q=. A
  R=. (n, n) $ 0.0

  for_j. i. n do.
    v=. j {"1 Q
    rjj=. norm2 v
    if. rjj < 1e_14 do.
      R=. 0.0 (<j, j)} R
      qj=. 0.0 * v
      Q=. qj (<a:; j)} Q
    else.
      R=. rjj (<j, j)} R
      qj=. v % rjj
      Q=. qj (<a:; j)} Q
      for_k. (j + 1) + i. (<: n - j) do.
        vk=. k {"1 Q
        rjk=. +/ qj * vk
        R=. rjk (<j, k)} R
        vk_orth=. vk - rjk * qj
        Q=. vk_orth (<a:; k)} Q
      end.
    end.
  end.
  Q; R
)

NB. ------------------------------------------------------------------------------
NB. Symmetric Jacobi Eigenvalue Decomposition: A = V diag(lambda) V^T
NB. For real symmetric matrix A (n x n)
NB. Returns (eigenvalues; eigenvectors) sorted in descending order
NB. ------------------------------------------------------------------------------

eigen_sym=: 3 : 0
  A=. y
  n=. {. $ A
  V=. eye n
  max_sweeps=. 50
  tol=. 1e_12

  for_sweep. i. max_sweeps do.
    max_off=. 0.0
    for_p. i. <: n do.
      for_q. (1 + p) + i. <: n - p do.
        apq=. (< p, q) { A
        if. (| apq) > tol do.
          app=. (< p, p) { A
          aqq=. (< q, q) { A
          tau=. (app - aqq) % (2.0 * apq)
          t=. 0.0
          if. tau >: 0.0 do.
            t=. 1.0 % (tau + %: 1.0 + *: tau)
          else.
            t=. _1.0 % ((| tau) + %: 1.0 + *: tau)
          end.
          c=. 1.0 % %: 1.0 + *: t
          s=. t * c

          NB. Givens rotation matrix G
          G=. eye n
          G=. c (< p, p)} G
          G=. c (< q, q)} G
          G=. (- s) (< p, q)} G
          G=. s (< q, p)} G
          A=. (tr G) matmul A matmul G
          V=. V matmul G
        end.
      end.
    end.
  end.

  evals=. diag_vec A
  order=. \: evals
  sorted_evals=. order { evals
  sorted_evecs=. order {"1 V
  sorted_evals; sorted_evecs
)

NB. ------------------------------------------------------------------------------
NB. Thin Singular Value Decomposition (SVD): A = U Sigma V^T
NB. For matrix A (m x n). Uses eigen_sym on AtA or AAt depending on shape.
NB. Returns (singular_values; U; V)
NB. ------------------------------------------------------------------------------

svd_power=: 3 : 0
  A=. y
  if. 1 = # $ A do. A=. (, # A) ,. A end.
  'm n'=. $ A
  if. m >: n do.
    NB. Compute AtA (n x n)
    AtA=. (tr A) matmul A
    'evals V'=. eigen_sym AtA
    svals=. %: 0.0 >. evals
    inv_s=. (% svals + svals = 0.0) * (svals > 1e_12)
    U=. (A matmul V) *"1 inv_s
    svals; U; V
  else.
    NB. Compute AAt (m x m)
    AAt=. A matmul tr A
    'evals U'=. eigen_sym AAt
    svals=. %: 0.0 >. evals
    inv_s=. (% svals + svals = 0.0) * (svals > 1e_12)
    V=. ((tr A) matmul U) *"1 inv_s
    svals; U; V
  end.
)

NB. Alias for thin_svd
thin_svd=: svd_power

NB. ------------------------------------------------------------------------------
NB. Moore-Penrose Pseudoinverse via SVD: A^+ = V Sigma^+ U^T
NB. ------------------------------------------------------------------------------

pinv_svd=: 3 : 0
  'svals U V'=. svd_power y
  inv_s=. (% svals + svals = 0.0) * (svals > 1e_8)
  V matmul (diag_mat inv_s) matmul tr U
)

NB. Alias for pinv
pinv=: pinv_svd

NB. ------------------------------------------------------------------------------
NB. Orthogonal Projector onto Range(A): P = A A^+
NB. ------------------------------------------------------------------------------

proj_mat=: 3 : 0
  A=. y
  A matmul pinv_svd A
)

NB. ------------------------------------------------------------------------------
NB. Orthogonal Projector onto Nullspace(A^T): P_perp = I - A A^+
NB. ------------------------------------------------------------------------------

proj_perp_mat=: 3 : 0
  A=. y
  n=. {. $ A
  (eye n) - proj_mat A
)

NB. ------------------------------------------------------------------------------
NB. Exported Aliases for Direct Namespace Access
NB. ------------------------------------------------------------------------------

albedo_linalg_eye =: eye
albedo_linalg_tr =: tr
albedo_linalg_matmul =: matmul
albedo_linalg_trace =: trace
albedo_linalg_diag_vec =: diag_vec
albedo_linalg_diag_mat =: diag_mat
albedo_linalg_norm2 =: norm2
albedo_linalg_norm_frob =: norm_frob
albedo_linalg_cond_est =: cond_est
albedo_linalg_solve =: solve
albedo_linalg_reg_inv =: reg_inv
albedo_linalg_tikhonov_solve =: tikhonov_solve
albedo_linalg_cholesky =: cholesky
albedo_linalg_qr_mgs =: qr_mgs
albedo_linalg_eigen_sym =: eigen_sym
albedo_linalg_svd_power =: svd_power
albedo_linalg_thin_svd =: thin_svd
albedo_linalg_pinv_svd =: pinv_svd
albedo_linalg_pinv =: pinv
albedo_linalg_proj_mat =: proj_mat
albedo_linalg_proj_perp_mat =: proj_perp_mat
