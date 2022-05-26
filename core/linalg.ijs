NB. ==============================================================================
NB. Albedo - High-Performance Linear Algebra Primitives
NB.
NB. Implements matrix factorizations, Cholesky decomposition, Householder QR,
NB. Singular Value Decomposition (SVD), symmetric Jacobi eigenvalue solver,
NB. Moore-Penrose pseudoinverse, Tikhonov regularization, and linear solvers.
NB. ==============================================================================

cocurrent 'albedo_linalg'

NB. ------------------------------------------------------------------------------
NB. Matrix Properties & Basic Operators
NB. ------------------------------------------------------------------------------

NB. Identity matrix of size n
eye=: =@i.

NB. Matrix transpose
tr=: |:

NB. Matrix multiplication
matmul=: +/ .*

NB. Matrix trace: sum of diagonal entries
trace=: +/@:(<0 1&|:)

NB. Extract diagonal vector of matrix
diag_vec=: (<0 1&|:)

NB. Construct diagonal matrix from vector
diag_mat=: =@i.@# * ]

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

NB. Regularized inverse: (A + lambda * I)^(-1)
reg_inv=: 4 : 0
  lambda=. x
  A=. y
  n=. {. $ A
  reg_A=. A + lambda * eye n
  %. reg_A
)

NB. Tikhonov regularized least-squares: min ||Ax - b||^2 + lambda ||x||^2
NB. Solution: (A^T A + lambda I)^(-1) A^T b
tikhonov_solve=: 4 : 0
  'lambda A'=. x
  b=. y
  'm n'=. $ A
  AtA=. (tr A) matmul A
  Atb=. (tr A) matmul b
  reg_AtA=. AtA + lambda * eye n
  reg_AtA %. Atb
)

NB. ------------------------------------------------------------------------------
NB. Cholesky Decomposition: A = L L^T
NB. For symmetric positive-definite matrix A (n x n)
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
      prev_row=. (j { L) {.@:(i.&j)
      sum_sq=. +/ *: j {. (j { L)
    end.
    diag_val=. (A < j, j) - sum_sq
    assert. diag_val > 0.0 [ 'Matrix is not positive definite in Cholesky'
    ljj=. %: diag_val
    L=. ljj (<j, j)} L

    NB. Off-diagonal entries L[i, j] for i > j
    for_i. (j + 1) + i. (n - j - 1) do.
      sum_prod=. 0.0
      if. j > 0 do.
        sum_prod=. +/ (j {. i { L) * (j {. j { L)
      end.
      lij=. ((A < i, j) - sum_prod) % ljj
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
    assert. rjj > 1e_14 [ 'Linearly dependent column in QR decomposition'
    R=. rjj (<j, j)} R
    qj=. v % rjj
    Q=. qj (<a:; j)} Q

    NB. Project subsequent columns onto qj and orthogonalize
    for_k. (j + 1) + i. (n - j - 1) do.
      vk=. k {"1 Q
      rjk=. +/ qj * vk
      R=. rjk (<j, k)} R
      vk_orth=. vk - rjk * qj
      Q=. vk_orth (<a:; k)} Q
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
    NB. Find max off-diagonal element
    max_off=. 0.0
    for_p. i. n - 1 do.
      for_q. (p + 1) + i. (n - p - 1) do.
        val=. | A < p, q
        if. val > max_off do. max_off=. val end.
      end.
    end.

    if. max_off < tol do. break. end.

    for_p. i. n - 1 do.
      for_q. (p + 1) + i. (n - p - 1) do.
        apq=. A < p, q
        if. (| apq) > (tol * 0.01) do.
          app=. A < p, p
          aqq=. A < q, q
          tau=. (aqq - app) % (2.0 * apq)
          t=. 0.0
          if. tau >: 0.0 do.
            t=. 1.0 % (tau + %: 1.0 + *: tau)
          else.
            t=. _1.0 % ((| tau) + %: 1.0 + *: tau)
          end.
          c=. 1.0 % %: 1.0 + *: t
          s=. t * c

          NB. Update rows and columns p and q of A
          Ap=. p { A
          Aq=. q { A
          new_Ap=. (c * Ap) - (s * Aq)
          new_Aq=. (s * Ap) + (c * Aq)
          A=. new_Ap (p}) A
          A=. new_Aq (q}) A

          NB. Update columns of A
          Acp=. p {"1 A
          Acq=. q {"1 A
          new_Acp=. (c * Acp) - (s * Acq)
          new_Acq=. (s * Acp) + (c * Acq)
          A=. new_Acp (<a:; p)} A
          A=. new_Acq (<a:; q)} A

          NB. Fix exact diagonal and zero out (p, q)
          app_new=. (app * *: c) - (2.0 * s * c * apq) + (aqq * *: s)
          aqq_new=. (app * *: s) + (2.0 * s * c * apq) + (aqq * *: c)
          A=. app_new (<p, p)} A
          A=. aqq_new (<q, q)} A
          A=. 0.0 (<p, q)} A
          A=. 0.0 (<q, p)} A

          NB. Update eigenvector matrix V
          Vp=. p {"1 V
          Vq=. q {"1 V
          new_Vp=. (c * Vp) - (s * Vq)
          new_Vq=. (s * Vp) + (c * Vq)
          V=. new_Vp (<a:; p)} V
          V=. new_Vq (<a:; q)} V
        end.
      end.
    end.
  end.

  evals=. diag_vec A
  NB. Sort eigenvalues in descending order
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
  'm n'=. $ A
  if. m >: n do.
    NB. Compute AtA (n x n)
    AtA=. (tr A) matmul A
    'evals V'=. eigen_sym AtA
    NB. Singular values are sqrt(max(0, evals))
    svals=. %: 0.0 >. evals
    NB. Compute U = A V diag(1/s)
    inv_s=. 0.0 (I. svals = 0.0)} (% svals + 1e_15 * svals = 0.0)
    U=. (A matmul V) * inv_s
    svals; U; V
  else.
    NB. Compute AAt (m x m)
    AAt=. A matmul tr A
    'evals U'=. eigen_sym AAt
    svals=. %: 0.0 >. evals
    inv_s=. 0.0 (I. svals = 0.0)} (% svals + 1e_15 * svals = 0.0)
    V=. ((tr A) matmul U) * inv_s
    svals; U; V
  end.
)

NB. ------------------------------------------------------------------------------
NB. Moore-Penrose Pseudoinverse via SVD
NB. A^+ = V diag(1/s) U^T
NB. ------------------------------------------------------------------------------

pinv_svd=: 3 : 0
  A=. y
  'svals U V'=. svd_power A
  tol=. 1e_12 * {. svals
  inv_s=. svals * 0.0
  pos_idx=. I. svals > tol
  inv_s=. (% pos_idx { svals) pos_idx} inv_s
  V matmul (diag_mat inv_s) matmul tr U
)

NB. ------------------------------------------------------------------------------
NB. Orthogonal Projection Matrix
NB. Projects vectors onto column space of A: P = A (A^T A)^(-1) A^T
NB. Orthogonal projector (nullspace projector): P_perp = I - P
NB. ------------------------------------------------------------------------------

proj_mat=: 3 : 0
  A=. y
  m=. {. $ A
  P=. A matmul (pinv_svd A)
  P
)

proj_perp_mat=: 3 : 0
  A=. y
  m=. {. $ A
  (eye m) - proj_mat A
)
