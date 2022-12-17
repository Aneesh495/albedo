NB. ==============================================================================
NB. Albedo - Non-Negative Least Squares (NNLS)
NB.
NB. Solves min ||A x - b||_2^2 subject to x >= 0 using the active set method
NB. (Lawson & Hanson, 1974). Guarantees exact convergence to global minimum.
NB. ==============================================================================

cocurrent 'base'

load 'core/linalg.ijs'

NB. ------------------------------------------------------------------------------
NB. Lawson-Hanson Active-Set NNLS Algorithm
NB. ------------------------------------------------------------------------------

NB. Solve min ||A x - b||^2 s.t. x >= 0
NB. A has shape (m, n), b has shape (m)
NB. Returns solution vector x of length n
solve_nnls=: 4 : 0
  A=. x
  b=. y
  'm n'=. $ A

  NB. Passive set P (indices of positive variables) and Active set Z (indices of zero variables)
  P=. 0 $ 0
  Z=. i. n
  x_vec=. n $ 0.0

  AtA=. (albedo_linalg_tr A) albedo_linalg_matmul A
  Atb=. (albedo_linalg_tr A) albedo_linalg_matmul b

  NB. Gradient vector w = At (b - A x)
  w=. Atb - (AtA albedo_linalg_matmul x_vec)

  max_iter=. 3 * n
  iter=. 0

  while. (0 < # Z) *. (iter < max_iter) do.
    iter=. iter + 1

    NB. Find index in Z with maximum positive gradient w_j
    w_Z=. Z { w
    max_w=. >./ w_Z
    if. max_w <: 1e_12 do. break. end.

    max_idx_in_Z=. (<: +/ w_Z = max_w) { (I. w_Z = max_w)
    t=. max_idx_in_Z { Z

    NB. Move t from Z to P
    P=. P , t
    Z=. Z -. t

    NB. Inner loop: unconstrained least squares on passive set P
    while. 1 do.
      AP=. P {"1 A
      AP_pinv=. albedo_linalg_pinv_svd AP
      zP=. AP_pinv albedo_linalg_matmul b

      NB. Check if all entries in zP are strictly positive
      if. *./ zP > 0.0 do.
        x_vec=. zP P} (n $ 0.0)
        break.
      end.

      NB. Find alpha = min { x_j / (x_j - z_j) } for z_j <= 0
      bad_indices=. I. zP <: 0.0
      alpha=. 1.0
      for_k. bad_indices do.
        idx_p=. k { P
        xk=. idx_p { x_vec
        zk=. k { zP
        cand_alpha=. xk % (xk - zk)
        if. (cand_alpha >: 0.0) *. (cand_alpha < alpha) do.
          alpha=. cand_alpha
        end.
      end.

      NB. Update x_vec = x_vec + alpha * (zP - xP)
      xP=. P { x_vec
      new_xP=. xP + alpha * (zP - xP)
      x_vec=. new_xP P} (n $ 0.0)

      NB. Move indices with x_j == 0 from P to Z
      zero_in_P=. P { I. x_vec <: 1e_12
      Z=. Z , zero_in_P
      P=. P -. zero_in_P
    end.

    NB. Update gradient w
    w=. Atb - (AtA albedo_linalg_matmul x_vec)
  end.

  0.0 >. x_vec
)

NB. Batch NNLS on data matrix Y of shape (m, N)
NB. Solves for each column of Y: A x = y_col
NB. Returns abundance matrix X of shape (n, N)
batch_nnls=: 4 : 0
  A=. x
  Y=. y
  'm N'=. $ Y
  'mA n'=. $ A
  X_res=. (n, N) $ 0.0

  for_col. i. N do.
    y_col=. col {"1 Y
    x_col=. A solve_nnls y_col
    X_res=. x_col (<a:; col)} X_res
  end.
  X_res
)


albedo_nnls_solve_nnls =: solve_nnls
