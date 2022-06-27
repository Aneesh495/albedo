NB. ==============================================================================
NB. Albedo - Canonical Polyadic (CP/PARAFAC) Tensor Decomposition via ALS
NB.
NB. Decomposes 3-way tensor X into R rank-1 components using Alternating Least
NB. Squares (ALS): X ~ sum_{r=1}^R lambda_r (a_r o. b_r o. c_r)
NB. Computes factor matrices A, B, C and component weights lambda.
NB. ==============================================================================

cocurrent 'albedo_cp_als'

load 'core/tensor.ijs'
load 'core/linalg.ijs'

NB. ------------------------------------------------------------------------------
NB. Alternating Least Squares (ALS) Solver
NB. ------------------------------------------------------------------------------

NB. Fit rank-R CP decomposition on 3D tensor X of shape (I, J, K)
NB. Dyadic: (R; max_iterations; tolerance) fit_cp_als X
fit_cp_als=: 4 : 0
  'R max_iter tol'=. x
  X=. y
  'I J K'=. $ X

  NB. Initialize factor matrices A, B, C with random values
  A=. ? I, R $ 1000
  B=. ? J, R $ 1000
  C=. ? K, R $ 1000

  NB. Normalize columns
  A=. A % ((albedo_linalg_norm2 "1 albedo_linalg_tr A) >. 1e_8)
  B=. B % ((albedo_linalg_norm2 "1 albedo_linalg_tr B) >. 1e_8)
  C=. C % ((albedo_linalg_norm2 "1 albedo_linalg_tr C) >. 1e_8)

  X0=. 0 albedo_tensor_tensor_unfold X
  X1=. 1 albedo_tensor_tensor_unfold X
  X2=. 2 albedo_tensor_tensor_unfold X

  lambda_weights=. R $ 1.0
  norm_X=. albedo_tensor_tensor_fnorm X

  for_iter. i. max_iter do.
    NB. --- Update A ---
    NB. V = (C^T C) * (B^T B) (Hadamard product)
    CtC=. (albedo_linalg_tr C) albedo_linalg_matmul C
    BtB=. (albedo_linalg_tr B) albedo_linalg_matmul B
    V_A=. CtC * BtB
    inv_VA=. %. (V_A + (1e_8 * albedo_linalg_eye R))
    NB. Khatri-Rao product C (*) B of shape (KJ, R)
    KR_CB=. C albedo_tensor_khatri_rao B
    A_raw=. (X0 albedo_linalg_matmul KR_CB) albedo_linalg_matmul inv_VA
    NB. Normalize columns of A
    col_norms_A=. albedo_linalg_norm2 "1 albedo_linalg_tr A_raw
    col_norms_A=. col_norms_A >. 1e_12
    A=. A_raw % col_norms_A

    NB. --- Update B ---
    AtA=. (albedo_linalg_tr A) albedo_linalg_matmul A
    V_B=. CtC * AtA
    inv_VB=. %. (V_B + (1e_8 * albedo_linalg_eye R))
    KR_CA=. C albedo_tensor_khatri_rao A
    B_raw=. (X1 albedo_linalg_matmul KR_CA) albedo_linalg_matmul inv_VB
    col_norms_B=. albedo_linalg_norm2 "1 albedo_linalg_tr B_raw
    col_norms_B=. col_norms_B >. 1e_12
    B=. B_raw % col_norms_B

    NB. --- Update C ---
    V_C=. BtB * AtA
    inv_VC=. %. (V_C + (1e_8 * albedo_linalg_eye R))
    KR_BA=. B albedo_tensor_khatri_rao A
    C_raw=. (X2 albedo_linalg_matmul KR_BA) albedo_linalg_matmul inv_VC
    col_norms_C=. albedo_linalg_norm2 "1 albedo_linalg_tr C_raw
    col_norms_C=. col_norms_C >. 1e_12
    C=. C_raw % col_norms_C

    NB. Component weights lambda = col_norms_C
    lambda_weights=. col_norms_C
  end.

  lambda_weights; A; B; C
)

NB. ------------------------------------------------------------------------------
NB. Reconstruction of Tensor from CP Factors
NB. X_hat = sum_{r=1}^R lambda_r (a_r o. b_r o. c_r)
NB. ------------------------------------------------------------------------------

reconstruct_cp=: 3 : 0
  'lambda_weights A B C'=. y
  R=. # lambda_weights
  'I R_A'=. $ A
  'J R_B'=. $ B
  'K R_C'=. $ C

  X_hat=. (I, J, K) $ 0.0

  for_r. i. R do.
    w_r=. r { lambda_weights
    ar=. r {"1 A
    br=. r {"1 B
    cr=. r {"1 C
    rank1=. (ar; br; cr) albedo_tensor_vec_outer3 ''
    X_hat=. X_hat + (w_r * rank1)
  end.

  X_hat
)
