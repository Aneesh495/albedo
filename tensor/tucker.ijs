NB. ==============================================================================
NB. Albedo - Tucker3 Higher-Order SVD (HOSVD) Tensor Decomposition
NB.
NB. Decomposes 3-way tensor X into core tensor G and orthonormal factor
NB. matrices U0, U1, U2:
NB. X ~ G x_0 U0 x_1 U1 x_2 U2
NB. Implements De Lathauwer et al. (SIAM J. Matrix Anal. Appl., 2000).
NB. ==============================================================================

cocurrent 'albedo_tucker'

load 'core/tensor.ijs'
load 'core/linalg.ijs'

NB. ------------------------------------------------------------------------------
NB. Higher-Order SVD (HOSVD)
NB. ------------------------------------------------------------------------------

NB. Compute HOSVD decomposition of 3D tensor X to ranks (r0, r1, r2)
NB. Dyadic: (r0; r1; r2) fit_hosvd X
fit_hosvd=: 4 : 0
  'r0 r1 r2'=. x
  X=. y
  'I J K'=. $ X

  NB. Mode-0 unfolding SVD
  X0=. 0 albedo_tensor_tensor_unfold X
  'svals0 U0_full V0'=. albedo_linalg_svd_power X0
  U0=. r0 {"1 U0_full

  NB. Mode-1 unfolding SVD
  X1=. 1 albedo_tensor_tensor_unfold X
  'svals1 U1_full V1'=. albedo_linalg_svd_power X1
  U1=. r1 {"1 U1_full

  NB. Mode-2 unfolding SVD
  X2=. 2 albedo_tensor_tensor_unfold X
  'svals2 U2_full V2'=. albedo_linalg_svd_power X2
  U2=. r2 {"1 U2_full

  NB. Compute core tensor G = X x_0 U0^T x_1 U1^T x_2 U2^T
  U0t=. albedo_linalg_tr U0
  U1t=. albedo_linalg_tr U1
  U2t=. albedo_linalg_tr U2

  G=. (U0t; U1t; U2t) albedo_tensor_ttm X

  G; U0; U1; U2
)

NB. ------------------------------------------------------------------------------
NB. Reconstruction of Tensor from Tucker Model
NB. X_hat = G x_0 U0 x_1 U1 x_2 U2
NB. ------------------------------------------------------------------------------

reconstruct_tucker=: 3 : 0
  'G U0 U1 U2'=. y
  (U0; U1; U2) albedo_tensor_ttm G
)
