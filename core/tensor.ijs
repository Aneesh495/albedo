NB. ==============================================================================
NB. Albedo - Core Tensor Algebra Module
NB.
NB. Provides multi-way array manipulations, tensor matricization (mode-n unfolding),
NB. tensor folding, mode-n matrix products, Kronecker products, Khatri-Rao products,
NB. tensor contractions, Frobenius norms, and rank-1 outer products.
NB. ==============================================================================

cocurrent 'albedo_tensor'

NB. ------------------------------------------------------------------------------
NB. Dimension & Rank Helpers
NB. ------------------------------------------------------------------------------

NB. Tensor rank (order): number of modes/dimensions
tensor_order=: #@$

NB. Total number of scalar entries in tensor
tensor_size=: */@$

NB. Tensor Frobenius norm: sqrt(sum(x^2))
tensor_fnorm=: %:@(+/@:(,^:2)@:|)

NB. Tensor inner product between two tensors of matching shape
tensor_inner=: 4 : 0
  'shape_x shape_y'=. ($ x);($ y)
  assert. shape_x -: shape_y [ 'Tensor shapes must match for inner product'
  +/ (, x) * (, y)
)

NB. Rank-1 outer product of two vectors
vec_outer=: */

NB. Rank-1 outer product of three vectors: a o. b o. c
NB. Produces 3D tensor of shape (#a, #b, #c)
vec_outer3=: 3 : 0
  'a b c'=. y
  ab=. a */ b
  ab */ c
)

NB. ------------------------------------------------------------------------------
NB. Kronecker & Khatri-Rao Products
NB. ------------------------------------------------------------------------------

NB. Kronecker product of two 2D matrices: A (x) B
NB. If A is (m x n) and B is (p x q), result is (mp x nq)
kronecker=: 4 : 0
  'm n'=. $ x
  'p q'=. $ y
  prod=. x */ y
  NB. prod shape is (m, n, p, q) -> transpose to (m, p, n, q) -> reshape to (mp, nq)
  (m * p, n * q) $ 1 3 0 2 |: prod
)

NB. Khatri-Rao product (column-wise Kronecker product): A (*) B
NB. A has shape (I, R) and B has shape (J, R)
NB. Result has shape (IJ, R)
khatri_rao=: 4 : 0
  'I R1'=. $ x
  'J R2'=. $ y
  assert. R1 = R2 [ 'Column dimensions must match for Khatri-Rao product'
  R=. R1
  NB. Compute for each column r: A[:, r] (x) B[:, r]
  res=. (I * J, R) $ 0
  for_r. i. R do.
    col_a=. r {"1 x
    col_b=. r {"1 y
    col_kr=. , col_a */ col_b
    res=. col_kr (<a:; r)} res
  end.
  res
)

NB. ------------------------------------------------------------------------------
NB. Tensor Matricization (Unfolding) & Folding
NB. ------------------------------------------------------------------------------

NB. Mode-n Unfolding (Matricization): X_(n)
NB. Arranges mode-n fibers as columns of the resulting matrix.
NB. Dyadic: n tensor_unfold X
NB. n is 0-indexed mode (0 for mode-0, 1 for mode-1, etc.)
NB. Output matrix has shape: (dim_n, product of all other dims)
tensor_unfold=: 4 : 0
  n=. x
  X=. y
  dims=. $ X
  order=. # dims
  assert. (n >: 0) *. (n < order) [ 'Mode index out of range'

  NB. Permute axes so mode n is first, followed by all other modes in natural order
  other_modes=. (i. order) -. n
  perm=. n, other_modes
  permuted=. perm |: X

  NB. Reshape to 2D matrix
  dim_n=. n { dims
  other_prod=. */ (other_modes { dims)
  (dim_n, other_prod) $ permuted
)

NB. Tensor Folding: reconstruct original tensor from mode-n unfolded matrix
NB. Dyadic: (n; original_dims) tensor_fold unfolded_matrix
tensor_fold=: 4 : 0
  'n dims'=. x
  M=. y
  order=. # dims
  dim_n=. n { dims
  other_modes=. (i. order) -. n

  NB. Reshape M back to permuted dimensions
  perm_dims=. dim_n, (other_modes { dims)
  reshaped=. perm_dims $ M

  NB. Invert permutation
  perm=. n, other_modes
  inv_perm=. /: perm
  inv_perm |: reshaped
)

NB. ------------------------------------------------------------------------------
NB. Mode-n Matrix-Tensor Product: Y = X x_n U
NB. ------------------------------------------------------------------------------

NB. Multiplies tensor X along mode n with matrix U
NB. If X has dimensions (I_0, I_1, ..., I_n, ..., I_{N-1})
NB. and U has shape (J, I_n)
NB. Result has dimensions (I_0, I_1, ..., J, ..., I_{N-1})
NB. Dyadic: (n; U) mode_n_prod X
mode_n_prod=: 4 : 0
  'n U'=. x
  X=. y
  dims=. $ X
  order=. # dims
  'J In'=. $ U
  assert. In = (n { dims) [ 'Inner dimension mismatch for mode-n product'

  NB. Unfold X along mode n
  Xn=. n tensor_unfold X

  NB. Matrix multiply: U (J x In) with Xn (In x rest)
  Yn=. U +/ .* Xn

  NB. Determine new tensor shape
  new_dims=. J (n}) dims

  NB. Fold Yn back into tensor of new_dims
  (n; new_dims) tensor_fold Yn
)

NB. ------------------------------------------------------------------------------
NB. Multi-mode Tensor Product: Y = X x_0 U_0 x_1 U_1 ... x_{N-1} U_{N-1}
NB. ------------------------------------------------------------------------------

NB. Applies a list of factor matrices along their respective modes
NB. Dyadic: factor_matrices ttm X
ttm=: 4 : 0
  factors=. x
  Y=. y
  for_n. i. # factors do.
    U=. n {:: factors
    if. 0 < # U do.
      Y=. (n; U) mode_n_prod Y
    end.
  end.
  Y
)

NB. ------------------------------------------------------------------------------
NB. Tensor Slicing & Fiber Extraction
NB. ------------------------------------------------------------------------------

NB. Extract mode-n fiber at index tuple
NB. For 3D tensor of shape (B, H, W):
NB. mode-0 fiber is spectral vector at (h, w)
extract_spectral_vector=: 3 : 0
  'X h w'=. y
  X <a:; h; w
)

NB. Extract spatial band (2D image slice) at spectral band b
extract_spatial_band=: 3 : 0
  'X b'=. y
  b { X
)

NB. Reshape 3D hyperspectral cube (B x H x W) to 2D matrix (B x HW)
cube_to_spectral_matrix=: 3 : 0
  X=. y
  'B H W'=. $ X
  (B, H * W) $ X
)

NB. Reshape 2D spectral matrix (B x HW) back to 3D cube (B x H x W)
spectral_matrix_to_cube=: 4 : 0
  'H W'=. x
  M=. y
  B=. {. $ M
  (B, H, W) $ M
)
