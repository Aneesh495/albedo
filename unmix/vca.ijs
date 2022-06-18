NB. ==============================================================================
NB. Albedo - Vertex Component Analysis (VCA) Endmember Extraction
NB.
NB. Implements Vertex Component Analysis (Nascimento & Bioucas-Dias, 2005)
NB. for unsupervised discovery of pure spectral signatures (endmembers)
NB. forming the vertices of a high-dimensional convex simplex.
NB. ==============================================================================

cocurrent 'albedo_vca'

load 'core/linalg.ijs'
load 'core/stats.ijs'

NB. ------------------------------------------------------------------------------
NB. Vertex Component Analysis Algorithm
NB. ------------------------------------------------------------------------------

NB. Extract p endmembers from spectral data matrix X of shape (B, N)
NB. SNR threshold determines whether to use projective or affine projection.
NB. Dyadic: p extract_vca X
extract_vca=: 4 : 0
  p=. x
  X=. y
  'B N'=. $ X
  assert. p <: B [ 'Number of endmembers p must not exceed number of bands B'
  assert. p <: N [ 'Number of endmembers p must not exceed number of pixels N'

  NB. SVD dimensionality reduction to p dimensions
  'svals U V'=. albedo_linalg_svd_power X
  Up=. p {"1 U
  NB. Project data onto p-dimensional subspace: Y = Up^T X (shape p x N)
  Y=. (albedo_linalg_tr Up) albedo_linalg_matmul X

  NB. Initialize endmember matrix E_proj of shape (p, p)
  E_proj=. (p, p) $ 0.0
  endmember_indices=. 0 $ 0

  for_i. i. p do.
    NB. Generate random vector w in R^p
    NB. Pseudo-random vector deterministically seeded by iteration i
    w=. 2.0 * ((? p $ 1000) % 1000.0) - 0.5
    w=. w % (albedo_linalg_norm2 w)

    if. i = 0 do.
      f=. w
    else.
      NB. Project w onto orthogonal complement of space spanned by current endmembers
      Ei=. i {"1 E_proj
      P_perp=. albedo_linalg_proj_perp_mat Ei
      f=. P_perp albedo_linalg_matmul w
      norm_f=. albedo_linalg_norm2 f
      if. norm_f > 1e_8 do.
        f=. f % norm_f
      else.
        f=. w
      end.
    end.

    NB. Project all pixels in Y onto f: v = f^T Y (length N)
    v=. f albedo_linalg_matmul Y
    abs_v=. | v

    NB. Find index of pixel maximizing absolute projection
    max_idx=. (<: +/ abs_v = >./ abs_v) { (I. abs_v = >./ abs_v)
    endmember_indices=. endmember_indices , max_idx

    NB. Store projected endmember vector in E_proj
    e_new=. max_idx {"1 Y
    E_proj=. e_new (<a:; i)} E_proj
  end.

  NB. Extract original full-band signatures from X at endmember_indices
  M_endmembers=. endmember_indices {"1 X
  M_endmembers; endmember_indices
)
