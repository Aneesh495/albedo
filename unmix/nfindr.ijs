NB. ==============================================================================
NB. Albedo - N-FINDR Simplex Volume Maximization Endmember Extraction
NB.
NB. Implements N-FINDR (Winter, 1999) to find the simplex of maximum volume
NB. inscribed within the hyperspectral data cloud. The vertices of the maximum
NB. volume simplex correspond to pure spectral endmembers.
NB. ==============================================================================

cocurrent 'albedo_nfindr'

load 'core/linalg.ijs'

NB. ------------------------------------------------------------------------------
NB. Simplex Volume Calculation
NB. V(E) = (1 / (p - 1)!) * |det [ 1^T ; E ]|
NB. ------------------------------------------------------------------------------

simplex_volume=: 3 : 0
  E_red=. y
  'p_minus_1 p'=. $ E_red
  NB. Augment top row of ones: shape (p, p)
  ones_row=. 1, p $ 1.0
  Aug=. (1.0 $~ p) ,: E_red
  NB. Compute determinant of Aug via QR factorization
  'Q R'=. albedo_linalg_qr_mgs Aug
  diag_R=. albedo_linalg_diag_vec R
  det_val=. */ diag_R
  fact=. ! (p - 1)
  (| det_val) % (fact >. 1.0)
)

NB. ------------------------------------------------------------------------------
NB. N-FINDR Iterative Algorithm
NB. ------------------------------------------------------------------------------

NB. Extract p endmembers from spectral data matrix X of shape (B, N)
NB. Dyadic: p extract_nfindr X
extract_nfindr=: 4 : 0
  p=. x
  X=. y
  'B N'=. $ X
  assert. p <: B [ 'p must not exceed B'
  assert. p <: N [ 'p must not exceed N'

  NB. SVD dimension reduction to (p - 1) dimensions
  'svals U V'=. albedo_linalg_svd_power X
  Up=. (p - 1) {"1 U
  NB. Project data onto (p - 1)-dimensional space: Y shape is ((p-1) x N)
  Y=. (albedo_linalg_tr Up) albedo_linalg_matmul X

  NB. Initialize p random endmember indices
  perm_idx=. /: ? N $ 1000
  endmember_indices=. p {. perm_idx
  E_red=. endmember_indices {"1 Y

  current_vol=. simplex_volume E_red
  max_iterations=. 20

  for_iter. i. max_iterations do.
    vol_changed=. 0

    for_j. i. p do.
      for_i. i. N do.
        NB. Candidate replacement
        cand_idx=. i
        if. cand_idx e. endmember_indices do. continue. end.

        NB. Test volume with pixel i replacing endmember j
        test_E=. (cand_idx {"1 Y) (<a:; j)} E_red
        test_vol=. simplex_volume test_E

        if. test_vol > current_vol do.
          current_vol=. test_vol
          E_red=. test_E
          endmember_indices=. cand_idx (<j)} endmember_indices
          vol_changed=. 1
        end.
      end.
    end.

    if. 0 = vol_changed do. break. end.
  end.

  M_endmembers=. endmember_indices {"1 X
  M_endmembers; endmember_indices; current_vol
)
