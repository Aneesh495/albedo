NB. ==============================================================================
NB. ALBEDO - High-Dimensional Hyperspectral Imaging Engine
NB. unmix/nmf.ijs : Non-negative Matrix Factorization (NMF) for Blind Unmixing
NB. Implements Lee-Seung multiplicative updates with sum-to-one projection
NB. and Frobenius reconstruction loss for simultaneous endmember and abundance retrieval.
NB. ==============================================================================

cocurrent 'albedo_nmf'

load 'core/linalg.ijs'

NB. Initialize M (B x P) and S (P x N) with random or uniform positive entries
NB. Input: (B ; N ; P)
init_nmf_factors =: 3 : 0
  'b n p' =. y
  NB. Uniform positive entries scaled to [0.1, 0.6]
  M_init =. 0.1 + 0.5 * ((b, p) $ (i. b * p) % (b * p))
  S_init =. (p, n) $ (1.0 % p)
  M_init ; S_init
)

NB. Compute Frobenius reconstruction loss: 0.5 * ||Y - M * S||_F^2
nmf_loss =: 3 : 0
  'Y M S' =. y
  resid =. Y - (M +/ .* S)
  0.5 * +/ , *: resid
)

NB. Multiplicative update step for NMF
NB. Y: (B x N), M: (B x P), S: (P x N)
nmf_update_step =: 3 : 0
  'Y M S' =. y
  eps =. 1e_9

  NB. 1. Update S: S = S * (M^T * Y) / (M^T * M * S + eps)
  Mt =. |: M
  MtY =. Mt +/ .* Y
  MtMS =. (Mt +/ .* M) +/ .* S
  S =. S * MtY % (MtMS + eps)

  NB. 2. Project abundances to satisfy ASC approximately: column normalization
  col_sum =. +/ S
  S =. S %"1 (col_sum + eps)

  NB. 3. Update M: M = M * (Y * S^T) / (M * S * S^T + eps)
  St =. |: S
  YSt =. Y +/ .* St
  MSSt =. (M +/ .* S) +/ .* St
  M =. M * YSt % (MSSt + eps)

  NB. Ensure strict non-negativity
  M =. 0.0 >. M
  S =. 0.0 >. S
  M ; S
)

NB. Run NMF blind spectral unmixing:
NB. x: (p_endmembers ; max_iterations)
NB. y: Y matrix (B bands x N pixels)
NB. Returns boxed: (M_endmembers ; S_abundances ; final_loss)
nmf_unmix =: 4 : 0
  'p max_iter' =. x
  Y =. y
  'b n' =. $ Y
  'M S' =. init_nmf_factors (b ; n ; p)

  for_iter. i. max_iter do.
    'M S' =. nmf_update_step (Y ; M ; S)
  end.

  loss =. nmf_loss (Y ; M ; S)
  M ; S ; loss
)
