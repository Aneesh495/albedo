NB. ==============================================================================
NB. Albedo - Pixel Purity Index (PPI) Endmember Extraction
NB.
NB. Identifies pure spectral pixels by projecting data onto thousands of
NB. random skewers (unit vectors in multidimensional space) and counting
NB. how frequently each pixel falls at extreme boundaries.
NB. ==============================================================================

cocurrent 'base'

load 'core/linalg.ijs'

NB. ------------------------------------------------------------------------------
NB. Pixel Purity Index Algorithm
NB. ------------------------------------------------------------------------------

NB. Compute PPI scores and extract p endmembers
NB. Dyadic: (p; num_skewers; threshold_percentile) extract_ppi X
extract_ppi=: 4 : 0
  'p num_skewers thresh_pct'=. x
  X=. y
  'B N'=. $ X

  NB. Dimension reduction to p dimensions via SVD
  'svals U V'=. albedo_linalg_svd_power X
  Up=. (i. p) {"1 U
  Y=. (albedo_linalg_tr Up) albedo_linalg_matmul X

  NB. Initialize pixel hit counts
  hit_counts=. N $ 0

  for_s. i. num_skewers do.
    NB. Generate random unit skewer vector w in R^p
    raw_w=. 2.0 * ((? p $ 1000) % 1000.0) - 0.5
    w=. raw_w % ((albedo_linalg_norm2 raw_w) >. 1e_8)

    NB. Project all pixels: proj = w^T Y (length N)
    proj=. w albedo_linalg_matmul Y

    NB. Find minimum and maximum extreme indices
    min_idx=. (<: +/ proj = <./ proj) { (I. proj = <./ proj)
    max_idx=. (<: +/ proj = >./ proj) { (I. proj = >./ proj)

    hit_counts=. (1 + min_idx { hit_counts) min_idx} hit_counts
    hit_counts=. (1 + max_idx { hit_counts) max_idx} hit_counts
  end.

  NB. Sort pixels by descending hit count
  sorted_indices=. \: hit_counts
  top_endmember_indices=. p {. sorted_indices

  M_endmembers=. top_endmember_indices {"1 X
  M_endmembers; top_endmember_indices; hit_counts
)


albedo_ppi_extract_ppi =: extract_ppi
