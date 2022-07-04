NB. ==============================================================================
NB. Albedo - Spatial-Spectral Regularization & Majority Voting
NB.
NB. Applies morphological filtering and 2D spatial majority voting to remove
NB. isolated salt-and-pepper classification noise and enforce spatial coherence.
NB. ==============================================================================

cocurrent 'albedo_spatial_spectral'

NB. Mode (most frequent element) of an integer vector
vector_mode=: 3 : 0
  v=. y
  unique_vals=. ~. v
  counts=. +/ "1 unique_vals =/ v
  max_count=. >./ counts
  max_idx=. (<: +/ counts = max_count) { (I. counts = max_count)
  max_idx { unique_vals
)

NB. Spatial majority voting filter with window size (e.g. 3x3 or 5x5)
NB. Dyadic: window_size majority_filter label_map
majority_filter=: 4 : 0
  w_size=. x
  labels=. y
  'H W'=. $ labels
  assert. 1 = 2 | w_size [ 'Majority filter window size must be odd'
  half_w=. <. w_size % 2

  cleaned=. labels

  for_h. half_w + i. (H - 2 * half_w) do.
    for_w. half_w + i. (W - 2 * half_w) do.
      rows=. (h - half_w) + i. w_size
      cols=. (w - half_w) + i. w_size
      neighborhood=. , (labels < rows; cols)
      NB. Filter out unknown/unclassified pixels (-1)
      valid_neighbors=. (#~ _1&~:) neighborhood
      if. 0 < # valid_neighbors do.
        m_val=. vector_mode valid_neighbors
        cleaned=. m_val (<h, w)} cleaned
      end.
    end.
  end.
  cleaned
)
