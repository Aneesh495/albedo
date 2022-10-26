NB. ==============================================================================
NB. ALBEDO - High-Dimensional Hyperspectral Imaging Engine
NB. signal/wavelet.ijs : 1D & 2D Discrete Wavelet Transform (DWT) & Wavelet Shrinkage
NB. Implements Haar and Daubechies-4 (db4) Mallat filter banks, multi-level pyramid
NB. decomposition, universal VisuShrink thresholding, and 2D spatial subband analysis.
NB. ==============================================================================

cocurrent 'albedo_wavelet'

NB. Filter coefficients
NB. Haar filter bank
SQRT2 =: %: 2
HAAR_H =: (1.0 , 1.0) % SQRT2
HAAR_G =: (1.0 , _1.0) % SQRT2

NB. Daubechies-4 (db4) scaling and wavelet coefficients
SQRT3 =: %: 3
DENOM =: 4 * SQRT2
DB4_H =: ((1 + SQRT3) , (3 + SQRT3) , (3 - SQRT3) , (1 - SQRT3)) % DENOM
DB4_G =: (DENOM %~ 1 - SQRT3) , (DENOM %~ - 3 - SQRT3) , (DENOM %~ 3 + SQRT3) , (DENOM %~ - 1 + SQRT3)

NB. 1D Haar forward transform step on even-length vector
NB. Yields boxed (approx ; detail)
dwt1_haar_step =: 3 : 0
  n =. # y
  evens =. (2 * i. -: n) { y
  odds  =. (1 + 2 * i. -: n) { y
  approx =. (evens + odds) % SQRT2
  detail =. (evens - odds) % SQRT2
  approx ; detail
)

NB. 1D Haar inverse transform step
NB. Input boxed (approx ; detail)
idwt1_haar_step =: 3 : 0
  'approx detail' =. y
  m =. # approx
  evens =. (approx + detail) % SQRT2
  odds  =. (approx - detail) % SQRT2
  NB. Interleave evens and odds
  , (evens ,: odds)
)

NB. Circular periodic convolution with step downsampling
NB. x: filter (reversed for reconstruction or forward), y: signal
conv_down2 =: 4 : 0
  h =. x
  sig =. y
  n =. # sig
  k =. # h
  half_n =. -: n
  out =. half_n # 0.0
  for_i. i. half_n do.
    idx =. 2 * i
    sum =. 0.0
    for_j. i. k do.
      shift_idx =. n | idx + j
      sum =. sum + (j { h) * shift_idx { sig
    end.
    out =. sum i} out
  end.
  out
)

NB. Up-sampling by 2 and circular convolution
upsample_conv =: 4 : 0
  h =. x
  sig =. y
  m =. # sig
  n =. 2 * m
  k =. # h
  NB. Upsample with zeros
  up =. , (sig ,: (m # 0.0))
  out =. n # 0.0
  for_i. i. n do.
    sum =. 0.0
    for_j. i. k do.
      shift_idx =. n | i - j
      sum =. sum + (j { h) * shift_idx { up
    end.
    out =. sum i} out
  end.
  out
)

NB. 1D db4 forward step
dwt1_db4_step =: 3 : 0
  approx =. DB4_H conv_down2 y
  detail =. DB4_G conv_down2 y
  approx ; detail
)

NB. Multi-level 1D DWT decomposition
NB. Input: signal y, levels: x
NB. Returns boxed list: (approx_L ; detail_L ; detail_(L-1) ; ... ; detail_1)
dwt1_multilevel =: 4 : 0
  levels =. x
  cur =. y
  details =. 0 $ < 0
  for_l. i. levels do.
    'a d' =. dwt1_haar_step cur
    details =. details , < d
    cur =. a
  end.
  (< cur) , details
)

NB. Multi-level 1D IDWT reconstruction
idwt1_multilevel =: 3 : 0
  levels =. <: # y
  cur =. > 0 { y
  for_l. i. levels do.
    d =. > (1 + l) { y
    cur =. idwt1_haar_step (cur ; d)
  end.
  cur
)

NB. Soft thresholding operator: sign(x) * max(0, |x| - tau)
soft_threshold =: 4 : 0
  tau =. x
  (* y) * 0.0 >. (| y) - tau
)

NB. Hard thresholding operator: x * (|x| >= tau)
hard_threshold =: 4 : 0
  tau =. x
  y * (| y) >: tau
)

NB. Universal VisuShrink threshold: sigma * sqrt(2 * ln(N))
NB. where sigma is estimated via MAD of finest detail coefficients: median(|d|) / 0.6745
visushrink_threshold =: 3 : 0
  d1 =. y
  n =. # d1
  NB. Median Absolute Deviation (MAD)
  sorted_abs =. /:~ | d1
  mid =. <. -: n
  mad =. mid { sorted_abs
  sigma =. mad % 0.6745
  sigma * %: 2 * ^. n
)

NB. 1D Wavelet Denoising via Haar / soft thresholding
wavelet_denoise_1d =: 3 : 0
  n =. # y
  'a d' =. dwt1_haar_step y
  tau =. visushrink_threshold d
  d_clean =. tau soft_threshold d
  idwt1_haar_step (a ; d_clean)
)

NB. 2D Single-Level DWT (Haar)
NB. Decomposes matrix M (rows x cols) into 4 subbands: LL (approx), LH (horiz), HL (vert), HH (diag)
dwt2_haar =: 3 : 0
  'nr nc' =. $ y
  NB. Apply 1D DWT on each row
  row_decomp =. dwt1_haar_step"1 y
  NB. Unbox rows: row_a has shape (nr, nc/2), row_d has shape (nr, nc/2)
  row_a =. > 0 {"1 row_decomp
  row_d =. > 1 {"1 row_decomp
  NB. Apply 1D DWT on columns (transpose, transform rows, transpose back)
  col_aa =. |: > 0 {"1 dwt1_haar_step"1 (|: row_a)
  col_ad =. |: > 1 {"1 dwt1_haar_step"1 (|: row_a)
  col_da =. |: > 0 {"1 dwt1_haar_step"1 (|: row_d)
  col_dd =. |: > 1 {"1 dwt1_haar_step"1 (|: row_d)
  NB. Return LL, LH, HL, HH
  col_aa ; col_ad ; col_da ; col_dd
)

NB. 2D Single-Level IDWT (Haar)
NB. Reconstructs matrix from (LL ; LH ; HL ; HH)
idwt2_haar =: 3 : 0
  'll lh hl hh' =. y
  NB. Reconstruct columns
  row_a =. |: idwt1_haar_step"1 (|: ll) ,: (|: lh)
  row_d =. |: idwt1_haar_step"1 (|: hl) ,: (|: hh)
  NB. Reconstruct rows
  idwt1_haar_step"1 row_a ,: row_d
)
