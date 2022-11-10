NB. ==============================================================================
NB. ALBEDO - High-Dimensional Hyperspectral Imaging Engine
NB. signal/fractional_derivative.ijs : Grünwald-Letnikov Fractional Spectral Differentiation
NB. Implements arbitrary real-order alpha in [0.0, 2.0] fractional calculus operators
NB. for spectral curve enhancement, slope extraction, and noise-attenuated inflection analysis.
NB. ==============================================================================

cocurrent 'albedo_frac'

NB. Generate Grünwald-Letnikov weights for order alpha and memory length m
NB. c_0 = 1
NB. c_j = c_(j-1) * (j - 1 - alpha) / j
gl_weights =: 4 : 0
  alpha =. x
  m =. y
  w =. m # 0.0
  cur =. 1.0
  w =. cur 0} w
  for_j. 1 + i. <: m do.
    cur =. cur * (j - 1.0 + alpha) % j NB. sign alternation handled
    NB. For Grünwald-Letnikov: c_j = (-1)^j * binom(alpha, j)
    NB. binom(alpha, j) = binom(alpha, j-1) * (alpha - j + 1) / j
    NB. (-1)^j * binom = (-1) * (-1)^(j-1) * binom = c_(j-1) * (j - 1 - alpha) / j
  end.
  NB. Recalculate cleanly:
  c =. m # 0.0
  val =. 1.0
  c =. val 0} c
  for_j. 1 + i. <: m do.
    val =. val * (j - 1.0 - alpha) % j
    c =. val j} c
  end.
  c
)

NB. Compute fractional derivative of 1D spectrum:
NB. x: (alpha ; step_size_h)
NB. y: spectral vector (length N)
fractional_derivative_1d =: 4 : 0
  'alpha h' =. x
  spec =. y
  n =. # spec
  h_factor =. h ^ alpha
  NB. Truncate filter memory to min(n, 32) for efficiency and stability
  m =. n <. 32
  w =. alpha gl_weights m
  out =. n # 0.0
  for_i. i. n do.
    k =. (1 + i) <. m
    sub_w =. k {. w
    sub_s =. |. (i - i. k) { spec
    out =. ((+/ sub_w * sub_s) % h_factor) i} out
  end.
  out
)

NB. Apply fractional derivative across all spectra in a 2D matrix (samples x bands)
fractional_derivative_matrix =: 4 : 0
  'alpha h' =. x
  mat =. y
  (alpha , h) & fractional_derivative_1d "1 mat
)
