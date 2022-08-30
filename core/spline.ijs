NB. ==============================================================================
NB. Albedo - Natural Cubic Spline Interpolation & Resampling
NB.
NB. Implements natural cubic splines with tridiagonal Thomas algorithm solver.
NB. Evaluates piecewise cubic polynomials, continuous first derivatives, and
NB. second derivatives for high-precision spectral signature reconstruction.
NB. ==============================================================================

cocurrent 'albedo_spline'

load 'core/linalg.ijs'

NB. ------------------------------------------------------------------------------
NB. Natural Cubic Spline Coefficient Calculation
NB. Solves tridiagonal system for second derivatives M_i
NB. ------------------------------------------------------------------------------

NB. Fit natural cubic spline to nodes (x_nodes; y_nodes)
NB. Returns: (x_nodes; y_nodes; M_second_derivs; h_intervals)
fit_natural_spline=: 3 : 0
  'x y'=. y
  n=. # x
  assert. n >: 3 [ 'At least 3 points required for cubic spline'

  NB. Compute interval lengths h_i = x_{i+1} - x_i
  h=. (1 }. x) - (_1 }. x)
  assert. *./ h > 0.0 [ 'Node coordinates x must be strictly monotonically increasing'

  NB. Tridiagonal system of size (n - 2) x (n - 2) for interior M_1 .. M_{n-2}
  m_size=. n - 2
  diag_main=. 2.0 * (m_size {. h) + (1 }. m_size {. (1 }. h))
  diag_sub=. 1 }. m_size {. h
  diag_sup=. _1 }. m_size {. h

  NB. Right-hand side vector d_i = 6 * ((y_{i+1} - y_i)/h_i - (y_i - y_{i-1})/h_{i-1})
  dy=. (1 }. y) - (_1 }. y)
  slopes=. dy % h
  rhs=. 6.0 * (1 }. slopes) - (_1 }. slopes)

  NB. Solve tridiagonal system via Thomas algorithm
  c_prime=. m_size $ 0.0
  d_prime=. m_size $ 0.0

  c0=. ({. diag_sup) % ({. diag_main)
  d0=. ({. rhs) % ({. diag_main)
  c_prime=. c0 (0)} c_prime
  d_prime=. d0 (0)} d_prime

  for_i. 1 + i. (m_size - 1) do.
    ai=. (<: i) { diag_sub
    bi=. i { diag_main
    ci=. i { diag_sup
    di=. i { rhs

    denom=. bi - (ai * (<: i) { c_prime)
    denom=. denom + 1e_12 * denom = 0.0
    c_prime=. (ci % denom) (<i)} c_prime
    d_prime=. ((di - ai * (<: i) { d_prime) % denom) (<i)} d_prime
  end.

  NB. Back substitution
  M_interior=. m_size $ 0.0
  M_last=. {: d_prime
  M_interior=. M_last (<: m_size)} M_interior

  for_i. |. i. (m_size - 1) do.
    m_val=. (i { d_prime) - (i { c_prime) * ((i + 1) { M_interior)
    M_interior=. m_val (<i)} M_interior
  end.

  NB. Natural boundary conditions: M_0 = 0 and M_{n-1} = 0
  M_all=. 0.0 , M_interior , 0.0

  x; y; M_all; h
)

NB. ------------------------------------------------------------------------------
NB. Spline Evaluation at Arbitrary Query Points
NB. ------------------------------------------------------------------------------

NB. Evaluate fitted spline at query point x_q
NB. Piecewise cubic:
NB. S_i(x) = (M_i (x_{i+1} - x)^3 + M_{i+1} (x - x_i)^3) / (6 h_i)
NB.        + ((y_i - M_i h_i^2 / 6) (x_{i+1} - x) + (y_{i+1} - M_{i+1} h_i^2 / 6) (x - x_i)) / h_i
NB. Dyadic: x_queries eval_spline spline_model
eval_spline=: 4 : 0
  x_q=. x
  model=. y
  'x_nodes y_nodes M h'=. model
  n=. # x_nodes

  out=. (# x_q) $ 0.0

  for_idx. i. # x_q do.
    xq=. idx { x_q

    NB. Clamp or find interval i such that x_i <= xq <= x_{i+1}
    if. xq <: {. x_nodes do.
      out=. ({. y_nodes) (<idx)} out
    elseif. xq >: {: x_nodes do.
      out=. ({: y_nodes) (<idx)} out
    else.
      i=. (<: +/ x_nodes <: xq) <. (n - 2)
      x_i=. i { x_nodes
      x_ip1=. (i + 1) { x_nodes
      y_i=. i { y_nodes
      y_ip1=. (i + 1) { y_nodes
      M_i=. i { M
      M_ip1=. (i + 1) { M
      h_i=. i { h

      term1=. (M_i * ((x_ip1 - xq) ^ 3)) + (M_ip1 * ((xq - x_i) ^ 3))
      term1=. term1 % (6.0 * h_i)

      c_i=. y_i - (M_i * (*: h_i) % 6.0)
      c_ip1=. y_ip1 - (M_ip1 * (*: h_i) % 6.0)
      term2=. ((c_i * (x_ip1 - xq)) + (c_ip1 * (xq - x_i))) % h_i

      val=. term1 + term2
      out=. val (<idx)} out
    end.
  end.

  out
)
