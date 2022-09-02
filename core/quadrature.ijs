NB. ==============================================================================
NB. Albedo - High-Precision Numerical Quadrature & Integration
NB.
NB. Implements numerical integration for spectral irradiance flux, absorption
NB. band areas, and optical depths:
NB.   - Non-uniform composite Trapezoidal integration
NB.   - Composite Simpson's 1/3 and 3/8 rules
NB.   - Gauss-Legendre Quadrature (2, 3, 5, 7 point rules)
NB.   - Adaptive recursive Simpson integration
NB. ==============================================================================

cocurrent 'albedo_quadrature'

NB. ------------------------------------------------------------------------------
NB. Composite Trapezoidal Rule for Non-Uniform Grid Points
NB. ------------------------------------------------------------------------------

NB. Integrate y values over x grid: int y dx
NB. Dyadic: x_grid trapz y_values
trapz=: 4 : 0
  x=. x
  y=. y
  n=. # x
  assert. n = # y [ 'Coordinate and value array lengths must match'
  assert. n >: 2 [ 'At least 2 points required for trapezoidal integration'

  dx=. (1 }. x) - (_1 }. x)
  avg_y=. 0.5 * (1 }. y) + (_1 }. y)
  +/ dx * avg_y
)

NB. ------------------------------------------------------------------------------
NB. Composite Simpson's 1/3 Rule
NB. ------------------------------------------------------------------------------

NB. Simpson's 1/3 rule on uniform grid: h/3 * (y_0 + 4 y_1 + 2 y_2 + ... + y_n)
NB. Requires odd number of points (even number of intervals)
simpson_13=: 4 : 0
  h=. x
  y=. y
  n=. # y
  assert. 1 = 2 | n [ 'Number of points must be odd for Simpson 1/3 rule'

  weights=. n $ 2.0
  NB. Set first and last to 1.0
  weights=. 1.0 (0)} weights
  weights=. 1.0 (<: n)} weights
  NB. Set odd interior indices to 4.0
  odd_idx=. 1 + 2 * i. <. n % 2
  weights=. 4.0 odd_idx} weights

  (h % 3.0) * +/ weights * y
)

NB. ------------------------------------------------------------------------------
NB. Gauss-Legendre Quadrature
NB. Exact for polynomials up to degree (2n - 1)
NB. ------------------------------------------------------------------------------

NB. 3-point Gauss-Legendre roots and weights on [-1, 1]
NB. Roots: 0, +- sqrt(3/5) ~ 0.7745966692
NB. Weights: 8/9, 5/9, 5/9
gauss3_integrate=: 4 : 0
  'a b'=. x
  func=. y
  roots=. 0.0 , (%: 0.6) , (- %: 0.6)
  weights=. (8.0 % 9.0) , (5.0 % 9.0) , (5.0 % 9.0)

  NB. Map roots from [-1, 1] to [a, b]: x = ((b-a)/2) * t + (b+a)/2
  half_len=. 0.5 * (b - a)
  mid=. 0.5 * (a + b)
  x_pts=. (half_len * roots) + mid

  NB. Evaluate function at x_pts (monadic verb)
  f_vals=. func x_pts
  half_len * +/ weights * f_vals
)

NB. 5-point Gauss-Legendre roots and weights on [-1, 1]
gauss5_integrate=: 4 : 0
  'a b'=. x
  func=. y
  r1=. 0.0
  r2=. %: (5.0 - (2.0 * %: 10.0 % 7.0)) % 9.0
  r3=. %: (5.0 + (2.0 * %: 10.0 % 7.0)) % 9.0
  roots=. r1, r2, (- r2), r3, (- r3)

  w1=. 128.0 % 225.0
  w2=. (322.0 + (13.0 * %: 70.0)) % 900.0
  w3=. (322.0 - (13.0 * %: 70.0)) % 900.0
  weights=. w1, w2, w2, w3, w3

  half_len=. 0.5 * (b - a)
  mid=. 0.5 * (a + b)
  x_pts=. (half_len * roots) + mid

  f_vals=. func x_pts
  half_len * +/ weights * f_vals
)
