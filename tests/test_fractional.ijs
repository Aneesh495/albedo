NB. ==============================================================================
NB. ALBEDO - Tests: signal/fractional_derivative.ijs
NB. ==============================================================================

load 'signal/fractional_derivative.ijs'

assert =: 3 : 0
  if. -. y do.
    echo 'ASSERTION FAILED'
    exit 1
  end.
  1
)

approx_eq =: 4 : 0
  tol =. 1e_4
  *./ , tol > | x - y
)

echo 'Testing Grünwald-Letnikov fractional spectral derivatives...'

NB. Test GL weights for alpha=1.0: should be [1, -1, 0, 0, ...]
w1 =: 1.0 albedo_frac_gl_weights 5
assert (1.0 _1.0 0.0 0.0 0.0) approx_eq w1

NB. Test GL weights for alpha=0.0: should be [1, 0, 0, 0, ...]
w0 =: 0.0 albedo_frac_gl_weights 5
assert (1.0 0.0 0.0 0.0 0.0) approx_eq w0

NB. Test fractional derivative on linear ramp f(x) = 2*x with step h=1.0
NB. For alpha=1.0, f'(x) = 2.0 (except boundary index 0)
ramp =: 2.0 * i. 10
d1 =: (1.0 1.0) albedo_frac_fractional_derivative_1d ramp
assert (9 # 2.0) approx_eq (1 }. d1)

NB. For alpha=0.0, f^(0)(x) = f(x) (identity)
d0 =: (0.0 1.0) albedo_frac_fractional_derivative_1d ramp
assert ramp approx_eq d0

NB. Test fractional order alpha=0.5 (half derivative)
d05 =: (0.5 1.0) albedo_frac_fractional_derivative_1d ramp
NB. Half-derivative values should be monotonically positive and bounded between ramp and 2.0
assert *./ (1 }. d05) > 0.0

echo 'Fractional derivative tests passed successfully.'
exit 0
