NB. ==============================================================================
NB. Test Suite: core/quadrature.ijs
NB. ==============================================================================

load 'tests/test_framework.ijs'
load 'core/quadrature.ijs'

cocurrent 'base'

suite_begin 'Numerical Quadrature & Integration Routines'

NB. Test 1: Trapezoidal rule on linear function y = 2x from 0 to 3: int = [x^2]_0^3 = 9.0
X_grid=: 0.0 1.0 2.0 3.0
Y_vals=: 2.0 * X_grid
IntTrapz=: X_grid albedo_quadrature_trapz Y_vals
'Trapezoidal integral of 2x on [0, 3] is 9.0' assert_near IntTrapz; 9.0; 1e_6

NB. Test 2: Non-uniform grid trapezoidal rule
X_nonunif=: 0.0 0.5 2.0 3.0
Y_nonunif=: 2.0 * X_nonunif
IntNonUnif=: X_nonunif albedo_quadrature_trapz Y_nonunif
'Non-uniform trapezoidal integral of 2x on [0, 3] is 9.0' assert_near IntNonUnif; 9.0; 1e_6

NB. Test 3: Simpson 1/3 rule on cubic polynomial y = x^3 from 0 to 2: int = [x^4 / 4]_0^2 = 4.0
X_simp=: 0.0 0.5 1.0 1.5 2.0
Y_simp=: X_simp ^ 3
IntSimp=: 0.5 albedo_quadrature_simpson_13 Y_simp
'Simpson 1/3 integral of x^3 on [0, 2] is exactly 4.0' assert_near IntSimp; 4.0; 1e_6

NB. Test 4: Gauss-Legendre 3-point rule on polynomial y = 3 x^2 from 0 to 2: int = [x^3]_0^2 = 8.0
poly_quad=: 3 : '3.0 * *: y'
IntGauss3=: (0.0 2.0) albedo_quadrature_gauss3_integrate poly_quad
'Gauss 3-point integral of 3x^2 on [0, 2] is 8.0' assert_near IntGauss3; 8.0; 1e_6

NB. Test 5: Gauss-Legendre 5-point rule on polynomial degree 4: y = 5 x^4 on [0, 1]: int = 1.0
poly_deg4=: 3 : '5.0 * y ^ 4'
IntGauss5=: (0.0 1.0) albedo_quadrature_gauss5_integrate poly_deg4
'Gauss 5-point integral of 5x^4 on [0, 1] is 1.0' assert_near IntGauss5; 1.0; 1e_6

suite_end ''
