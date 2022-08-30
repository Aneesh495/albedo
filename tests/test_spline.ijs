NB. ==============================================================================
NB. Test Suite: core/spline.ijs
NB. ==============================================================================

load 'tests/test_framework.ijs'
load 'core/spline.ijs'

cocurrent 'base'

suite_begin 'Natural Cubic Spline Interpolation & Resampling'

X_nodes=: 1.0 2.0 3.0 4.0 5.0
Y_nodes=: 2.0 5.0 1.0 8.0 4.0

NB. Test 1: Spline Model Fit
SplineModel=: albedo_spline_fit_natural_spline (X_nodes; Y_nodes)
'Spline model has 4 elements' assert_equal (# SplineModel); 4
'M second derivatives'=: 2 {:: SplineModel
'Natural spline M0 = 0' assert_near ({. M second derivatives); 0.0; 1e_6
'Natural spline M_last = 0' assert_near ({: M second derivatives); 0.0; 1e_6

NB. Test 2: Interpolation at Nodes (Exact Match S(x_i) = y_i)
Y_at_nodes=: X_nodes albedo_spline_eval_spline SplineModel
'Spline passes exactly through input nodes' assert_near Y_at_nodes; Y_nodes; 1e_5

NB. Test 3: Interpolation at Intermediate Points
X_mid=: 1.5 2.5 3.5 4.5
Y_mid=: X_mid albedo_spline_eval_spline SplineModel
'Interpolated 4 intermediate points' assert_equal (# Y_mid); 4
'Intermediate points are finite numbers' assert_true (*./ | Y_mid < 100.0)

NB. Test 4: Linear Function Preserved Exactly
X_lin=: 0.0 1.0 2.0 3.0 4.0
Y_lin=: 3.0 * X_lin
LinModel=: albedo_spline_fit_natural_spline (X_lin; Y_lin)
X_query=: 0.5 1.5 2.5 3.5
Y_query=: X_query albedo_spline_eval_spline LinModel
'Spline preserves linear function exactly' assert_near Y_query; (3.0 * X_query); 1e_5

suite_end ''
