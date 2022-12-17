NB. ==============================================================================
NB. Test Suite: unmix/endmember_bundles.ijs
NB. ==============================================================================

load 'tests/test_framework.ijs'
load 'unmix/endmember_bundles.ijs'
load 'sim/synthetic_scene.ijs'
load 'core/linalg.ijs'

cocurrent 'base'

suite_begin 'Endmember Bundles & Extended Linear Mixing Model (ELMM)'

'RadCube ReflNoisy A_gt M_true Wvls Plume'=: albedo_sim_generate_scene (16; 16; 30.0)

NB. Test 1: Bundle Extraction across 2x2 tiles
Bundles=: (2; 2; 2) albedo_bundles_extract_bundles ReflNoisy
'Bundles band dimension is 224' assert_equal ({. $ Bundles); 224
'Extracted at least 4 bundle candidate signatures' assert_true (({: $ Bundles) >: 4)

NB. Test 2: ELMM Scaling Invariance
NB. Test pixel with 50% shading factor: y_shaded = 0.5 * (M s)
s_true=: 0.6 0.2 0.2 0.0 0.0
y_clean=: M_true albedo_linalg_matmul s_true
y_shaded=: 0.5 * y_clean

's_est psi_est'=: M_true albedo_bundles_solve_elmm y_shaded
'ELMM recovers shading factor psi near 0.5' assert_near psi_est; 0.5; 0.05
'ELMM recovers abundance fractions' assert_near s_est; s_true; 0.05

NB. Test 3: Batch ELMM Unmixing on 3D Cube
'Abundances ScaleMap'=: M_true albedo_bundles_unmix_elmm_cube ((<a:; (i. 4); (i. 4)) { ReflNoisy)
'Abundances shape is 5 4 4' assert_shape Abundances; 5 4 4
'ScaleMap shape is 4 4' assert_shape ScaleMap; 4 4
'Scaling factors are strictly positive' assert_true (*./ (*./ ScaleMap > 0.1))
'All pixel abundances sum to 1.0' assert_near (+/ Abundances); (4 4 $ 1.0); 1e_4

suite_end ''
