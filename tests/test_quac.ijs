NB. ==============================================================================
NB. Test Suite: atmos/quac.ijs
NB. ==============================================================================

load 'tests/test_framework.ijs'
load 'atmos/quac.ijs'
load 'sim/synthetic_scene.ijs'

cocurrent 'base'

suite_begin 'Quick Atmospheric Correction (QUAC)'

'RadCube ReflTrue A_gt M_true Wvls Plume'=: albedo_sim_generate_scene (16; 16; 30.0)

NB. Test 1: Run QUAC with 4 endmembers
'QuacRefl Gains Offsets'=: 4 albedo_quac_quac_correction RadCube
'QUAC reflectance shape is 224 16 16' assert_shape QuacRefl; 224 16 16
'QUAC gains count is 224' assert_equal (# Gains); 224
'QUAC offsets count is 224' assert_equal (# Offsets); 224
'QUAC gains are positive' assert_true (*./ Gains > 0.0)
'QUAC offsets are non-negative' assert_true (*./ Offsets >: 0.0)
'QUAC reflectance strictly bounded in [0, 1]' assert_true (*./ (*./ , QuacRefl >: 0.0) *. (*./ , QuacRefl <: 1.0))

suite_end ''
