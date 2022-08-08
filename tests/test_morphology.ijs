NB. ==============================================================================
NB. Test Suite: core/morphology.ijs
NB. ==============================================================================

load 'tests/test_framework.ijs'
load 'core/morphology.ijs'

cocurrent 'base'

suite_begin 'Grayscale Mathematical Morphology & Morphological Profiles'

NB. Create 7x7 image with background 0 and single bright impulse at (3, 3)
ImageImpulse=: 7 7 $ 0.0
ImageImpulse=: 10.0 (<3 3)} ImageImpulse

NB. Test 1: Dilation expands impulse
Dilated1=: 1 albedo_morphology_dilate_2d ImageImpulse
'Dilated shape is 7 7' assert_shape Dilated1; 7 7
'Dilated neighborhood has 10.0' assert_equal (Dilated1 < 2 3); 10.0
'Dilated neighborhood has 10.0' assert_equal (Dilated1 < 3 4); 10.0
'Dilated far corner is still 0.0' assert_equal (Dilated1 < 0 0); 0.0

NB. Test 2: Erosion removes isolated impulse
Eroded1=: 1 albedo_morphology_erode_2d ImageImpulse
'Eroded center is 0.0' assert_equal (Eroded1 < 3 3); 0.0

NB. Test 3: Opening removes isolated bright impulse
Opened1=: 1 albedo_morphology_open_2d ImageImpulse
'Opening removes isolated impulse' assert_near Opened1; (7 7 $ 0.0); 1e_6

NB. Test 4: White Top-Hat isolates impulse
TopHatW=: 1 albedo_morphology_top_hat_white ImageImpulse
'White Top-Hat extracts original impulse' assert_near TopHatW; ImageImpulse; 1e_6

NB. Test 5: Black Top-Hat on dark hole
ImageHole=: 7 7 $ 10.0
ImageHole=: 0.0 (<3 3)} ImageHole
TopHatB=: 1 albedo_morphology_top_hat_black ImageHole
'Black Top-Hat extracts hole depth 10.0 at (3,3)' assert_near (TopHatB < 3 3); 10.0; 1e_6
'Black Top-Hat is 0.0 away from hole' assert_near (TopHatB < 0 0); 0.0; 1e_6

NB. Test 6: Multiscale Morphological Profile
MP=: (1 2) albedo_morphology_build_mp_2d ImageImpulse
'MP shape is 5 7 7 (2*2 + 1 bands)' assert_shape MP; 5 7 7
'MP band 2 is original image' assert_equal (2 { MP); ImageImpulse

suite_end ''
