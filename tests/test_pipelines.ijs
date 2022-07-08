NB. ==============================================================================
NB. Test Suite: sim/synthetic_scene.ijs, pipeline/pipeline_unmix.ijs, pipeline/pipeline_detect.ijs
NB. ==============================================================================

load 'tests/test_framework.ijs'
load 'sim/synthetic_scene.ijs'
load 'pipeline/pipeline_unmix.ijs'
load 'pipeline/pipeline_detect.ijs'

cocurrent 'base'

suite_begin 'Synthetic Scene Generation & End-to-End Processing Pipelines'

NB. Test 1: Synthetic Scene Generation
'RadCube ReflNoisy A_gt M_true Wvls PlumeCoord'=: albedo_sim_generate_scene (16; 16; 30.0)
'RadCube shape is 224 16 16' assert_shape RadCube; 224 16 16
'ReflNoisy shape is 224 16 16' assert_shape ReflNoisy; 224 16 16
'Ground-truth abundance shape is 5 16 16' assert_shape A_gt; 5 16 16
'Ground-truth abundances sum to 1.0' assert_near (+/ A_gt); (16 16 $ 1.0); 1e_4

NB. Test 2: End-to-End Unmixing Pipeline
'BoaCube M_est Abundances dim_k rmse'=: (Wvls; 120; 35.0) albedo_pipeline_unmix_run_unmixing_pipeline RadCube
'BoaCube shape is 224 16 16' assert_shape BoaCube; 224 16 16
'Identified dimension k > 0' assert_true (dim_k > 0)
'Abundances shape has 16 16 spatial dimension' assert_equal (1 2 { $ Abundances); 16 16
'All estimated abundances satisfy non-negativity (ANC)' assert_true (*./ (*./ , Abundances >: _1e_12))
'All estimated abundances satisfy sum-to-one (ASC)' assert_near (+/ Abundances); (16 16 $ 1.0); 1e_4
'Unmixing RMSE is low (< 0.15)' assert_true (rmse < 0.15)

NB. Test 3: End-to-End Target Detection Pipeline (Methane Search)
Methane_sig=: albedo_spectral_lib_ref_methane_plume Wvls
'CemMap AceMap RxMap DetMask DetCoords'=: Methane_sig albedo_pipeline_detect_run_detection_pipeline BoaCube
'CEM map shape is 16 16' assert_shape CemMap; 16 16
'ACE map shape is 16 16' assert_shape AceMap; 16 16
'At least one target hit detected' assert_true (0 < # DetCoords)
NB. Verify detected coordinate contains the injected plume (around row 8, col 4)
'Plume area detected by pipeline' assert_true (+/ (DetMask < 8 4) , (DetMask < 8 5) , (DetMask < 9 4))

suite_end ''
