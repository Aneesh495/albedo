NB. ==============================================================================
NB. Test Suite: detect/rx.ijs, detect/local_rx.ijs, detect/matched_filter.ijs, detect/sam.ijs
NB. ==============================================================================

load 'tests/test_framework.ijs'
load 'detect/rx.ijs'
load 'detect/local_rx.ijs'
load 'detect/matched_filter.ijs'
load 'detect/sam.ijs'

cocurrent 'base'

suite_begin 'Target & Anomaly Detection Algorithms'

NB. Construct background scene (8 bands, 16 x 16) with normal random distribution
Cube_bg=: ? 8 16 16 $ 100
NB. Inject rare anomalous target signature at pixel (8, 8)
Target_sig=: 8 $ 500.0
Cube_injected=: Cube_bg
Cube_injected=: Target_sig (<a:; 8; 8)} Cube_injected

NB. Test 1: Global RX Anomaly Detector
Rx_map=: albedo_rx_detect_global_rx Cube_injected
'RX map shape is 16 16' assert_shape Rx_map; 16 16
'Anomaly pixel (8,8) has highest RX score' assert_true ((((< 8 8) { Rx_map)) = >./ , Rx_map)

NB. Threshold top 1% anomalies
Mask_anom=: 0.99 albedo_rx_threshold_anomalies Rx_map
'Injected pixel (8,8) is flagged as anomaly' assert_true (((< 8 8) { Mask_anom))

NB. Test 2: Dual-Window Local RX Detector
Local_rx=: (3; 7) albedo_local_rx_detect_local_rx Cube_injected
'Local RX map shape is 16 16' assert_shape Local_rx; 16 16
'Anomaly pixel (8,8) detected locally' assert_true ((((< 8 8) { Local_rx)) > 10.0)

NB. Test 3: Constrained Energy Minimization (CEM)
Cem_scores=: Target_sig albedo_target_detect_detect_cem Cube_injected
'CEM map shape is 16 16' assert_shape Cem_scores; 16 16
'CEM response on target pixel near 1.0' assert_near (((< 8 8) { Cem_scores)); 1.0; 0.1
'CEM response on target is highest in scene' assert_true ((((< 8 8) { Cem_scores)) = >./ , Cem_scores)

NB. Test 4: Spectral Matched Filter (MF)
Mf_scores=: Target_sig albedo_target_detect_detect_matched_filter Cube_injected
'MF response on target is highest in scene' assert_true ((((< 8 8) { Mf_scores)) = >./ , Mf_scores)

NB. Test 5: Adaptive Cosine Estimator (ACE)
Ace_scores=: Target_sig albedo_target_detect_detect_ace Cube_injected
'ACE scores strictly in [0, 1]' assert_true (*./ (*./ , Ace_scores >: 0.0) *. (*./ , Ace_scores <: 1.0))
'ACE score on target pixel near 1.0' assert_near (((< 8 8) { Ace_scores)); 1.0; 0.1

NB. Test 6: Spectral Angle Mapper (SAM)
VecA=: 1.0 2.0 3.0 4.0
VecB=: 2.0 4.0 6.0 8.0
AngleCollinear=: VecA albedo_similarity_spectral_angle VecB
'Collinear vectors have 0 angle' assert_near AngleCollinear; 0.0; 1e_6

VecOrth=: 0.0 0.0 0.0 5.0
AngleOrth=: VecA albedo_similarity_spectral_angle VecOrth
'Orthogonal angle is positive' assert_true (AngleOrth > 0.0)

Sam_map=: Target_sig albedo_similarity_sam_cube Cube_injected
'SAM map shape is 16 16' assert_shape Sam_map; 16 16
'Target pixel has minimal angle' assert_true ((((< 8 8) { Sam_map)) = <./ , Sam_map)

NB. Test 7: Spectral Information Divergence (SID)
Sid_self=: VecA albedo_similarity_spectral_information_divergence VecA
'SID of vector with itself is 0' assert_near Sid_self; 0.0; 1e_6

Sid_diff=: VecA albedo_similarity_spectral_information_divergence (4 $ 1.0)
'SID between different distributions is positive' assert_true (Sid_diff > 0.0)

suite_end ''
