NB. ==============================================================================
NB. Albedo - Automated Target & Anomaly Detection Pipeline
NB.
NB. Executes automated target search and anomaly triage:
NB.   1. Global RX anomaly detection across scene background.
NB.   2. Dual-window local RX detection for isolated subpixel anomalies.
NB.   3. Constrained Energy Minimization (CEM) matching against library targets.
NB.   4. Constant False Alarm Rate (CFAR) adaptive thresholding.
NB.   5. Extraction of pixel coordinates for verified target detections.
NB. ==============================================================================

cocurrent 'base'

load 'detect/rx.ijs'
load 'detect/local_rx.ijs'
load 'detect/matched_filter.ijs'
load 'detect/sam.ijs'

NB. Run target detection pipeline on surface reflectance cube
NB. Dyadic: target_signature run_detection_pipeline refl_cube
run_detection_pipeline=: 4 : 0
  target_sig=. x
  refl_cube=. y
  'B H W'=. $ refl_cube

  smoutput '[Albedo Detection] Step 1: Evaluating Global RX Anomalies...'
  rx_map=. albedo_rx_detect_global_rx refl_cube

  smoutput '[Albedo Detection] Step 2: Executing Constrained Energy Minimization...'
  cem_map=. target_sig albedo_target_detect_detect_cem refl_cube

  smoutput '[Albedo Detection] Step 3: Executing Adaptive Cosine Estimator...'
  ace_map=. target_sig albedo_target_detect_detect_ace refl_cube

  smoutput '[Albedo Detection] Step 4: CFAR Thresholding Top 0.5% Target Hits...'
  thresh_mask=. 0.995 albedo_rx_threshold_anomalies cem_map

  NB. Identify detected pixel coordinates (row, col)
  det_indices=. I. , thresh_mask
  det_coords=. 0 $ 0 0
  for_idx. det_indices do.
    r=. <. idx % W
    c=. W | idx
    det_coords=. det_coords , r, c
  end.
  smoutput '[Albedo Detection]   Identified ', (": # det_indices), ' Target Candidate Pixels.'

  NB. Returns: (cem_map; ace_map; rx_map; thresh_mask; det_coords)
  cem_map; ace_map; rx_map; thresh_mask; det_coords
)


albedo_pipeline_detect_run_detection_pipeline =: run_detection_pipeline
