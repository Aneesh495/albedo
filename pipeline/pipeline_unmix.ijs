NB. ==============================================================================
NB. Albedo - Automated End-to-End Hyperspectral Unmixing Pipeline
NB.
NB. Executes automated scientific processing pipeline:
NB.   1. Radiometric calibration and DOS-1 atmospheric haze correction.
NB.   2. Intrinsic subspace dimension estimation via HySime.
NB.   3. Pure endmember extraction via Vertex Component Analysis (VCA).
NB.   4. Fully Constrained Least Squares (FCLS) abundance inversion.
NB.   5. Quality assessment (reconstruction RMSE and unmixing diagnostics).
NB. ==============================================================================

cocurrent 'base'

load 'atmos/dos.ijs'
load 'subspace/hysime.ijs'
load 'unmix/vca.ijs'
load 'unmix/fcls.ijs'
load 'io/cube.ijs'

NB. Run automated unmixing pipeline on raw radiance cube
NB. Dyadic: (wvl; doy; sza_deg) run_unmixing_pipeline rad_cube
run_unmixing_pipeline=: 4 : 0
  'wvl doy sza_deg'=. x
  rad_cube=. y
  'B H W'=. $ rad_cube
  N=. H * W

  smoutput '[Albedo Pipeline] Step 1: Atmospheric DOS-1 Correction...'
  'boa_cube Lp'=. (wvl; doy; sza_deg) albedo_dos_dos1_correction rad_cube

  smoutput '[Albedo Pipeline] Step 2: HySime Subspace Dimension Estimation...'
  X_mat=. (B, N) $ boa_cube
  'dim_k costs U_hy'=. albedo_hysime_hysime_subspace_dim X_mat
  smoutput '[Albedo Pipeline]   Identified Intrinsic Dimension k = ', (": dim_k)

  smoutput '[Albedo Pipeline] Step 3: Vertex Component Analysis (VCA)...'
  'M_est endmember_idx'=. dim_k albedo_vca_extract_vca X_mat
  smoutput '[Albedo Pipeline]   Extracted ', (": dim_k), ' Endmember Vertices.'

  smoutput '[Albedo Pipeline] Step 4: Fully Constrained Least Squares (FCLS)...'
  abundances=. M_est albedo_fcls_unmix_cube boa_cube
  smoutput '[Albedo Pipeline]   Computed Abundance Fractions (ANC & ASC satisfied).'

  smoutput '[Albedo Pipeline] Step 5: Reconstruction Error Assessment...'
  rmse=. (M_est; abundances) albedo_fcls_unmix_rmse boa_cube
  smoutput '[Albedo Pipeline]   Unmixing Root Mean Square Error: ', (": rmse)

  NB. Returns pipeline result bundle
  boa_cube; M_est; abundances; dim_k; rmse
)


albedo_pipeline_unmix_run_unmixing_pipeline =: run_unmixing_pipeline
