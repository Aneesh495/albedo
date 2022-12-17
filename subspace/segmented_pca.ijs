NB. ==============================================================================
NB. Albedo - Segmented Spectral Principal Component Analysis
NB.
NB. Divides continuous spectral channels into physically distinct optical
NB. windows (VNIR: 400-1000 nm, SWIR-1: 1000-1800 nm, SWIR-2: 2000-2500 nm)
NB. and executes partitioned PCA within each segment to prevent dominant
NB. VNIR vegetation signals from obscuring diagnostic SWIR mineral absorption features.
NB. ==============================================================================

cocurrent 'base'

load 'subspace/pca.ijs'

NB. ------------------------------------------------------------------------------
NB. Segment Definition & Partitioning
NB. ------------------------------------------------------------------------------

NB. Partition wavelength vector into 3 segments:
NB. Segment 0 (VNIR):   400 - 1000 nm
NB. Segment 1 (SWIR-1): 1000 - 1800 nm
NB. Segment 2 (SWIR-2): 1800 - 2500 nm
NB. Returns boxed list of band index vectors: (idx_vnir; idx_swir1; idx_swir2)
partition_bands=: 3 : 0
  wvl=. y
  idx_vnir=.  I. (wvl >: 400.0)  *. (wvl < 1000.0)
  idx_swir1=. I. (wvl >: 1000.0) *. (wvl < 1800.0)
  idx_swir2=. I. (wvl >: 1800.0) *. (wvl <: 2500.0)
  idx_vnir; idx_swir1; idx_swir2
)

NB. ------------------------------------------------------------------------------
NB. Segmented PCA Training & Transform
NB. ------------------------------------------------------------------------------

NB. Train segmented PCA model
NB. Dyadic: (k_per_segment; wvl) train_segmented_pca X_mat
train_segmented_pca=: 4 : 0
  'k_list wvl'=. x
  X_mat=. y
  'idx_vnir idx_swir1 idx_swir2'=. partition_bands wvl

  models=. 0 $ <''
  segments=. idx_vnir; idx_swir1; idx_swir2

  for_s. i. 3 do.
    idx_s=. s {:: segments
    if. 0 < # idx_s do.
      X_s=. idx_s { X_mat
      model_s=. albedo_pca_train_pca X_s
      models=. models , < model_s
    else.
      models=. models , < ''
    end.
  end.

  k_list; segments; models
)

NB. Forward Segmented PCA Projection
NB. Dyadic: seg_pca_model project_segmented_pca X_mat
project_segmented_pca=: 4 : 0
  'k_list segments models'=. x
  X_mat=. y

  proj_blocks=. 0 $ <''
  for_s. i. 3 do.
    idx_s=. s {:: segments
    k_s=. s { k_list
    model_s=. s {:: models
    if. (0 < # idx_s) *. (0 < # model_s) do.
      X_s=. idx_s { X_mat
      Y_s=. (k_s; model_s) albedo_pca_project_pca X_s
      proj_blocks=. proj_blocks , < Y_s
    end.
  end.

  NB. Vertically concatenate projected feature blocks
  ; proj_blocks
)


albedo_segmented_pca_partition_bands =: partition_bands
albedo_segmented_pca_train_segmented_pca =: train_segmented_pca
albedo_segmented_pca_project_segmented_pca =: project_segmented_pca
