NB. ==============================================================================
NB. Test Suite: indices/vegetation.ijs, indices/minerals_water.ijs,
NB.             classify/sam_classifier.ijs, classify/kmeans.ijs,
NB.             classify/spatial_spectral.ijs
NB. ==============================================================================

load 'tests/test_framework.ijs'
load 'indices/vegetation.ijs'
load 'indices/minerals_water.ijs'
load 'classify/sam_classifier.ijs'
load 'classify/kmeans.ijs'
load 'classify/spatial_spectral.ijs'
load 'io/spectral_lib.ijs'
load 'io/cube.ijs'

cocurrent 'base'

suite_begin 'Spectral Indices & Classification Modules'

W224=: albedo_spectral_lib_standard_wavelengths_224 ''
Veg_spec=: albedo_spectral_lib_ref_green_vegetation W224
Water_spec=: albedo_spectral_lib_ref_clear_water W224
Soil_spec=: albedo_spectral_lib_ref_kaolinite_soil W224

NB. Construct a synthetic 4x4 scene with known quadrants:
NB. Top-left: vegetation, Top-right: water, Bottom: soil
Scene_3d=: (224, 4, 4) $ 0.0
for_h. 0 1 do.
  for_w. 0 1 do. Scene_3d=: Veg_spec (<a:; h; w)} Scene_3d end.
  for_w. 2 3 do. Scene_3d=: Water_spec (<a:; h; w)} Scene_3d end.
end.
for_h. 2 3 do.
  for_w. i. 4 do. Scene_3d=: Soil_spec (<a:; h; w)} Scene_3d end.
end.

CubeScene=: albedo_cube_create_cube (Scene_3d; W224; '')

NB. Test 1: Vegetation Indices (NDVI, EVI, NDRE)
Ndvi=: albedo_indices_veg_calc_ndvi CubeScene
'NDVI shape is 4 4' assert_shape Ndvi; 4 4
'Vegetation quadrant (0,0) NDVI is strongly positive (> 0.6)' assert_true ((Ndvi < 0 0) > 0.6)
'Water quadrant (0,2) NDVI is negative (< 0.0)' assert_true ((Ndvi < 0 2) < 0.0)

Evi=: albedo_indices_veg_calc_evi CubeScene
'EVI on vegetation is positive' assert_true ((Evi < 0 0) > 0.3)

Ndre=: albedo_indices_veg_calc_ndre CubeScene
'NDRE on vegetation is positive' assert_true ((Ndre < 0 0) > 0.2)

NB. Test 2: Mineral & Water Indices (NDWI, Clay)
Ndwi=: albedo_indices_minerals_calc_ndwi CubeScene
'Water quadrant (0,2) NDWI is positive (> 0.2)' assert_true ((Ndwi < 0 2) > 0.2)
'Vegetation quadrant (0,0) NDWI is negative (< 0.0)' assert_true ((Ndwi < 0 0) < 0.0)

ClayIdx=: albedo_indices_minerals_calc_clay_index CubeScene
'Clay index on soil is higher than on vegetation' assert_true ((ClayIdx < 2 2) > (ClayIdx < 0 0))

NB. Test 3: Supervised SAM Classifier
Library_M=: |: Veg_spec ,: Water_spec ,: Soil_spec
'Labels MinAngles'=: (Library_M; 0.5) albedo_sam_classifier_classify_cube Scene_3d
'Classification map shape is 4 4' assert_shape Labels; 4 4
'Top-left pixel classified as 0 (vegetation)' assert_equal (Labels < 0 0); 0
'Top-right pixel classified as 1 (water)' assert_equal (Labels < 0 2); 1
'Bottom pixel classified as 2 (soil)' assert_equal (Labels < 2 2); 2

NB. Test rejection threshold on distant outlier
'LabelsRej MinAngRej'=: (Library_M; 0.001) albedo_sam_classifier_classify_cube Scene_3d
'Stringent rejection threshold rejects pixels as _1' assert_true (+/ , LabelsRej = _1)

NB. Test 4: Vectorized K-Means Clustering
X_data=: (224, 16) $ Scene_3d
'K_labels Centroids'=: (3; 20) albedo_kmeans_fit_kmeans X_data
'K-Means labels count is 16' assert_equal (# K_labels); 16
'Centroids shape is 224 3' assert_shape Centroids; 224 3

NB. Test 5: Spatial Majority Filter
DirtyLabels=: 5 5 $ 0
NB. Inject single noise pixel
DirtyLabels=: 1 (<2 2)} DirtyLabels
CleanedLabels=: 3 albedo_spatial_spectral_majority_filter DirtyLabels
'Majority filter removed isolated noise pixel at (2,2)' assert_equal (CleanedLabels < 2 2); 0

suite_end ''
