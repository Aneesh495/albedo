NB. ==============================================================================
NB. Test Suite: data/usgs_minerals.ijs, data/usgs_vegetation.ijs, data/usgs_urban.ijs
NB. ==============================================================================

load 'tests/test_framework.ijs'
load 'io/cube.ijs'
load 'data/usgs_minerals.ijs'
load 'data/usgs_vegetation.ijs'
load 'data/usgs_urban.ijs'

cocurrent 'base'

suite_begin 'Calibrated USGS Spectral Libraries (Minerals, Flora, Urban)'

W224=: albedo_usgs_minerals_W224

NB. Test 1: Mineral Library
MineralLib=: albedo_usgs_minerals_build_usgs_mineral_library W224
'Mineral library shape is 224 15' assert_shape MineralLib; 224 15
'All mineral reflectance values in [0, 1]' assert_true (*./ (*./ MineralLib >: 0.0) *. (*./ MineralLib <: 1.0))

NB. Verify Calcite absorption at 2340 nm
CalciteSpec=: 3 {"1 MineralLib
idx_2340=: 2340.0 albedo_cube_find_nearest_band W224
idx_2100=: 2100.0 albedo_cube_find_nearest_band W224
'Calcite has strong absorption trough at 2340 nm' assert_true ((idx_2340 { CalciteSpec) < (idx_2100 { CalciteSpec))

NB. Verify Kaolinite Al-OH doublet at 2205 nm
KaoliniteSpec=: 0 {"1 MineralLib
idx_2205=: 2205.0 albedo_cube_find_nearest_band W224
'Kaolinite has Al-OH doublet absorption at 2205 nm' assert_true ((idx_2205 { KaoliniteSpec) < (idx_2100 { KaoliniteSpec))

NB. Test 2: Botanical Library
FloraLib=: albedo_usgs_vegetation_build_botanical_library W224
'Botanical library shape is 224 10' assert_shape FloraLib; 224 10
'All botanical reflectance values in [0, 1]' assert_true (*./ (*./ FloraLib >: 0.0) *. (*./ FloraLib <: 1.0))

NB. Verify vegetation red edge jump
DecidSpec=: 0 {"1 FloraLib
idx_red=: 670.0 albedo_cube_find_nearest_band W224
idx_nir=: 800.0 albedo_cube_find_nearest_band W224
RedEdgeJump=: (idx_nir { DecidSpec) - (idx_red { DecidSpec)
'Broadleaf exhibits sharp red edge jump (> 0.3)' assert_true (RedEdgeJump > 0.3)

NB. Test 3: Urban & Anthropogenic Materials Library
UrbanLib=: albedo_usgs_urban_build_urban_library W224
'Urban library shape is 224 10' assert_shape UrbanLib; 224 10
'All urban reflectance values in [0, 1]' assert_true (*./ (*./ UrbanLib >: 0.0) *. (*./ UrbanLib <: 1.0))

NB. Verify Polyethylene plastic C-H absorption at 1730 nm
PeSpec=: 6 {"1 UrbanLib
idx_1730=: 1730.0 albedo_cube_find_nearest_band W224
idx_1600=: 1600.0 albedo_cube_find_nearest_band W224
'Polyethylene plastic exhibits diagnostic C-H dip at 1730 nm' assert_true ((idx_1730 { PeSpec) < (idx_1600 { PeSpec))

suite_end ''
