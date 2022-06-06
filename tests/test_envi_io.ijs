NB. ==============================================================================
NB. Test Suite: io/envi.ijs, io/spectral_lib.ijs, io/cube.ijs
NB. ==============================================================================

load 'tests/test_framework.ijs'
load 'io/envi.ijs'
load 'io/spectral_lib.ijs'
load 'io/cube.ijs'

cocurrent 'base'

suite_begin 'ENVI I/O, Spectral Library & Cube Slicing'

NB. Test 1: Interleave Conversions
Cube_bsq=: i. 4 3 2
Cube_bil=: albedo_envi_bsq_to_bil Cube_bsq
'BIL shape should be 3 4 2' assert_shape Cube_bil; 3 4 2
Cube_bsq_rec1=: albedo_envi_bil_to_bsq Cube_bil
'BIL roundtrip matches BSQ' assert_equal Cube_bsq_rec1; Cube_bsq

Cube_bip=: albedo_envi_bsq_to_bip Cube_bsq
'BIP shape should be 3 2 4' assert_shape Cube_bip; 3 2 4
Cube_bsq_rec2=: albedo_envi_bip_to_bsq Cube_bip
'BIP roundtrip matches BSQ' assert_equal Cube_bsq_rec2; Cube_bsq

NB. Test 2: ENVI Header Parsing
HdrText=: 'ENVI', LF, 'samples = 64', LF, 'lines = 64', LF, 'bands = 10', LF, 'interleave = bsq', LF, 'data type = 4', LF
ParsedHdr=: albedo_envi_parse_envi_header HdrText
'Header samples lookup' assert_equal ('samples' albedo_envi_get_header_field ParsedHdr); '64'
'Header lines lookup' assert_equal ('lines' albedo_envi_get_header_field ParsedHdr); '64'
'Header bands lookup' assert_equal ('bands' albedo_envi_get_header_field ParsedHdr); '10'
'Header interleave lookup' assert_equal ('interleave' albedo_envi_get_header_field ParsedHdr); 'bsq'

NB. Test 3: ENVI Header Generation
GenHdr=: albedo_envi_generate_envi_header (64; 64; 10; 4; 'bsq'; 0; (i. 10); (10 $ 10.0))
'Generated header starts with ENVI' assert_true ('ENVI' -: 4 {. GenHdr)
'Generated header contains samples = 64' assert_true (+/ 'samples = 64' E. GenHdr)

NB. Test 4: Spectral Library Resampling
SrcWvl=: 400.0 500.0 600.0 700.0
SrcRefl=: 0.1 0.2 0.3 0.4
TargetWvl=: 450.0 550.0 650.0
Resampled=: TargetWvl albedo_spectral_lib_resample_spectrum (SrcWvl; SrcRefl)
'Resampled length is 3' assert_equal (# Resampled); 3
'Resampled at 450 is 0.15' assert_near (0 { Resampled); 0.15; 1e_6
'Resampled at 550 is 0.25' assert_near (1 { Resampled); 0.25; 1e_6
'Resampled at 650 is 0.35' assert_near (2 { Resampled); 0.35; 1e_6

NB. Test 5: Standard Endmember Library Matrix
W224=: albedo_spectral_lib_standard_wavelengths_224 ''
'W224 length is 224' assert_equal (# W224); 224
StdLib=: albedo_spectral_lib_build_standard_library W224
'Standard library shape should be 224 5' assert_shape StdLib; 224 5
'Reflectance values strictly in [0, 1]' assert_true (*./ (*./ StdLib >: 0.0) *. (*./ StdLib <: 1.0))

NB. Test 6: Cube Container & Slicing
TestCubeData=: ? 224 16 16 $ 1000
CubeBox=: albedo_cube_create_cube (TestCubeData; W224; ParsedHdr)
'Cube shape is 224 16 16' assert_shape (albedo_cube_cube_data CubeBox); 224 16 16

NB. Spatial slice [0, 8, 0, 8]
SpatialSub=: (0 8 0 8) albedo_cube_slice_spatial CubeBox
'Spatial subcube shape is 224 8 8' assert_shape (albedo_cube_cube_data SpatialSub); 224 8 8

NB. Spectral slice [500, 900 nm]
SpecSub=: (500.0 900.0) albedo_cube_slice_spectral CubeBox
'Spectral subcube bands within range' assert_true (*./ ((albedo_cube_cube_wvl SpecSub) >: 500.0) *. ((albedo_cube_cube_wvl SpecSub) <: 900.0))

NB. False Color RGB Composite
RGB_cir=: albedo_cube_color_infrared_cir CubeBox
'CIR RGB shape is 3 16 16' assert_shape RGB_cir; 3 16 16
'CIR RGB values normalized to [0, 1]' assert_true (*./ (*./ , RGB_cir >: 0.0) *. (*./ , RGB_cir <: 1.0))

suite_end ''
