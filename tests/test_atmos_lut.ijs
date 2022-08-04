NB. ==============================================================================
NB. Test Suite: data/atmospheric_lut.ijs
NB. ==============================================================================

load 'tests/test_framework.ijs'
load 'data/atmospheric_lut.ijs'

cocurrent 'base'

suite_begin 'Atmospheric Radiative Transfer Look-Up Table (LUT)'

W224=: albedo_atmos_lut_W224

NB. Test 1: Rayleigh Optical Depth
tau_r_400=: albedo_atmos_lut_rayleigh_optical_depth 400.0
tau_r_2000=: albedo_atmos_lut_rayleigh_optical_depth 2000.0
'Rayleigh tau at 400 nm is near 0.35' assert_near tau_r_400; 0.35; 0.05
'Rayleigh tau at 2000 nm is near 0' assert_near tau_r_2000; 0.001; 0.005
'Rayleigh decreases monotonically with wavelength' assert_true (tau_r_400 > tau_r_2000)

NB. Test 2: Oxygen A-Band Absorption
tau_o2_760=: albedo_atmos_lut_oxygen_optical_depth 760.5
tau_o2_700=: albedo_atmos_lut_oxygen_optical_depth 700.0
'Oxygen A-band peak at 760.5 nm' assert_true (tau_o2_760 > 0.7)
'Oxygen absorption off-peak is low' assert_true (tau_o2_700 < 0.1)

NB. Test 3: Standard Atmosphere Transmissions
T_us=: albedo_atmos_lut_trans_us_standard_1976 W224
'US Standard transmittance length is 224' assert_equal (# T_us); 224
'All transmittance values bounded in (0, 1]' assert_true (*./ (T_us > 0.0) *. (T_us <: 1.0))

NB. Test 4: Tropical vs Subarctic Winter Water Vapor Opacity
T_tropical=: albedo_atmos_lut_trans_tropical W224
T_subarctic_w=: albedo_atmos_lut_trans_subarctic_winter W224

idx_940=: 940.0 albedo_cube_find_nearest_band W224
idx_1400=: 1400.0 albedo_cube_find_nearest_band W224

'Tropical has more water vapor absorption (lower transmission) at 940 nm than subarctic winter' assert_true ((idx_940 { T_tropical) < (idx_940 { T_subarctic_w))
'Both models exhibit near total extinction at 1400 nm (< 0.05)' assert_true ((idx_1400 { T_tropical) < 0.05)

suite_end ''
