NB. ==============================================================================
NB. Test Suite: atmos/solar.ijs, atmos/radiance.ijs, atmos/dos.ijs, atmos/continuum.ijs
NB. ==============================================================================

load 'tests/test_framework.ijs'
load 'atmos/solar.ijs'
load 'atmos/radiance.ijs'
load 'atmos/dos.ijs'
load 'atmos/continuum.ijs'

cocurrent 'base'

suite_begin 'Atmospheric Radiative Transfer & Calibration'

NB. Test 1: Solar Distance and Declination
NB. Perihelion around DOY 4: d^2 ~ 0.967
d2_peri=: albedo_solar_earth_sun_dist_sq 4
'Perihelion d^2 is near 0.967' assert_near d2_peri; 0.967; 0.005
NB. Aphelion around DOY 185: d^2 ~ 1.034
d2_aph=: albedo_solar_earth_sun_dist_sq 185
'Aphelion d^2 is near 1.034' assert_near d2_aph; 1.034; 0.005

NB. Equinox declination ~ 0 deg (DOY 80)
decl_eq=: albedo_solar_solar_declination_deg 80
'Equinox declination near 0 deg' assert_near decl_eq; 0.0; 1.5

NB. Test 2: Solar Zenith Angle & Extraterrestrial Irradiance
NB. At equator at equinox solar noon: SZA = 0 deg
sza_noon=: albedo_solar_solar_zenith_angle 0.0; 0.0; 0.0
'Noon SZA is 0 deg' assert_near sza_noon; 0.0; 1e_4
cos_noon=: albedo_solar_cos_solar_zenith sza_noon
'cos(SZA) at noon is 1.0' assert_near cos_noon; 1.0; 1e_4

E0_500=: albedo_solar_solar_irradiance_e0 500.0
'E0 at 500 nm near 2.0 W/(m^2*nm)' assert_near E0_500; 2.0; 0.1
E0_2000=: albedo_solar_solar_irradiance_e0 2000.0
'E0 at 2000 nm is lower than at 500 nm' assert_true (E0_2000 < E0_500)

NB. Test 3: Radiometric DN to Radiance & TOA Reflectance
RawDN=: 2 4 4 $ 1000
Gains=: 0.05 0.08
Offsets=: _1.0 _2.0
RadCube=: (Gains; Offsets) albedo_radiance_dn_to_radiance RawDN
'Radiance band 0 value is 49.0' assert_near (RadCube < 0 0 0); 49.0; 1e_5
'Radiance band 1 value is 78.0' assert_near (RadCube < 1 0 0); 78.0; 1e_5

NB. Convert to TOA Reflectance
Wvls=: 550.0 850.0
ToaCube=: (Wvls; 80; 30.0) albedo_radiance_radiance_to_toa_reflectance RadCube
'TOA cube shape is 2 4 4' assert_shape ToaCube; 2 4 4
'TOA values within valid range [0, 1]' assert_true (*./ (*./ , ToaCube >: 0.0) *. (*./ , ToaCube <: 1.0))

NB. Test 4: Atmospheric Bad Bands Masking
W_all=: 1000.0 1400.0 1600.0 1900.0 2200.0
Mask=: albedo_radiance_get_valid_band_mask W_all
'Band 1000 nm is valid' assert_true (0 { Mask)
'Band 1400 nm (water vapor) is invalid' assert_false (1 { Mask)
'Band 1600 nm is valid' assert_true (2 { Mask)
'Band 1900 nm (water vapor) is invalid' assert_false (3 { Mask)
'Band 2200 nm is valid' assert_true (4 { Mask)

NB. Test 5: Dark Object Subtraction (DOS-1)
Lp_est=: albedo_dos_estimate_path_radiance RadCube
'Path radiance count matches bands' assert_equal (# Lp_est); 2
'Path radiance >= 0' assert_true (*./ Lp_est >: 0.0)
'BOACube Lp'=: (Wvls; 80; 30.0) albedo_dos_dos1_correction RadCube
'BOA cube shape is 2 4 4' assert_shape BOACube; 2 4 4

NB. Test 6: Convex Hull Continuum Removal
NB. Synthetic spectrum with Gaussian absorption trough at 2200 nm
W_scan=: 2000.0 + (i. 50) * 10.0
Refl_trough=: 0.5 - 0.25 * ^ - *: (W_scan - 2200.0) % 30.0
Hull=: W_scan albedo_continuum_upper_convex_hull Refl_trough
'Hull is always >= original reflectance' assert_true (*./ Hull >: (Refl_trough - 1e_7))

CR=: W_scan albedo_continuum_remove_continuum Refl_trough
'Continuum removed values in [0, 1]' assert_true (*./ (CR >: 0.0) *. (CR <: (1.0 + 1e_6)))
'Minimum CR is near 2200 nm' assert_near (<./ CR); 0.5; 0.05

'center max_d area fwhd'=: (2100.0 2300.0) albedo_continuum_characterize_absorption (W_scan; Refl_trough)
'Absorption center near 2200 nm' assert_near center; 2200.0; 10.0
'Max absorption depth near 0.5' assert_near max_d; 0.5; 0.05
'Absorption area > 0' assert_true (area > 0.0)
'FWHD > 0' assert_true (fwhd > 0.0)

suite_end ''
