NB. ==============================================================================
NB. Albedo - Radiometric Calibration & TOA Reflectance
NB.
NB. Converts raw sensor Digital Numbers (DN) to at-sensor spectral radiance
NB. L(lambda), and converts spectral radiance to Top-of-Atmosphere (TOA)
NB. reflectance rho_TOA(lambda). Identifies bad atmospheric absorption bands.
NB. ==============================================================================

cocurrent 'base'

load 'atmos/solar.ijs'

PI=: 3.14159265358979323846

NB. ------------------------------------------------------------------------------
NB. Digital Number (DN) to Radiance Calibration
NB. L(lambda) = gain(lambda) * DN + offset(lambda)
NB. ------------------------------------------------------------------------------

NB. Calibrate 3D raw DN cube (B x H x W) to Radiance (W / (m^2 * sr * nm))
NB. Dyadic: (gain_vec; offset_vec) dn_to_radiance dn_cube
dn_to_radiance=: 4 : 0
  'gain offset'=. x
  dn_cube=. y
  'B H W'=. $ dn_cube
  assert. B = # gain [ 'Gain vector length must match band count'
  assert. B = # offset [ 'Offset vector length must match band count'

  rad_cube=. (B, H, W) $ 0.0
  for_b. i. B do.
    g=. b { gain
    o=. b { offset
    rad_cube=. (o + g * b { dn_cube) (<b)} rad_cube
  end.
  rad_cube
)

NB. ------------------------------------------------------------------------------
NB. Radiance to Top-of-Atmosphere (TOA) Reflectance
NB. rho_TOA = (pi * L * d^2) / (E0 * cos(theta_s))
NB. ------------------------------------------------------------------------------

NB. Converts 3D Radiance cube to TOA Reflectance cube
NB. Dyadic: (wvl; doy; sza_deg) radiance_to_toa_reflectance rad_cube
radiance_to_toa_reflectance=: 4 : 0
  'wvl doy sza_deg'=. x
  rad_cube=. y
  'B H W'=. $ rad_cube

  d2=. albedo_solar_earth_sun_dist_sq doy
  cos_sza=. albedo_solar_cos_solar_zenith sza_deg
  assert. cos_sza > 0.01 [ 'Solar zenith angle too steep for daytime reflectance'

  E0=. albedo_solar_solar_irradiance_e0 wvl

  refl_cube=. (B, H, W) $ 0.0
  for_b. i. B do.
    denom=. (b { E0) * cos_sza
    scale=. (PI * d2) % denom
    refl_band=. scale * (b { rad_cube)
    refl_cube=. (0.0 >. 1.5 <. refl_band) (<b)} refl_cube
  end.
  refl_cube
)

NB. ------------------------------------------------------------------------------
NB. Atmospheric Bad Bands Masking
NB. Opaque atmospheric water vapor absorption windows:
NB. 1350 - 1450 nm and 1800 - 1950 nm
NB. ------------------------------------------------------------------------------

NB. Identifies indices of valid atmospheric transmission bands
get_valid_band_mask=: 3 : 0
  wvl=. y
  NB. Atmospheric absorption windows to exclude
  water1=. (wvl >: 1350.0) *. (wvl <: 1440.0)
  water2=. (wvl >: 1800.0) *. (wvl <: 1940.0)
  bad=. water1 +. water2
  0 = bad
)

NB. Filter cube to only valid atmospheric transmission bands
filter_valid_bands=: 3 : 0
  'cube wvl meta'=. y
  mask=. get_valid_band_mask wvl
  valid_indices=. I. mask
  sub_cube=. valid_indices { cube
  sub_wvl=. valid_indices { wvl
  sub_cube; sub_wvl; meta
)


albedo_radiance_dn_to_radiance =: dn_to_radiance
albedo_radiance_radiance_to_toa_reflectance =: radiance_to_toa_reflectance
albedo_radiance_get_valid_band_mask =: get_valid_band_mask
