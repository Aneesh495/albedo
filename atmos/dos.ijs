NB. ==============================================================================
NB. Albedo - Dark Object Subtraction (DOS) Atmospheric Haze Correction
NB.
NB. Implements DOS-1, DOS-2, and Chavez COST atmospheric correction algorithms
NB. to estimate path radiance Lp(lambda) from shadow/dark pixels and retrieve
NB. Bottom-of-Atmosphere (BOA) surface reflectance.
NB. ==============================================================================

cocurrent 'base'

load 'atmos/solar.ijs'

PI=: 3.14159265358979323846

NB. ------------------------------------------------------------------------------
NB. Dark Object Detection & Path Radiance Estimation
NB. ------------------------------------------------------------------------------

NB. Estimate path radiance Lp(lambda) for each band as the 1st percentile of radiance
estimate_path_radiance=: 3 : 0
  rad_cube=. y
  'B H W'=. $ rad_cube
  Lp=. B $ 0.0

  for_b. i. B do.
    band=. b { rad_cube
    vals=. , band
    s_vals=. /:~ vals
    n=. # s_vals
    NB. Take 1% percentile as dark object radiance
    p01=. (<. 0.01 * n) { s_vals
    Lp=. (0.0 >. p01) (<b)} Lp
  end.
  Lp
)

NB. ------------------------------------------------------------------------------
NB. DOS-1 Surface Reflectance Retrieval
NB. rho_BOA = pi * (L - Lp) * d^2 / (E0 * cos(theta_s))
NB. ------------------------------------------------------------------------------

dos1_correction=: 4 : 0
  'wvl doy sza_deg'=. x
  rad_cube=. y
  'B H W'=. $ rad_cube

  d2=. albedo_solar_earth_sun_dist_sq doy
  cos_sza=. albedo_solar_cos_solar_zenith sza_deg
  E0=. albedo_solar_solar_irradiance_e0 wvl
  Lp=. estimate_path_radiance rad_cube

  boa_cube=. (B, H, W) $ 0.0
  for_b. i. B do.
    lp_val=. b { Lp
    e0_val=. b { E0
    scale=. (PI * d2) % (e0_val * cos_sza)
    band_clean=. scale * (0.0 >. (b { rad_cube) - lp_val)
    boa_cube=. (0.0 >. 1.0 <. band_clean) (<b)} boa_cube
  end.
  boa_cube; Lp
)

NB. ------------------------------------------------------------------------------
NB. Chavez COST Atmospheric Model
NB. Incorporates cosine of solar zenith angle for atmospheric transmittance:
NB. tau_v = 1.0, tau_z = cos(theta_s)
NB. rho_BOA = pi * (L - Lp) * d^2 / (E0 * cos^2(theta_s))
NB. ------------------------------------------------------------------------------

cost_correction=: 4 : 0
  'wvl doy sza_deg'=. x
  rad_cube=. y
  'B H W'=. $ rad_cube

  d2=. albedo_solar_earth_sun_dist_sq doy
  cos_sza=. albedo_solar_cos_solar_zenith sza_deg
  cos2_sza=. *: cos_sza
  assert. cos2_sza > 0.001 [ 'Solar zenith angle too steep for COST model'

  E0=. albedo_solar_solar_irradiance_e0 wvl
  Lp=. estimate_path_radiance rad_cube

  boa_cube=. (B, H, W) $ 0.0
  for_b. i. B do.
    lp_val=. b { Lp
    e0_val=. b { E0
    scale=. (PI * d2) % (e0_val * cos2_sza)
    band_clean=. scale * (0.0 >. (b { rad_cube) - lp_val)
    boa_cube=. (0.0 >. 1.0 <. band_clean) (<b)} boa_cube
  end.
  boa_cube; Lp
)


albedo_dos_estimate_path_radiance =: estimate_path_radiance
albedo_dos_dos1_correction =: dos1_correction
