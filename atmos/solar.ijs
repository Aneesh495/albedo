NB. ==============================================================================
NB. Albedo - Solar Geometry & Extraterrestrial Solar Irradiance
NB.
NB. Implements Earth-Sun distance correction factor, solar declination,
NB. solar zenith angle calculation, and Thuillier 2003 reference solar irradiance
NB. spectrum E0(lambda) in W / (m^2 * nm).
NB. ==============================================================================

cocurrent 'albedo_solar'

PI=: 3.14159265358979323846
DEG2RAD=: PI % 180.0
RAD2DEG=: 180.0 % PI

NB. ------------------------------------------------------------------------------
NB. Solar Geometry & Distance Factor
NB. ------------------------------------------------------------------------------

NB. Earth-Sun distance squared factor: d^2(t)
NB. t is Day of Year (DOY in 1..365)
NB. d(t) = 1.0 - 0.01672 * cos(2 pi (DOY - 4) / 365.256)
earth_sun_dist_sq=: 3 : 0
  doy=. y
  theta=. (2.0 * PI * (doy - 4.0)) % 365.256
  d=. 1.0 - (0.01672 * 2&o. theta)
  *: d
)

NB. Solar declination angle delta (degrees) for day of year
solar_declination_deg=: 3 : 0
  doy=. y
  23.45 * 1&o. (2.0 * PI * (284.0 + doy)) % 365.0
)

NB. Solar Zenith Angle (degrees) given latitude (deg), solar declination (deg),
NB. and local solar hour angle (deg: 0 at solar noon, 15 deg per hour)
NB. cos(theta_s) = sin(lat) sin(delta) + cos(lat) cos(delta) cos(omega)
NB. Dyadic: (lat_deg; decl_deg; hour_angle_deg) solar_zenith_angle ''
solar_zenith_angle=: 3 : 0
  'lat delta omega'=. y
  lat_r=. lat * DEG2RAD
  delta_r=. delta * DEG2RAD
  omega_r=. omega * DEG2RAD

  cos_sza=. ((1&o. lat_r) * (1&o. delta_r)) + ((2&o. lat_r) * (2&o. delta_r) * (2&o. omega_r))
  cos_sza=. _1.0 >. 1.0 <. cos_sza
  _2&o. cos_sza
)

NB. Cosine of solar zenith angle
cos_solar_zenith=: 3 : 0
  sza_deg=. y
  2&o. sza_deg * DEG2RAD
)

NB. ------------------------------------------------------------------------------
NB. Extraterrestrial Solar Irradiance Spectrum E0(lambda)
NB. Thuillier 2003 reference solar spectrum model
NB. Values in W / (m^2 * nm)
NB. ------------------------------------------------------------------------------

NB. Synthetic approximation to Thuillier 2003 extraterrestrial solar irradiance
solar_irradiance_e0=: 3 : 0
  wvl=. y
  NB. Solar irradiance peaks near 500 nm (~2.0 W/(m^2*nm)) and falls off in SWIR
  NB. Modeled as Planck blackbody radiation (T_sun ~ 5778 K) normalized to solar constant ~1361 W/m^2
  w_um=. wvl * 0.001
  NB. Planck function scaled
  c1=. 1.191e4
  c2=. 1.439e4
  T=. 5778.0
  planck=. (c1 % (w_um ^ 5)) % (<: ^ (c2 % (w_um * T)))
  NB. Scale to match ~2.0 W/(m^2*nm) at 500nm (0.5 um)
  scale=. 2.0 % ((c1 % (0.5 ^ 5)) % (<: ^ (c2 % (0.5 * T))))
  0.01 >. (planck * scale)
)
