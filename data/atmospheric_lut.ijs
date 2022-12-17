NB. ==============================================================================
NB. Albedo - Atmospheric Radiative Transfer Look-Up Table (LUT)
NB.
NB. Physical atmospheric optical transmission spectra, gaseous absorption
NB. cross-sections, Rayleigh optical depth, and aerosol extinction models
NB. across 224 spectral channels (400 - 2500 nm) for 6 standard atmospheres:
NB.   1. US Standard 1976
NB.   2. Tropical (high columnar water vapor 4.1 g/cm^2)
NB.   3. Midlatitude Summer (columnar water vapor 2.9 g/cm^2)
NB.   4. Midlatitude Winter (columnar water vapor 0.85 g/cm^2)
NB.   5. Subarctic Summer (columnar water vapor 2.1 g/cm^2)
NB.   6. Subarctic Winter (dry atmosphere 0.42 g/cm^2)
NB. ==============================================================================

cocurrent 'base'

load 'io/spectral_lib.ijs'

W224=: albedo_spectral_lib_standard_wavelengths_224 ''

NB. ------------------------------------------------------------------------------
NB. Rayleigh Optical Depth Model (Bodhaine et al., 1999)
NB. tau_R(lambda) = 0.008569 * lambda^(-4) * (1 + 0.0113 * lambda^(-2) + 0.00013 * lambda^(-4))
NB. where lambda is in micrometers (um)
NB. ------------------------------------------------------------------------------

rayleigh_optical_depth=: 3 : 0
  wvl_nm=. y
  w_um=. wvl_nm * 0.001
  term1=. 0.008569 * (w_um ^ _4.0)
  term2=. 1.0 + (0.0113 * (w_um ^ _2.0)) + (0.00013 * (w_um ^ _4.0))
  term1 * term2
)

NB. ------------------------------------------------------------------------------
NB. Gaseous Absorption Models
NB. ------------------------------------------------------------------------------

NB. Atmospheric Ozone (O3) Chappuis absorption band in visible (500 - 700 nm)
ozone_optical_depth=: 3 : 0
  wvl=. y
  NB. Peak around 600 nm with bandwidth ~100 nm
  0.035 * ^ - *: (wvl - 600.0) % 55.0
)

NB. Atmospheric Oxygen (O2) A-band (760 nm) and B-band (690 nm)
oxygen_optical_depth=: 3 : 0
  wvl=. y
  tau_a=. 0.85 * ^ - *: (wvl - 760.5) % 6.0
  tau_b=. 0.15 * ^ - *: (wvl - 688.0) % 5.0
  tau_a + tau_b
)

NB. Carbon Dioxide (CO2) absorption in SWIR-2 (2010 nm & 2060 nm)
carbon_dioxide_optical_depth=: 3 : 0
  wvl=. y
  tau_co2_1=. 0.45 * ^ - *: (wvl - 2010.0) % 12.0
  tau_co2_2=. 0.55 * ^ - *: (wvl - 2060.0) % 14.0
  tau_co2_1 + tau_co2_2
)

NB. Water Vapor (H2O) absorption cross-sections (720, 820, 940, 1130, 1400, 1900 nm)
water_vapor_absorption_profile=: 3 : 0
  wvl=. y
  tau=. wvl * 0.0
  tau=. tau + 0.08 * ^ - *: (wvl - 720.0) % 10.0
  tau=. tau + 0.15 * ^ - *: (wvl - 820.0) % 15.0
  tau=. tau + 0.65 * ^ - *: (wvl - 940.0) % 25.0
  tau=. tau + 0.75 * ^ - *: (wvl - 1130.0) % 30.0
  tau=. tau + 4.50 * ^ - *: (wvl - 1400.0) % 45.0
  tau=. tau + 5.20 * ^ - *: (wvl - 1900.0) % 55.0
  tau
)

NB. Aerosol Extinction Optical Depth (Angstrom law: tau_aer = beta * lambda^(-alpha))
NB. Continental aerosol: alpha ~ 1.3, beta ~ 0.08 at 550 nm
aerosol_optical_depth=: 3 : 0
  wvl_nm=. y
  w_um=. wvl_nm * 0.001
  0.08 * (w_um % 0.55) ^ _1.3
)

NB. ------------------------------------------------------------------------------
NB. Total Atmospheric Direct Transmittance
NB. T(lambda) = exp(-(tau_R + tau_aer + tau_O3 + tau_O2 + tau_CO2 + cw * tau_H2O) / cos(theta_s))
NB. ------------------------------------------------------------------------------

NB. Compute total atmospheric direct transmission spectrum
NB. Dyadic: (sza_deg; columnar_water_vapor_scale) atmospheric_transmission wvl
atmospheric_transmission=: 4 : 0
  'sza_deg cw_scale'=. x
  wvl=. y
  cos_sza=. 2&o. sza_deg * (3.141592653589793 % 180.0)
  cos_sza=. cos_sza >. 0.05

  tau_R=.   rayleigh_optical_depth wvl
  tau_aer=. aerosol_optical_depth wvl
  tau_o3=.  ozone_optical_depth wvl
  tau_o2=.  oxygen_optical_depth wvl
  tau_co2=. carbon_dioxide_optical_depth wvl
  tau_h2o=. cw_scale * water_vapor_absorption_profile wvl

  tau_total=. tau_R + tau_aer + tau_o3 + tau_o2 + tau_co2 + tau_h2o
  0.001 >. 1.0 <. ^ - (tau_total % cos_sza)
)

NB. ------------------------------------------------------------------------------
NB. Standard Atmosphere Models (Water Vapor Column Scalers)
NB. ------------------------------------------------------------------------------

NB. US Standard 1976 (cw = 1.42 g/cm^2)
trans_us_standard_1976=: 3 : 0
  (30.0; 1.0) atmospheric_transmission y
)

NB. Tropical (cw = 4.11 g/cm^2 -> scale = 2.89)
trans_tropical=: 3 : 0
  (30.0; 2.89) atmospheric_transmission y
)

NB. Midlatitude Summer (cw = 2.92 g/cm^2 -> scale = 2.05)
trans_midlatitude_summer=: 3 : 0
  (30.0; 2.05) atmospheric_transmission y
)

NB. Midlatitude Winter (cw = 0.85 g/cm^2 -> scale = 0.60)
trans_midlatitude_winter=: 3 : 0
  (30.0; 0.60) atmospheric_transmission y
)

NB. Subarctic Summer (cw = 2.08 g/cm^2 -> scale = 1.46)
trans_subarctic_summer=: 3 : 0
  (30.0; 1.46) atmospheric_transmission y
)

NB. Subarctic Winter (cw = 0.42 g/cm^2 -> scale = 0.30)
trans_subarctic_winter=: 3 : 0
  (30.0; 0.30) atmospheric_transmission y
)


albedo_atmos_lut_W224 =: W224
albedo_atmos_lut_rayleigh_optical_depth =: rayleigh_optical_depth
albedo_atmos_lut_oxygen_optical_depth =: oxygen_optical_depth
albedo_atmos_lut_trans_tropical =: trans_tropical
albedo_atmos_lut_trans_us_standard_1976 =: trans_us_standard_1976
albedo_atmos_lut_trans_subarctic_winter =: trans_subarctic_winter
