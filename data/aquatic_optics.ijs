NB. ==============================================================================
NB. ALBEDO - High-Dimensional Hyperspectral Imaging Engine
NB. data/aquatic_optics.ijs : Bio-Optical Ocean and Inland Water Reflectance Models
NB. Implements Case-1 (open ocean) and Case-2 (turbid/coastal/estuarine) water models
NB. with pure seawater absorption, phytoplankton pigment, CDOM, and sediment scattering.
NB. ==============================================================================

cocurrent 'albedo_aquatic'

NB. 224 Wavelength centers in micrometers (0.400 to 2.500 um)
WAVELENGTHS =: 0.400 + 0.009375 * i. 224

NB. Pure water absorption coefficient a_w (m^-1) across 0.400 to 2.500 um
NB. Smith & Baker / Pope & Fry empirical parameterization
calc_aw =: 3 : 0
  NB. Base exponential rise from blue minimum (0.005 m^-1 at 0.420 um) to red/SWIR
  w =. WAVELENGTHS
  NB. Visible regime (0.400 to 0.700 um)
  vis_aw =. 0.005 + 0.03 * *: (w - 0.420) % 0.150
  NB. NIR steps: 0.750 um (2.5 m^-1), 0.970 um (35 m^-1)
  nir_step =. 2.5 % 1.0 + ^ - (w - 0.720) % 0.030
  nir_peak =. 35.0 * ^. - 0.5 * *: (w - 0.970) % 0.050
  NB. SWIR regime: heavy absorption > 1.2 um (100 to 10000 m^-1)
  swir_step =. 250.0 % 1.0 + ^ - (w - 1.150) % 0.080
  swir_huge =. 1200.0 % 1.0 + ^ - (w - 1.400) % 0.100
  vis_aw + nir_step + nir_peak + swir_step + swir_huge
)

NB. Pure water backscattering coefficient b_bw (m^-1): lambda^-4.32 power law
calc_bbw =: 3 : 0
  0.00193 * (0.500 % WAVELENGTHS) ^ 4.32
)

NB. Phytoplankton specific absorption a_star_ph (m^2 / mg Chl-a):
NB. Diagnostic blue peak at 0.440 um, red peak at 0.675 um, fluorescence at 0.685 um
calc_aph =: 3 : 0
  w =. WAVELENGTHS
  blue_peak =. 0.045 * ^. - 0.5 * *: (w - 0.440) % 0.035
  red_peak  =. 0.025 * ^. - 0.5 * *: (w - 0.675) % 0.020
  blue_peak + red_peak
)

NB. Colored Dissolved Organic Matter (CDOM) absorption:
NB. a_cdom(lambda) = a_cdom(0.440) * exp(- S * (lambda - 0.440))
NB. where S is spectral slope factor (~15.0 um^-1)
calc_acdom =: 3 : 0
  'a440 s_slope' =. y
  a440 * ^ - s_slope * (WAVELENGTHS - 0.440)
)

NB. Non-Algal Particulate (NAP) / suspended sediment backscattering:
NB. b_bp(lambda) = b_bp(0.550) * (0.550 / lambda)^eta
calc_bbp =: 3 : 0
  'bbp550 eta' =. y
  bbp550 * (0.550 % WAVELENGTHS) ^ eta
)

A_W   =: calc_aw 0
BB_W  =: calc_bbw 0
A_STAR_PH =: calc_aph 0

NB. Forward bio-optical model: compute subsurface remote-sensing reflectance r_rs(lambda)
NB. and above-water remote-sensing reflectance R_rs(lambda)
NB. Input y: (chl_mg_m3 ; a_cdom_440 ; bbp_550)
forward_water_reflectance =: 3 : 0
  'chl a_cdom bbp' =. y
  NB. Total absorption: a = a_w + a_ph + a_cdom
  a_ph =. chl * A_STAR_PH
  a_cd =. calc_acdom (a_cdom ; 16.0)
  a_tot =. A_W + a_ph + a_cd

  NB. Total backscattering: bb = bb_w + bb_p
  bb_p =. calc_bbp (bbp ; 0.8)
  bb_tot =. BB_W + bb_p

  NB. Gordon subsurface model: r_rs = (g0 + g1 * u) * u, where u = bb / (a + bb)
  u =. bb_tot % a_tot + bb_tot
  r_rs =. (0.089 + 0.125 * u) * u

  NB. Air-water interface transfer: R_rs = 0.52 * r_rs / (1 - 1.7 * r_rs)
  R_rs =. (0.52 * r_rs) % 1.0 - (1.7 * r_rs)
  0.0 >. R_rs
)

NB. Case 1 Open Ocean Clear Water (Low Chl=0.05, Low CDOM=0.01, Low Sed=0.001)
clear_ocean =: forward_water_reflectance 0.05 0.01 0.001

NB. Eutrophic Algal Bloom Water (High Chl=15.0, CDOM=0.08, Sed=0.005)
algal_bloom =: forward_water_reflectance 15.0 0.08 0.005

NB. Turbid Coastal / Estuarine Water (High Sed=0.05, Moderate Chl=2.0, High CDOM=0.15)
turbid_coastal =: forward_water_reflectance 2.0 0.15 0.05

NB. Blackwater / Humic Lake (Extreme CDOM=0.50, Chl=1.0, Sed=0.002)
humic_lake =: forward_water_reflectance 1.0 0.50 0.002

NB. Shape: 4 rows by 224 columns
AQUATIC_MATRIX =: > clear_ocean ; algal_bloom ; turbid_coastal ; humic_lake
AQUATIC_NAMES  =: <;._1 ' ClearOcean' , ' AlgalBloom' , ' TurbidCoastal' , ' HumicLake'
