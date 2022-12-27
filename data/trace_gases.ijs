NB. ==============================================================================
NB. ALBEDO - High-Dimensional Hyperspectral Imaging Engine
NB. data/trace_gases.ijs : Atmospheric Trace Gas Cross-Sections and Plume Radiance
NB. Absorption cross-sections for CH4, CO2, CO, SO2, NH3, NO2, N2O, and H2O
NB. from 0.400 to 2.500 micrometers with Beer-Lambert plume transmission models.
NB. ==============================================================================

cocurrent 'base'

NB. 224 Wavelength centers in micrometers (0.400 to 2.500 um)
WAVELENGTHS =: 0.400 + 0.009375 * i. 224

NB. 8 Trace Gas Names
GAS_NAMES =: <;._1 ' CH4' , ' CO2' , ' CO' , ' SO2' , ' NH3' , ' NO2' , ' N2O' , ' H2O'

NB. Helper to construct characteristic absorption cross-section (in 1e_21 cm^2 / molecule):
NB. lines: list of (center_um, peak_cross_section, fwhm_um)
calc_cross_section =: 3 : 0
  lines =. y
  sigma =. 224 $ 0.0
  for_idx. i. # lines do.
    'c0 amp fwhm' =. idx { lines
    w =. (fwhm % 2) % %: 2 * ^. 2
    lorentz =. amp % 1.0 + *: (WAVELENGTHS - c0) % (fwhm % 2)
    gauss   =. amp * ^ - 0.5 * *: (WAVELENGTHS - c0) % w
    voigt   =. 0.5 * (lorentz + gauss)
    sigma =. sigma + voigt
  end.
  0.0 >. sigma
)

NB. HITRAN-derived diagnostic absorption manifolds:
NB. 1. Methane (CH4): 1.66 um (2v3 overtone), 2.20 - 2.45 um (v3 + v4 octad)
ch4_lines =: 4 3 $ 1.665 4.5 0.025  2.205 12.0 0.035  2.315 28.5 0.045  2.370 18.0 0.040

NB. 2. Carbon Dioxide (CO2): 1.57 um, 1.60 um, 2.01 um, 2.06 um
co2_lines =: 4 3 $ 1.575 3.2 0.020  1.605 4.8 0.022  2.010 32.0 0.030  2.060 26.0 0.030

NB. 3. Carbon Monoxide (CO): 2.33 um (first overtone band)
co_lines  =: 2 3 $ 2.330 14.0 0.035  2.360 8.0 0.035

NB. 4. Sulfur Dioxide (SO2): 0.400-0.420 um, 2.48 um
so2_lines =: 2 3 $ 0.405 2.5 0.025  2.480 6.0 0.040

NB. 5. Ammonia (NH3): 1.51 um, 1.98 um, 2.25 um
nh3_lines =: 3 3 $ 1.510 5.0 0.025  1.980 8.5 0.030  2.250 16.0 0.035

NB. 6. Nitrogen Dioxide (NO2): visible 0.400 to 0.500 um electronic transition
no2_lines =: 2 3 $ 0.415 15.0 0.040  0.450 10.0 0.045

NB. 7. Nitrous Oxide (N2O): 2.14 um, 2.26 um
n2o_lines =: 2 3 $ 2.140 11.0 0.030  2.260 22.0 0.035

NB. 8. Water Vapor (H2O): 0.94 um, 1.13 um, 1.40 um, 1.90 um
h2o_lines =: 4 3 $ 0.940 18.0 0.035  1.135 35.0 0.040  1.400 95.0 0.050  1.900 120.0 0.055

sigma_ch4 =: calc_cross_section ch4_lines
sigma_co2 =: calc_cross_section co2_lines
sigma_co  =: calc_cross_section co_lines
sigma_so2 =: calc_cross_section so2_lines
sigma_nh3 =: calc_cross_section nh3_lines
sigma_no2 =: calc_cross_section no2_lines
sigma_n2o =: calc_cross_section n2o_lines
sigma_h2o =: calc_cross_section h2o_lines

NB. Cross-section matrix (8 rows x 224 bands, in 1e_21 cm^2 / molecule)
SIGMA_MATRIX =: > sigma_ch4 ; sigma_co2 ; sigma_co ; sigma_so2 ; sigma_nh3 ; sigma_no2 ; sigma_n2o ; sigma_h2o

NB. Query gas index by name
get_gas_index =: 3 : 0
  idx =. GAS_NAMES i. < y
  if. idx = # GAS_NAMES do. _1 else. idx end.
)

NB. Get cross section vector by gas name
get_gas_cross_section =: 3 : 0
  idx =. get_gas_index y
  if. idx = _1 do.
    224 # 0
  else.
    idx { SIGMA_MATRIX
  end.
)

NB. Compute Beer-Lambert transmission through gas plume:
NB. y: (gas_name; column_density_ppm_m)
NB. 1 ppm-m corresponds to approx 2.45e15 molecules / cm^2 at STP.
NB. With cross section in 1e_21 cm^2 / molec:
NB. tau = sigma * 1e_21 * u_ppm_m * 2.45e15 = sigma * u_ppm_m * 2.45e_6
gas_transmittance =: 3 : 0
  'gas_name u_ppm_m' =. y
  sigma =. get_gas_cross_section gas_name
  tau =. sigma * u_ppm_m * 2.45e_6
  ^ - tau
)

NB. Linearized plume radiance absorption perturbation:
NB. delta_L(lambda) = - L_0(lambda) * sigma(lambda) * u_ppm_m * 2.45e_6
gas_radiance_perturbation =: 3 : 0
  'L_0 gas_name u_ppm_m' =. y
  sigma =. get_gas_cross_section gas_name
  - L_0 * sigma * u_ppm_m * 2.45e_6
)


NB. Exported Aliases
albedo_gases_GAS_NAMES =: GAS_NAMES
albedo_gases_WAVELENGTHS =: WAVELENGTHS
albedo_gases_SIGMA_MATRIX =: SIGMA_MATRIX
albedo_gases_get_gas_index =: get_gas_index
albedo_gases_get_gas_cross_section =: get_gas_cross_section
albedo_gases_gas_transmittance =: gas_transmittance
albedo_gases_gas_radiance_perturbation =: gas_radiance_perturbation
