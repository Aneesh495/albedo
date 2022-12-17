NB. ==============================================================================
NB. ALBEDO - Tests: data/trace_gases.ijs
NB. ==============================================================================

load 'data/trace_gases.ijs'

assert =: 3 : 0
  if. -. y do.
    echo 'ASSERTION FAILED'
    exit 1
  end.
  1
)

echo 'Testing atmospheric trace gas cross-sections and plume models...'

NB. Test dimensions
assert 8 = # albedo_gases_GAS_NAMES
assert 224 = # albedo_gases_WAVELENGTHS
assert (8 224) -: $ albedo_gases_SIGMA_MATRIX

NB. Test lookups
ch4_idx =: albedo_gases_get_gas_index 'CH4'
assert 0 = ch4_idx
co2_idx =: albedo_gases_get_gas_index 'CO2'
assert 1 = co2_idx
h2o_idx =: albedo_gases_get_gas_index 'H2O'
assert 7 = h2o_idx

NB. Non-negative cross sections
assert *./ , albedo_gases_SIGMA_MATRIX >: 0.0

NB. CH4 has primary peak near 2.315 um
NB. 2.315 um corresponds to index: round((2.315 - 0.400) / 0.009375) = 204
ch4_sigma =: albedo_gases_get_gas_cross_section 'CH4'
idx_2315 =: 204
assert (idx_2315 { ch4_sigma) > 20.0
assert (idx_2315 { ch4_sigma) > (50 { ch4_sigma)

NB. Test Beer-Lambert transmittance
NB. Zero column density -> transmission = 1.0 everywhere
t_zero =: albedo_gases_gas_transmittance ('CH4' ; 0)
assert *./ t_zero = 1.0

NB. Strong column density -> transmission < 1.0 at absorption band
t_plume =: albedo_gases_gas_transmittance ('CH4' ; 5000)
assert *./ (t_plume >: 0.0) *. (t_plume <: 1.0)
assert (idx_2315 { t_plume) < 0.90

NB. Radiance perturbation is negative in absorption region
L0 =: 224 # 50.0
dL =: albedo_gases_gas_radiance_perturbation (L0 ; 'CH4' ; 2000)
assert (idx_2315 { dL) < 0.0

echo 'Atmospheric trace gas tests passed successfully.'
exit 0
