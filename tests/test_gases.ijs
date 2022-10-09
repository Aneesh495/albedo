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
assert 8 = # GAS_NAMES_albedo_gases_
assert 224 = # WAVELENGTHS_albedo_gases_
assert (8 224) -: $ SIGMA_MATRIX_albedo_gases_

NB. Test lookups
ch4_idx =: get_gas_index_albedo_gases_ 'CH4'
assert 0 = ch4_idx
co2_idx =: get_gas_index_albedo_gases_ 'CO2'
assert 1 = co2_idx
h2o_idx =: get_gas_index_albedo_gases_ 'H2O'
assert 7 = h2o_idx

NB. Non-negative cross sections
assert *./ , SIGMA_MATRIX_albedo_gases_ >: 0.0

NB. CH4 has primary peak near 2.315 um
NB. 2.315 um corresponds to index: round((2.315 - 0.400) / 0.009375) = 204
ch4_sigma =: get_gas_cross_section_albedo_gases_ 'CH4'
idx_2315 =: 204
assert (idx_2315 { ch4_sigma) > 20.0
assert (idx_2315 { ch4_sigma) > (50 { ch4_sigma)

NB. Test Beer-Lambert transmittance
NB. Zero column density -> transmission = 1.0 everywhere
t_zero =: gas_transmittance_albedo_gases_ ('CH4' ; 0)
assert *./ t_zero = 1.0

NB. Strong column density -> transmission < 1.0 at absorption band
t_plume =: gas_transmittance_albedo_gases_ ('CH4' ; 5000)
assert *./ (t_plume >: 0.0) *. (t_plume <: 1.0)
assert (idx_2315 { t_plume) < 0.90

NB. Radiance perturbation is negative in absorption region
L0 =: 224 # 50.0
dL =: gas_radiance_perturbation_albedo_gases_ (L0 ; 'CH4' ; 2000)
assert (idx_2315 { dL) < 0.0

echo 'Atmospheric trace gas tests passed successfully.'
exit 0
