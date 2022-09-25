NB. ==============================================================================
NB. ALBEDO - Tests: data/ecostress_soils.ijs
NB. ==============================================================================

load 'data/ecostress_soils.ijs'

assert =: 3 : 0
  if. -. y do.
    echo 'ASSERTION FAILED'
    exit 1
  end.
  1
)

echo 'Testing USDA and ECOSTRESS soil library...'

NB. Test dimensions
assert 25 = # SOIL_NAMES_albedo_soils_
assert 224 = # WAVELENGTHS_albedo_soils_
assert (25 224) -: $ SOIL_MATRIX_albedo_soils_

NB. Test lookups
sand_idx =: get_soil_index_albedo_soils_ 'Sand'
assert 3 = sand_idx
peat_idx =: get_soil_index_albedo_soils_ 'Peat'
assert 8 = peat_idx
ox_idx   =: get_soil_index_albedo_soils_ 'Oxisol'
assert 10 = ox_idx

NB. Test physical bounds: all reflectance values in [0.02, 0.85]
assert *./ , (SOIL_MATRIX_albedo_soils_ >: 0.02) *. (SOIL_MATRIX_albedo_soils_ <: 0.85)

NB. Peat (organic rich) has much lower reflectance than Sand or Aridisol
peat_spec =: get_soil_spectrum_albedo_soils_ 'Peat'
sand_spec =: get_soil_spectrum_albedo_soils_ 'Sand'
arid_spec =: get_soil_spectrum_albedo_soils_ 'Aridisol'
assert (+/ peat_spec) < (+/ sand_spec)
assert (+/ peat_spec) < (+/ arid_spec)

NB. Iron-rich soils (Laterite, Oxisol) exhibit Fe3+ crystal field absorption near 0.88 um
NB. 0.880 um corresponds to index: round((0.880 - 0.400) / 0.009375) = 51
fe_idx =: 51
laterite_spec =: get_soil_spectrum_albedo_soils_ 'Laterite'
assert (fe_idx { laterite_spec) < ((fe_idx - 15) { laterite_spec)

NB. Clay-rich soil (Vertisol) exhibits strong 2.206 um Al-OH absorption
NB. 2.206 um corresponds to index: round((2.206 - 0.400) / 0.009375) = 193
clay_idx =: 193
vertisol_spec =: get_soil_spectrum_albedo_soils_ 'Vertisol'
assert (clay_idx { vertisol_spec) < ((clay_idx - 10) { vertisol_spec)

echo 'USDA and ECOSTRESS soil library tests passed successfully.'
exit 0
