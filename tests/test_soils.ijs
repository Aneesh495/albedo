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
assert 25 = # albedo_soils_SOIL_NAMES
assert 224 = # albedo_soils_WAVELENGTHS
assert (25 224) -: $ albedo_soils_SOIL_MATRIX

NB. Test lookups
sand_idx =: albedo_soils_get_soil_index 'Sand'
assert 3 = sand_idx
peat_idx =: albedo_soils_get_soil_index 'Peat'
assert 8 = peat_idx
ox_idx   =: albedo_soils_get_soil_index 'Oxisol'
assert 10 = ox_idx

NB. Test physical bounds: all reflectance values in [0.02, 0.85]
assert *./ , (albedo_soils_SOIL_MATRIX >: 0.02) *. (albedo_soils_SOIL_MATRIX <: 0.85)

NB. Peat (organic rich) has much lower reflectance than Sand or Aridisol
peat_spec =: albedo_soils_get_soil_spectrum 'Peat'
sand_spec =: albedo_soils_get_soil_spectrum 'Sand'
arid_spec =: albedo_soils_get_soil_spectrum 'Aridisol'
assert (+/ peat_spec) < (+/ sand_spec)
assert (+/ peat_spec) < (+/ arid_spec)

NB. Iron-rich soils (Laterite, Oxisol) exhibit Fe3+ crystal field absorption near 0.88 um
NB. 0.880 um corresponds to index: round((0.880 - 0.400) / 0.009375) = 51
fe_idx =: 51
laterite_spec =: albedo_soils_get_soil_spectrum 'Laterite'
assert (fe_idx { laterite_spec) < ((fe_idx - 15) { laterite_spec)

NB. Clay-rich soil (Vertisol) exhibits strong 2.206 um Al-OH absorption
NB. 2.206 um corresponds to index: round((2.206 - 0.400) / 0.009375) = 193
clay_idx =: 193
vertisol_spec =: albedo_soils_get_soil_spectrum 'Vertisol'
assert (clay_idx { vertisol_spec) < ((clay_idx - 10) { vertisol_spec)

echo 'USDA and ECOSTRESS soil library tests passed successfully.'
exit 0
