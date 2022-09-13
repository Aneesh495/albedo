NB. ==============================================================================
NB. ALBEDO - Tests: data/ecostress_rocks.ijs
NB. ==============================================================================

load 'data/ecostress_rocks.ijs'

assert =: 3 : 0
  if. -. y do.
    echo 'ASSERTION FAILED'
    exit 1
  end.
  1
)

approx_eq =: 4 : 0
  tol =. 1e_4
  *./ , tol > | x - y
)

echo 'Testing ECOSTRESS rock library...'

NB. Test dimensions
assert 30 = # ROCK_NAMES_albedo_rocks_
assert 224 = # WAVELENGTHS_albedo_rocks_
assert (30 224) -: $ ROCK_MATRIX_albedo_rocks_

NB. Test geological classes: 10 igneous, 10 metamorphic, 10 sedimentary
assert 10 = # get_rocks_by_class_albedo_rocks_ 0
assert 10 = # get_rocks_by_class_albedo_rocks_ 1
assert 10 = # get_rocks_by_class_albedo_rocks_ 2

NB. Test rock lookup
g_idx =: get_rock_index_albedo_rocks_ 'Granite'
assert 0 = g_idx
b_idx =: get_rock_index_albedo_rocks_ 'Basalt'
assert 1 = b_idx
m_idx =: get_rock_index_albedo_rocks_ 'Marble'
assert 10 = m_idx
l_idx =: get_rock_index_albedo_rocks_ 'Limestone'
assert 20 = l_idx

NB. Test non-existent lookup returns _1
bad_idx =: get_rock_index_albedo_rocks_ 'NonexistentStone'
assert _1 = bad_idx

NB. Test physical bounds: all reflectance values in [0.01, 1.0]
assert *./ , (ROCK_MATRIX_albedo_rocks_ >: 0.01) *. (ROCK_MATRIX_albedo_rocks_ <: 1.0)

NB. Test carbonate absorption feature at ~2.335 um for Limestone and Dolomite
NB. 2.335 um corresponds to index: round((2.335 - 0.400) / 0.009375) = round(1.935 / 0.009375) = 206
carb_idx =: 206
limestone_spec =: get_rock_spectrum_albedo_rocks_ 'Limestone'
dolomite_spec  =: get_rock_spectrum_albedo_rocks_ 'Dolomite'
assert (carb_idx { limestone_spec) < ((carb_idx - 15) { limestone_spec)
assert (carb_idx { dolomite_spec) < ((carb_idx - 15) { dolomite_spec)

NB. Basalt is mafic/dark, Chalk is highly reflective
basalt_spec =: get_rock_spectrum_albedo_rocks_ 'Basalt'
chalk_spec  =: get_rock_spectrum_albedo_rocks_ 'Chalk'
assert (+/ basalt_spec) < (+/ chalk_spec)

echo 'ECOSTRESS rock library tests passed successfully.'
exit 0
