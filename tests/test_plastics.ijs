NB. ==============================================================================
NB. ALBEDO - Tests: data/polymers_plastics.ijs
NB. ==============================================================================

load 'data/polymers_plastics.ijs'

assert =: 3 : 0
  if. -. y do.
    echo 'ASSERTION FAILED'
    exit 1
  end.
  1
)

echo 'Testing synthetic polymers and plastics library...'

NB. Test dimensions
assert 20 = # PLASTIC_NAMES_albedo_plastics_
assert 224 = # WAVELENGTHS_albedo_plastics_
assert (20 224) -: $ PLASTIC_MATRIX_albedo_plastics_

NB. Test lookups
pet_idx  =: get_plastic_index_albedo_plastics_ 'PET'
assert 0 = pet_idx
hdpe_idx =: get_plastic_index_albedo_plastics_ 'HDPE'
assert 1 = hdpe_idx
eps_idx  =: get_plastic_index_albedo_plastics_ 'EPSFoam'
assert 19 = eps_idx

NB. Test physical bounds: all reflectance values in [0.01, 0.98]
assert *./ , (PLASTIC_MATRIX_albedo_plastics_ >: 0.01) *. (PLASTIC_MATRIX_albedo_plastics_ <: 0.98)

NB. Polyethylene (HDPE) has intense C-H overtone at 1.73 um (index 142)
hdpe_spec =: get_plastic_spectrum_albedo_plastics_ 'HDPE'
idx_173 =: 142
assert (idx_173 { hdpe_spec) < ((idx_173 - 10) { hdpe_spec)
assert (idx_173 { hdpe_spec) < ((idx_173 + 10) { hdpe_spec)

NB. Test plastic hydrocarbon index on HDPE vs PTFE (Teflon has no C-H bonds)
hdpe_fpi =: plastic_hydrocarbon_index_albedo_plastics_ hdpe_spec
ptfe_spec =: get_plastic_spectrum_albedo_plastics_ 'PTFE'
ptfe_fpi =: plastic_hydrocarbon_index_albedo_plastics_ ptfe_spec
assert hdpe_fpi > 0.15
assert ptfe_fpi < 0.05

echo 'Synthetic polymers and plastics library tests passed successfully.'
exit 0
