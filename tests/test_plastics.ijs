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
assert 20 = # albedo_plastics_PLASTIC_NAMES
assert 224 = # albedo_plastics_WAVELENGTHS
assert (20 224) -: $ albedo_plastics_PLASTIC_MATRIX

NB. Test lookups
pet_idx  =: albedo_plastics_get_plastic_index 'PET'
assert 0 = pet_idx
hdpe_idx =: albedo_plastics_get_plastic_index 'HDPE'
assert 1 = hdpe_idx
eps_idx  =: albedo_plastics_get_plastic_index 'EPSFoam'
assert 19 = eps_idx

NB. Test physical bounds: all reflectance values in [0.01, 0.98]
assert *./ , (albedo_plastics_PLASTIC_MATRIX >: 0.01) *. (albedo_plastics_PLASTIC_MATRIX <: 0.98)

NB. Polyethylene (HDPE) has intense C-H overtone at 1.73 um (index 142)
hdpe_spec =: albedo_plastics_get_plastic_spectrum 'HDPE'
idx_173 =: 142
assert (idx_173 { hdpe_spec) < ((idx_173 - 10) { hdpe_spec)
assert (idx_173 { hdpe_spec) < ((idx_173 + 10) { hdpe_spec)

NB. Test plastic hydrocarbon index on HDPE vs PTFE (Teflon has no C-H bonds)
hdpe_fpi =: albedo_plastics_plastic_hydrocarbon_index hdpe_spec
ptfe_spec =: albedo_plastics_get_plastic_spectrum 'PTFE'
ptfe_fpi =: albedo_plastics_plastic_hydrocarbon_index ptfe_spec
assert hdpe_fpi > 0.15
assert ptfe_fpi < 0.05

echo 'Synthetic polymers and plastics library tests passed successfully.'
exit 0
