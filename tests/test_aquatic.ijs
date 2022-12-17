NB. ==============================================================================
NB. ALBEDO - Tests: data/aquatic_optics.ijs
NB. ==============================================================================

load 'data/aquatic_optics.ijs'

assert =: 3 : 0
  if. -. y do.
    echo 'ASSERTION FAILED'
    exit 1
  end.
  1
)

echo 'Testing bio-optical aquatic reflectance models...'

NB. Test dimensions
assert 4 = # albedo_aquatic_AQUATIC_NAMES
assert 224 = # albedo_aquatic_WAVELENGTHS
assert (4 224) -: $ albedo_aquatic_AQUATIC_MATRIX

NB. Clear ocean reflectance peaks in blue (near 0.440 um, index ~4)
clear_rrs =: 0 { albedo_aquatic_AQUATIC_MATRIX
blue_idx  =: 4
swir_idx  =: 100
assert (blue_idx { clear_rrs) > (swir_idx { clear_rrs)

NB. Water reflectance in SWIR is virtually zero due to immense water absorption
assert (150 { clear_rrs) < 1e_4

NB. Algal bloom has green peak (near 0.550 um, index 16) and red absorption trough at 0.675 um (index 29)
algal_rrs =: 1 { albedo_aquatic_AQUATIC_MATRIX
green_idx =: 16
red_idx   =: 29
assert (green_idx { algal_rrs) > (red_idx { algal_rrs)

NB. Turbid coastal water with high suspended sediment has much higher red backscattering
turbid_rrs =: 2 { albedo_aquatic_AQUATIC_MATRIX
assert (red_idx { turbid_rrs) > (red_idx { clear_rrs)

echo 'Bio-optical aquatic reflectance tests passed successfully.'
exit 0
