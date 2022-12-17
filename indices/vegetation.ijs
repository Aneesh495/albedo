NB. ==============================================================================
NB. Albedo - Biochemical Vegetation Spectral Indices
NB.
NB. Computes standard biophysical and biochemical indices from hyperspectral cubes:
NB.   - NDVI: Normalized Difference Vegetation Index (Rouse et al., 1974)
NB.   - EVI:  Enhanced Vegetation Index (Huete et al., 2002)
NB.   - NDRE: Normalized Difference Red Edge (Gitelson & Merzlyak, 1994)
NB.   - PRI:  Photochemical Reflectance Index (Gamon et al., 1992)
NB.   - MCARI: Modified Chlorophyll Absorption in Reflectance Index
NB. ==============================================================================

cocurrent 'base'

load 'io/cube.ijs'

NB. Helper to extract nearest band to target wavelength
get_band=: 4 : 0
  target_nm=. x
  cube=. y
  'data wvl meta'=. cube
  idx=. target_nm albedo_cube_find_nearest_band wvl
  idx { data
)

NB. Normalized Difference helper: (b2 - b1) / (b2 + b1)
norm_diff=: 4 : 0
  b1=. x
  b2=. y
  diff=. b2 - b1
  sum=. b2 + b1
  denom=. sum + 1e_6 * sum = 0.0
  _1.0 >. 1.0 <. (diff % denom)
)

NB. Normalized Difference Vegetation Index (NDVI)
NB. Red: 660 nm, NIR: 840 nm
calc_ndvi=: 3 : 0
  cube=. y
  red=. 660.0 get_band cube
  nir=. 840.0 get_band cube
  red norm_diff nir
)

NB. Normalized Difference Red Edge (NDRE)
NB. Red Edge: 705 nm, NIR: 780 nm
calc_ndre=: 3 : 0
  cube=. y
  re=. 705.0 get_band cube
  nir=. 780.0 get_band cube
  re norm_diff nir
)

NB. Photochemical Reflectance Index (PRI)
NB. Xanthophyll pigment cycle proxy: 531 nm vs 570 nm
calc_pri=: 3 : 0
  cube=. y
  b531=. 531.0 get_band cube
  b570=. 570.0 get_band cube
  b570 norm_diff b531
)

NB. Enhanced Vegetation Index (EVI)
NB. 2.5 * (NIR - Red) / (NIR + 6*Red - 7.5*Blue + 1)
calc_evi=: 3 : 0
  cube=. y
  blue=. 470.0 get_band cube
  red=.  660.0 get_band cube
  nir=.  840.0 get_band cube
  num=. nir - red
  denom=. (nir + (6.0 * red) - (7.5 * blue) + 1.0) >. 1e_6
  _1.0 >. 1.0 <. (2.5 * num % denom)
)

NB. Modified Chlorophyll Absorption in Reflectance Index (MCARI)
NB. ((R700 - R670) - 0.2*(R700 - R550)) * (R700 / R670)
calc_mcari=: 3 : 0
  cube=. y
  r550=. 550.0 get_band cube
  r670=. 670.0 get_band cube
  r700=. 700.0 get_band cube
  part1=. (r700 - r670) - 0.2 * (r700 - r550)
  part2=. r700 % (r670 >. 1e_6)
  part1 * part2
)


albedo_indices_veg_norm_diff =: norm_diff
albedo_indices_veg_calc_ndvi =: calc_ndvi
albedo_indices_veg_calc_ndre =: calc_ndre
albedo_indices_veg_calc_evi =: calc_evi
