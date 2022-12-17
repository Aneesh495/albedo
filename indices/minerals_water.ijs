NB. ==============================================================================
NB. Albedo - Mineralogical, Hydrological & Burn Spectral Indices
NB.
NB. Computes diagnostic remote sensing indices for water, moisture, hydrothermal
NB. alteration minerals (alunite, kaolinite, iron oxides), and wildfire burn scars:
NB.   - NDWI:  Normalized Difference Water Index (McFeeters, 1996)
NB.   - MNDWI: Modified Normalized Difference Water Index (Xu, 2006)
NB.   - NDMI:  Normalized Difference Moisture Index (Gao, 1996)
NB.   - NBR:   Normalized Burn Ratio (Key & Benson, 2006)
NB.   - Clay Mineral Index (SWIR-1 / SWIR-2)
NB.   - Ferrous Iron Index (SWIR / NIR)
NB.   - Ferric Iron Index (Red / Blue)
NB. ==============================================================================

cocurrent 'base'

load 'io/cube.ijs'

get_band=: 4 : 0
  target_nm=. x
  cube=. y
  'data wvl meta'=. cube
  idx=. target_nm albedo_cube_find_nearest_band wvl
  idx { data
)

norm_diff=: 4 : 0
  b1=. x
  b2=. y
  diff=. b2 - b1
  sum=. b2 + b1
  denom=. sum + 1e_6 * sum = 0.0
  _1.0 >. 1.0 <. (diff % denom)
)

NB. Normalized Difference Water Index (NDWI)
NB. Green: 560 nm, NIR: 860 nm
calc_ndwi=: 3 : 0
  cube=. y
  green=. 560.0 get_band cube
  nir=.   860.0 get_band cube
  nir norm_diff green
)

NB. Modified NDWI (MNDWI)
NB. Green: 560 nm, SWIR-1: 1610 nm
calc_mndwi=: 3 : 0
  cube=. y
  green=. 560.0 get_band cube
  swir1=. 1610.0 get_band cube
  swir1 norm_diff green
)

NB. Normalized Difference Moisture Index (NDMI)
NB. NIR: 860 nm, SWIR-1: 1610 nm
calc_ndmi=: 3 : 0
  cube=. y
  nir=.   860.0 get_band cube
  swir1=. 1610.0 get_band cube
  swir1 norm_diff nir
)

NB. Normalized Burn Ratio (NBR)
NB. NIR: 860 nm, SWIR-2: 2200 nm
calc_nbr=: 3 : 0
  cube=. y
  nir=.   860.0 get_band cube
  swir2=. 2200.0 get_band cube
  swir2 norm_diff nir
)

NB. Clay Minerals Index (Al-OH absorption): SWIR-1 (1600 nm) / SWIR-2 (2200 nm)
calc_clay_index=: 3 : 0
  cube=. y
  swir1=. 1600.0 get_band cube
  swir2=. 2200.0 get_band cube
  swir1 % (swir2 >. 1e_6)
)

NB. Ferrous Iron (Fe2+) Index: SWIR-1 (1600 nm) / NIR (860 nm)
calc_ferrous_iron=: 3 : 0
  cube=. y
  swir1=. 1600.0 get_band cube
  nir=.    860.0 get_band cube
  swir1 % (nir >. 1e_6)
)

NB. Ferric Iron (Fe3+) Index: Red (660 nm) / Blue (480 nm)
calc_ferric_iron=: 3 : 0
  cube=. y
  red=.  660.0 get_band cube
  blue=. 480.0 get_band cube
  red % (blue >. 1e_6)
)


albedo_indices_minerals_calc_clay_index =: calc_clay_index
albedo_indices_minerals_calc_ndwi =: calc_ndwi
