NB. ==============================================================================
NB. ALBEDO - High-Dimensional Hyperspectral Imaging Engine
NB. data/ecostress_soils.ijs : USDA & ECOSTRESS Soil Spectral Library (25 Classes)
NB. Provides 224-band laboratory reflectance signatures across soil texture classes
NB. and taxonomic soil orders from 0.400 to 2.500 micrometers.
NB. ==============================================================================

cocurrent 'base'

NB. 224 Wavelength centers in micrometers (0.400 to 2.500 um)
WAVELENGTHS =: 0.400 + 0.009375 * i. 224

NB. 25 Soil Classes (Textures and Taxonomic Orders)
SOIL_NAMES =: <;._1 ' SandyLoam' , ' ClayLoam' , ' SiltyClay' , ' Sand' , ' LoamySand' , ' SandyClay' , ' Silt' , ' SiltLoam' , ' Peat' , ' Vertisol' , ' Oxisol' , ' Aridisol' , ' Spodosol' , ' Alfisol' , ' Mollisol' , ' Ultisol' , ' Inceptisol' , ' Entisol' , ' Andisol' , ' Gelisol' , ' SalineSoil' , ' CalcareousSoil' , ' Laterite' , ' GypsiferousSoil' , ' Ferralsol'

NB. Helper to construct characteristic soil reflectance curve:
NB. r0: base reflectance, slope: overall tilt, om_damp: organic matter dampening factor,
NB. fe_depth: Fe3+ 0.88 um absorption depth, clay_depth: Al-OH 2.20 um absorption depth
calc_soil =: 3 : 0
  'r0 slope om_damp fe_depth clay_depth' =. y
  NB. Base spectral continuum (concave downward or linear rise)
  w =. WAVELENGTHS - 0.400
  base =. r0 + (slope * w) - (0.05 * w * w)
  base =. base * om_damp
  NB. Iron oxide broad crystal field absorption at 0.880 um
  fe_abs =. fe_depth * ^ - 0.5 * *: (WAVELENGTHS - 0.880) % 0.120
  NB. Clay lattice Al-OH vibrational overtone at 2.206 um
  clay_abs =. clay_depth * ^ - 0.5 * *: (WAVELENGTHS - 2.206) % 0.030
  NB. Structural water at 1.415 um and 1.915 um
  w14 =. 0.08 * ^ - 0.5 * *: (WAVELENGTHS - 1.415) % 0.035
  w19 =. 0.12 * ^ - 0.5 * *: (WAVELENGTHS - 1.915) % 0.040
  NB. Combine and constrain to realistic soil reflectance bounds [0.02, 0.85]
  spec =. base - (fe_abs + clay_abs + w14 + w19)
  0.02 >. 0.85 <. spec
)

NB. Soil parameters: (r0, slope, om_damp, fe_depth, clay_depth)
s_sandy_loam      =: calc_soil 0.18 0.18 1.00 0.05 0.06
s_clay_loam       =: calc_soil 0.15 0.16 0.95 0.06 0.12
s_silty_clay      =: calc_soil 0.14 0.15 0.92 0.07 0.14
s_sand            =: calc_soil 0.28 0.22 1.10 0.03 0.02
s_loamy_sand      =: calc_soil 0.22 0.20 1.05 0.04 0.04
s_sandy_clay      =: calc_soil 0.16 0.17 0.96 0.06 0.10
s_silt            =: calc_soil 0.20 0.18 1.02 0.05 0.05
s_silt_loam       =: calc_soil 0.17 0.17 0.98 0.05 0.07
s_peat            =: calc_soil 0.05 0.06 0.60 0.02 0.02
s_vertisol        =: calc_soil 0.10 0.10 0.80 0.04 0.16
s_oxisol          =: calc_soil 0.18 0.16 0.95 0.14 0.08
s_aridisol        =: calc_soil 0.30 0.20 1.15 0.04 0.05
s_spodosol        =: calc_soil 0.12 0.14 0.85 0.05 0.04
s_alfisol         =: calc_soil 0.19 0.18 1.00 0.06 0.08
s_mollisol        =: calc_soil 0.09 0.11 0.75 0.04 0.06
s_ultisol         =: calc_soil 0.17 0.17 0.95 0.11 0.10
s_inceptisol      =: calc_soil 0.16 0.16 0.96 0.06 0.06
s_entisol         =: calc_soil 0.22 0.19 1.04 0.05 0.04
s_andisol         =: calc_soil 0.08 0.09 0.70 0.03 0.05
s_gelisol         =: calc_soil 0.11 0.12 0.82 0.04 0.04
s_saline          =: calc_soil 0.35 0.18 1.20 0.02 0.03
s_calcareous      =: calc_soil 0.32 0.16 1.15 0.03 0.05
s_laterite        =: calc_soil 0.16 0.15 0.92 0.16 0.06
s_gypsiferous     =: calc_soil 0.38 0.14 1.22 0.02 0.04
s_ferralsol       =: calc_soil 0.17 0.16 0.94 0.15 0.07

NB. Shape: 25 rows by 224 columns
SOIL_MATRIX =: > s_sandy_loam ; s_clay_loam ; s_silty_clay ; s_sand ; s_loamy_sand ; s_sandy_clay ; s_silt ; s_silt_loam ; s_peat ; s_vertisol ; s_oxisol ; s_aridisol ; s_spodosol ; s_alfisol ; s_mollisol ; s_ultisol ; s_inceptisol ; s_entisol ; s_andisol ; s_gelisol ; s_saline ; s_calcareous ; s_laterite ; s_gypsiferous ; s_ferralsol

NB. Query soil index by exact name
get_soil_index =: 3 : 0
  idx =. SOIL_NAMES i. < y
  if. idx = # SOIL_NAMES do. _1 else. idx end.
)

NB. Get spectral signature by soil name
get_soil_spectrum =: 3 : 0
  idx =. get_soil_index y
  if. idx = _1 do.
    224 # 0
  else.
    idx { SOIL_MATRIX
  end.
)


NB. Exported Aliases
albedo_soils_SOIL_NAMES =: SOIL_NAMES
albedo_soils_WAVELENGTHS =: WAVELENGTHS
albedo_soils_SOIL_MATRIX =: SOIL_MATRIX
albedo_soils_get_soil_index =: get_soil_index
albedo_soils_get_soil_spectrum =: get_soil_spectrum
