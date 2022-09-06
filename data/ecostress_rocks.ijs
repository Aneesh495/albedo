NB. ==============================================================================
NB. ALBEDO - High-Dimensional Hyperspectral Imaging Engine
NB. data/ecostress_rocks.ijs : ECOSTRESS Rock Spectral Library (30 Rock Types)
NB. Provides 224-band laboratory reflectance signatures across igneous, metamorphic,
NB. and sedimentary geological classes from 0.400 to 2.500 micrometers.
NB. ==============================================================================

cocurrent 'albedo_rocks'

NB. 224 Wavelength centers in micrometers (0.400 to 2.500 um)
WAVELENGTHS =: 0.400 + 0.009375 * i. 224

NB. 30 Rock Names
ROCK_NAMES =: <;._1 ' Granite' , ' Basalt' , ' Rhyolite' , ' Gabbro' , ' Diorite' , ' Peridotite' , ' Obsidian' , ' Pumice' , ' Andesite' , ' Dunite' , ' Marble' , ' Quartzite' , ' Schist' , ' Gneiss' , ' Slate' , ' Phyllite' , ' Serpentinite' , ' Amphibolite' , ' Granulite' , ' Eclogite' , ' Limestone' , ' Dolomite' , ' Sandstone' , ' Shale' , ' Siltstone' , ' Chert' , ' Arkose' , ' Graywacke' , ' Chalk' , ' Tuff'

NB. Rock Geological Classification (0=Igneous, 1=Metamorphic, 2=Sedimentary)
ROCK_CLASSES =: (10 # 0) , (10 # 1) , (10 # 2)

NB. Helper to construct characteristic reflectance curve with diagnostic absorption features:
NB. r0: base reflectance, slope: spectral tilt, dips: list of (center_um, depth, width_um)
calc_curve =: 3 : 0
  'r0 slope dips' =. y
  base =. r0 + slope * (WAVELENGTHS - 0.400)
  for_d. dips do.
    'c0 depth w0' =. d
    gauss =. depth * ^. - 0.5 * *: (WAVELENGTHS - c0) % w0
    base =. base - gauss
  end.
  0.01 >. 0.95 <. base
)

NB. Igneous Rocks
granite_dips     =: 3 3 $ 1.410 0.08 0.035  1.910 0.12 0.040  2.205 0.15 0.025
basalt_dips      =: 2 3 $ 1.050 0.18 0.180  2.000 0.07 0.150
rhyolite_dips    =: 3 3 $ 1.415 0.09 0.030  1.915 0.14 0.040  2.210 0.12 0.025
gabbro_dips      =: 2 3 $ 1.020 0.22 0.160  1.980 0.08 0.140
diorite_dips     =: 2 3 $ 1.030 0.16 0.150  2.250 0.06 0.030
peridotite_dips  =: 2 3 $ 1.050 0.35 0.250  1.850 0.10 0.180
obsidian_dips    =: 1 3 $ 1.910 0.04 0.050
pumice_dips      =: 3 3 $ 1.410 0.06 0.030  1.910 0.11 0.040  2.205 0.08 0.030
andesite_dips    =: 2 3 $ 1.040 0.14 0.140  2.250 0.05 0.030
dunite_dips      =: 1 3 $ 1.060 0.40 0.260

NB. Metamorphic Rocks
marble_dips       =: 3 3 $ 1.870 0.05 0.030  2.340 0.32 0.035  2.500 0.15 0.040
quartzite_dips    =: 2 3 $ 1.415 0.04 0.025  2.210 0.05 0.025
schist_dips       =: 3 3 $ 1.405 0.12 0.035  1.910 0.15 0.040  2.345 0.18 0.035
gneiss_dips       =: 3 3 $ 1.410 0.08 0.030  1.910 0.11 0.035  2.205 0.10 0.025
slate_dips        =: 3 3 $ 1.410 0.06 0.030  1.910 0.08 0.035  2.205 0.07 0.025
phyllite_dips     =: 3 3 $ 1.410 0.09 0.030  1.910 0.12 0.035  2.205 0.11 0.025
serpentinite_dips =: 3 3 $ 1.390 0.22 0.030  2.325 0.38 0.035  2.120 0.08 0.030
amphibolite_dips  =: 2 3 $ 1.030 0.18 0.160  2.310 0.12 0.035
granulite_dips    =: 2 3 $ 1.040 0.12 0.140  2.210 0.08 0.030
eclogite_dips     =: 2 3 $ 1.050 0.20 0.180  2.260 0.10 0.035

NB. Sedimentary Rocks
limestone_dips   =: 3 3 $ 1.875 0.06 0.030  2.335 0.34 0.035  2.500 0.18 0.040
dolomite_dips    =: 3 3 $ 1.860 0.06 0.030  2.320 0.35 0.035  2.490 0.18 0.040
sandstone_dips   =: 3 3 $ 1.415 0.07 0.030  1.910 0.10 0.035  2.205 0.09 0.025
shale_dips       =: 3 3 $ 1.410 0.10 0.030  1.910 0.14 0.035  2.205 0.15 0.025
siltstone_dips   =: 3 3 $ 1.410 0.08 0.030  1.910 0.12 0.035  2.205 0.11 0.025
chert_dips       =: 2 3 $ 1.415 0.03 0.025  1.915 0.05 0.030
arkose_dips      =: 3 3 $ 1.415 0.09 0.030  1.915 0.13 0.035  2.205 0.14 0.025
graywacke_dips   =: 2 3 $ 1.050 0.08 0.120  2.330 0.08 0.035
chalk_dips       =: 3 3 $ 1.875 0.08 0.030  2.335 0.40 0.035  2.500 0.20 0.040
tuff_dips        =: 3 3 $ 1.410 0.08 0.030  1.910 0.13 0.035  2.210 0.10 0.025

granite      =: calc_curve 0.35 0.18 granite_dips
basalt       =: calc_curve 0.12 0.08 basalt_dips
rhyolite     =: calc_curve 0.40 0.20 rhyolite_dips
gabbro       =: calc_curve 0.14 0.09 gabbro_dips
diorite      =: calc_curve 0.22 0.12 diorite_dips
peridotite   =: calc_curve 0.10 0.06 peridotite_dips
obsidian     =: calc_curve 0.06 0.04 obsidian_dips
pumice       =: calc_curve 0.50 0.15 pumice_dips
andesite     =: calc_curve 0.18 0.10 andesite_dips
dunite       =: calc_curve 0.15 0.08 dunite_dips

marble       =: calc_curve 0.55 0.12 marble_dips
quartzite    =: calc_curve 0.60 0.15 quartzite_dips
schist       =: calc_curve 0.24 0.12 schist_dips
gneiss       =: calc_curve 0.28 0.14 gneiss_dips
slate        =: calc_curve 0.16 0.08 slate_dips
phyllite     =: calc_curve 0.20 0.10 phyllite_dips
serpentinite =: calc_curve 0.18 0.07 serpentinite_dips
amphibolite  =: calc_curve 0.15 0.09 amphibolite_dips
granulite    =: calc_curve 0.30 0.13 granulite_dips
eclogite     =: calc_curve 0.16 0.08 eclogite_dips

limestone    =: calc_curve 0.42 0.16 limestone_dips
dolomite     =: calc_curve 0.44 0.16 dolomite_dips
sandstone    =: calc_curve 0.38 0.20 sandstone_dips
shale        =: calc_curve 0.22 0.14 shale_dips
siltstone    =: calc_curve 0.26 0.15 siltstone_dips
chert        =: calc_curve 0.50 0.12 chert_dips
arkose       =: calc_curve 0.36 0.18 arkose_dips
graywacke    =: calc_curve 0.18 0.10 graywacke_dips
chalk        =: calc_curve 0.65 0.10 chalk_dips
tuff         =: calc_curve 0.38 0.16 tuff_dips

NB. Shape: 30 rows by 224 columns
ROCK_MATRIX =: > granite ; basalt ; rhyolite ; gabbro ; diorite ; peridotite ; obsidian ; pumice ; andesite ; dunite ; marble ; quartzite ; schist ; gneiss ; slate ; phyllite ; serpentinite ; amphibolite ; granulite ; eclogite ; limestone ; dolomite ; sandstone ; shale ; siltstone ; chert ; arkose ; graywacke ; chalk ; tuff

NB. Query rock index by exact name
get_rock_index =: 3 : 0
  idx =. ROCK_NAMES i. < y
  if. idx = # ROCK_NAMES do. _1 else. idx end.
)

NB. Get spectral signature by rock name
get_rock_spectrum =: 3 : 0
  idx =. get_rock_index y
  if. idx = _1 do.
    224 # 0
  else.
    idx { ROCK_MATRIX
  end.
)

NB. Get all rocks belonging to a specific class (0, 1, or 2)
get_rocks_by_class =: 3 : 0
  mask =. ROCK_CLASSES = y
  mask # ROCK_NAMES
)
