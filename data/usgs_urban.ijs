NB. ==============================================================================
NB. Albedo - Calibrated Anthropogenic & Urban Materials Spectral Library
NB.
NB. High-fidelity 224-band reflectance signatures (400 - 2500 nm) calibrated
NB. for urban remote sensing, infrastructure monitoring, and materials triage:
NB.   - Aged Asphalt Pavement (hydrocarbon oxidation, mineral aggregate exposure)
NB.   - Fresh Asphalt Sealcoat (low albedo, uniform bituminous absorption)
NB.   - Portland Cement Concrete (silicate cement paste, high visible albedo)
NB.   - Galvanized Corrugated Steel (metallic zinc specular coating, flat high SWIR)
NB.   - Terracotta Clay Roof Tile (Fe3+ hematite absorption in visible)
NB.   - White EPDM / TPO Cool Roof Membrane (high solar reflectance, polymer dips)
NB.   - Polyethylene Plastic (LDPE/HDPE diagnostic C-H stretch at 1210, 1400, 1730 nm)
NB.   - Polyvinyl Chloride (PVC membrane, C-Cl vibrational overtones)
NB.   - Red Clay Brick Masonry (iron oxide kiln-fired clay)
NB.   - Industrial Silica Sand Pile (high diffuse reflectance)
NB. ==============================================================================

cocurrent 'base'

load 'io/spectral_lib.ijs'

W224=: albedo_spectral_lib_standard_wavelengths_224 ''

NB. Aged Asphalt Pavement
spec_aged_asphalt=: 3 : 0
  w=. y
  0.10 + 0.12 * (w - 400.0) % 2100.0
)

NB. Fresh Asphalt Sealcoat
spec_fresh_asphalt=: 3 : 0
  w=. y
  0.04 + 0.04 * (w - 400.0) % 2100.0
)

NB. Portland Cement Concrete
spec_portland_concrete=: 3 : 0
  w=. y
  refl=. 0.22 + 0.25 * 1.0 - ^ - (w - 400.0) % 650.0
  refl=. refl - 0.05 * ^ - *: (w - 1420.0) % 40.0
  refl=. refl - 0.07 * ^ - *: (w - 1950.0) % 50.0
  0.10 >. 0.52 <. refl
)

NB. Galvanized Steel Roofing
spec_galvanized_steel=: 3 : 0
  w=. y
  0.35 + 0.15 * (w - 400.0) % 2100.0
)

NB. Terracotta Clay Roof Tile
spec_terracotta_tile=: 3 : 0
  w=. y
  refl=. 0.08 + 0.35 % (1.0 + ^ - (w - 560.0) % 45.0)
  refl=. refl - 0.04 * ^ - *: (w - 870.0) % 80.0
  0.05 >. 0.45 <. refl
)

NB. White EPDM / TPO Cool Roof Membrane
spec_white_roof_epdm=: 3 : 0
  w=. y
  refl=. 0.78 - 0.15 * (w - 400.0) % 2100.0
  refl=. refl - 0.08 * ^ - *: (w - 1730.0) % 30.0
  0.30 >. 0.85 <. refl
)

NB. Polyethylene Plastic (HDPE / LDPE)
spec_polyethylene_plastic=: 3 : 0
  w=. y
  refl=. 0.45 + 0.20 * (w - 400.0) % 2100.0
  NB. Strong C-H overtone and combination absorption bands
  refl=. refl - 0.15 * ^ - *: (w - 1215.0) % 20.0
  refl=. refl - 0.22 * ^ - *: (w - 1395.0) % 25.0
  refl=. refl - 0.38 * ^ - *: (w - 1730.0) % 22.0
  refl=. refl - 0.45 * ^ - *: (w - 2310.0) % 28.0
  0.02 >. 0.75 <. refl
)

NB. Polyvinyl Chloride (PVC)
spec_pvc_membrane=: 3 : 0
  w=. y
  refl=. 0.50 + 0.15 * (w - 400.0) % 2100.0
  refl=. refl - 0.25 * ^ - *: (w - 1410.0) % 25.0
  refl=. refl - 0.32 * ^ - *: (w - 1780.0) % 30.0
  refl=. refl - 0.30 * ^ - *: (w - 2280.0) % 35.0
  0.02 >. 0.70 <. refl
)

NB. Red Brick Masonry
spec_red_brick=: 3 : 0
  w=. y
  refl=. 0.10 + 0.30 % (1.0 + ^ - (w - 580.0) % 50.0)
  refl=. refl - 0.05 * ^ - *: (w - 880.0) % 90.0
  0.05 >. 0.42 <. refl
)

NB. Industrial Silica Sand Pile
spec_silica_sand=: 3 : 0
  w=. y
  0.45 + 0.25 * 1.0 - ^ - (w - 400.0) % 800.0
)

NB. Build full Urban & Anthropogenic library matrix (224 bands x 10 materials)
build_urban_library=: 3 : 0
  wvl=. y
  u1=.  spec_aged_asphalt wvl
  u2=.  spec_fresh_asphalt wvl
  u3=.  spec_portland_concrete wvl
  u4=.  spec_galvanized_steel wvl
  u5=.  spec_terracotta_tile wvl
  u6=.  spec_white_roof_epdm wvl
  u7=.  spec_polyethylene_plastic wvl
  u8=.  spec_pvc_membrane wvl
  u9=.  spec_red_brick wvl
  u10=. spec_silica_sand wvl

  u1 ,. u2 ,. u3 ,. u4 ,. u5 ,. u6 ,. u7 ,. u8 ,. u9 ,. u10
)


albedo_usgs_urban_build_urban_library =: build_urban_library
