NB. ==============================================================================
NB. Albedo - Calibrated Botanical & Agricultural Spectral Library
NB.
NB. High-fidelity 224-band reflectance signatures (400 - 2500 nm) calibrated
NB. for forestry, precision agriculture, and aquatic vegetation:
NB.   - Deciduous Broadleaf (chlorophyll a/b absorption, steep red edge, NIR plateau)
NB.   - Conifer Pine Needle (dense mesophyll, reduced NIR scattering)
NB.   - Dry Senescent Corn / Maize (cellulose, lignin, dry biomass)
NB.   - Soybean Crop Canopy (foliar nitrogen, chlorophyll reflectance peak)
NB.   - Vineyard Grapevine (anthocyanin and carotenoid pigments)
NB.   - Wetland Sedge / Reed (high internal moisture, standing water background)
NB.   - Sphagnum Peat Moss (sponge-like cellular water retention)
NB.   - Freshwater Microalgae (chlorophyll-a at 675 nm and phycocyanin at 620 nm)
NB.   - Biological Soil Crust / Lichen (fungal-algal symbiosis, dry pigments)
NB.   - Turf Grass / Golf Fairway (uniform managed turf canopy)
NB. ==============================================================================

cocurrent 'base'

load 'io/spectral_lib.ijs'

W224=: albedo_spectral_lib_standard_wavelengths_224 ''

NB. Deciduous Broadleaf (Oak / Maple)
spec_deciduous_broadleaf=: 3 : 0
  w=. y
  refl=. 0.04 + 0.08 * ^ - *: (w - 550.0) % 25.0
  refl=. refl + 0.48 % (1.0 + ^ - (w - 710.0) % 22.0)
  nir_mask=. (w >: 760.0) *. (w <: 1300.0)
  refl=. (0.52 + 0.03 * 2&o. w * 0.015) nir_mask} refl
  refl=. refl - 0.22 * ^ - *: (w - 1450.0) % 55.0
  refl=. refl - 0.28 * ^ - *: (w - 1940.0) % 65.0
  0.01 >. 0.65 <. refl
)

NB. Conifer Pine Needle (Douglas Fir / Loblolly)
spec_conifer_needle=: 3 : 0
  w=. y
  refl=. 0.03 + 0.05 * ^ - *: (w - 550.0) % 25.0
  refl=. refl + 0.32 % (1.0 + ^ - (w - 715.0) % 25.0)
  nir_mask=. (w >: 760.0) *. (w <: 1300.0)
  refl=. (0.35 + 0.02 * 2&o. w * 0.015) nir_mask} refl
  refl=. refl - 0.16 * ^ - *: (w - 1450.0) % 50.0
  refl=. refl - 0.20 * ^ - *: (w - 1940.0) % 60.0
  0.01 >. 0.45 <. refl
)

NB. Dry Senescent Corn / Stover (Lignin & Cellulose)
spec_senescent_maize=: 3 : 0
  w=. y
  refl=. 0.12 + 0.38 * (w - 400.0) % 2100.0
  refl=. refl - 0.09 * ^ - *: (w - 2100.0) % 75.0
  refl=. refl - 0.07 * ^ - *: (w - 2310.0) % 65.0
  0.03 >. 0.58 <. refl
)

NB. Soybean Crop Canopy
spec_soybean_canopy=: 3 : 0
  w=. y
  refl=. 0.05 + 0.09 * ^ - *: (w - 555.0) % 25.0
  refl=. refl + 0.50 % (1.0 + ^ - (w - 708.0) % 20.0)
  nir_mask=. (w >: 750.0) *. (w <: 1300.0)
  refl=. (0.55 + 0.03 * 2&o. w * 0.012) nir_mask} refl
  refl=. refl - 0.24 * ^ - *: (w - 1450.0) % 55.0
  refl=. refl - 0.30 * ^ - *: (w - 1940.0) % 65.0
  0.01 >. 0.68 <. refl
)

NB. Vineyard Grapevine Canopy
spec_vineyard_grape=: 3 : 0
  w=. y
  refl=. 0.04 + 0.07 * ^ - *: (w - 548.0) % 22.0
  refl=. refl + 0.42 % (1.0 + ^ - (w - 712.0) % 24.0)
  nir_mask=. (w >: 750.0) *. (w <: 1300.0)
  refl=. (0.46 + 0.02 * 2&o. w * 0.01) nir_mask} refl
  refl=. refl - 0.18 * ^ - *: (w - 1450.0) % 50.0
  refl=. refl - 0.24 * ^ - *: (w - 1940.0) % 60.0
  0.01 >. 0.58 <. refl
)

NB. Wetland Sedge (Carex stricta)
spec_wetland_sedge=: 3 : 0
  w=. y
  refl=. 0.04 + 0.06 * ^ - *: (w - 550.0) % 25.0
  refl=. refl + 0.38 % (1.0 + ^ - (w - 710.0) % 22.0)
  nir_mask=. (w >: 750.0) *. (w <: 1300.0)
  refl=. (0.42 + 0.02 * 2&o. w * 0.01) nir_mask} refl
  refl=. refl - 0.26 * ^ - *: (w - 1450.0) % 55.0
  refl=. refl - 0.32 * ^ - *: (w - 1940.0) % 65.0
  0.01 >. 0.50 <. refl
)

NB. Sphagnum Peat Moss
spec_sphagnum_moss=: 3 : 0
  w=. y
  refl=. 0.06 + 0.07 * ^ - *: (w - 560.0) % 30.0
  refl=. refl + 0.25 % (1.0 + ^ - (w - 715.0) % 30.0)
  refl=. refl - 0.22 * ^ - *: (w - 1450.0) % 60.0
  refl=. refl - 0.28 * ^ - *: (w - 1940.0) % 70.0
  0.02 >. 0.38 <. refl
)

NB. Freshwater Microalgae (Chlorella / Microcystis)
spec_freshwater_algae=: 3 : 0
  w=. y
  refl=. 0.02 + 0.04 * ^ - *: (w - 550.0) % 20.0
  NB. Phycocyanin dip at 620 nm, chlorophyll-a at 675 nm
  refl=. refl - 0.015 * ^ - *: (w - 620.0) % 15.0
  refl=. refl + 0.08 % (1.0 + ^ - (w - 700.0) % 15.0)
  NB. Water absorption dominates beyond 750 nm
  refl=. refl * ^ - (w - 700.0) % 180.0
  0.005 >. 0.18 <. refl
)

NB. Biological Soil Crust / Lichen
spec_soil_crust_lichen=: 3 : 0
  w=. y
  refl=. 0.08 + 0.22 * (w - 400.0) % 2100.0
  refl=. refl + 0.03 * ^ - *: (w - 680.0) % 35.0
  refl=. refl - 0.06 * ^ - *: (w - 2100.0) % 60.0
  0.04 >. 0.35 <. refl
)

NB. Turf Grass (Kentucky Bluegrass)
spec_turf_grass=: 3 : 0
  w=. y
  refl=. 0.04 + 0.08 * ^ - *: (w - 550.0) % 25.0
  refl=. refl + 0.45 % (1.0 + ^ - (w - 710.0) % 22.0)
  nir_mask=. (w >: 750.0) *. (w <: 1300.0)
  refl=. (0.50 + 0.02 * 2&o. w * 0.01) nir_mask} refl
  refl=. refl - 0.20 * ^ - *: (w - 1450.0) % 55.0
  refl=. refl - 0.26 * ^ - *: (w - 1940.0) % 65.0
  0.01 >. 0.62 <. refl
)

NB. Build full Botanical library matrix (224 bands x 10 species)
build_botanical_library=: 3 : 0
  wvl=. y
  b1=.  spec_deciduous_broadleaf wvl
  b2=.  spec_conifer_needle wvl
  b3=.  spec_senescent_maize wvl
  b4=.  spec_soybean_canopy wvl
  b5=.  spec_vineyard_grape wvl
  b6=.  spec_wetland_sedge wvl
  b7=.  spec_sphagnum_moss wvl
  b8=.  spec_freshwater_algae wvl
  b9=.  spec_soil_crust_lichen wvl
  b10=. spec_turf_grass wvl

  b1 ,. b2 ,. b3 ,. b4 ,. b5 ,. b6 ,. b7 ,. b8 ,. b9 ,. b10
)


albedo_usgs_vegetation_build_botanical_library =: build_botanical_library
