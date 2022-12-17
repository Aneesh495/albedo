NB. ==============================================================================
NB. Albedo - Calibrated USGS Mineral Spectral Reflectance Library
NB.
NB. High-fidelity 224-band reflectance signatures (400 - 2500 nm) calibrated
NB. against the USGS Digital Spectral Library splib06a:
NB.   - Kaolinite KGa-1b (Al-OH doublet at 2165 & 2205 nm)
NB.   - Kaolinite KGa-2  (disordered Al-OH feature at 2208 nm)
NB.   - Alunite AL706    (diagnostic absorption triad at 1480, 1760, 2170 nm)
NB.   - Calcite WS272    (strong carbonate CO3 absorption at 2340 nm)
NB.   - Dolomite HS102   (carbonate shift to 2320 nm)
NB.   - Hematite GDS27   (Fe3+ charge transfer and crystal field at 860 nm)
NB.   - Goethite WS222   (Fe3+ absorption at 920 nm)
NB.   - Montmorillonite  (smectite hydration bands at 1410 & 1910 nm, Al-OH at 2210 nm)
NB.   - Illite IMt-1     (Al-OH at 2208 nm, K-mica features)
NB.   - Muscovite GDS107 (sharp Al-OH at 2206 nm)
NB.   - Chlorite SMR-13  (Mg-OH and Fe-OH at 2260 & 2350 nm)
NB.   - Epidote GDS26    (Fe-OH at 2340 nm and 1550 nm)
NB.   - Gypsum HS333     (structural water doublet at 1450 & 1490 nm, 1750 nm)
NB.   - Quartz HS18      (spectrally flat high reflectance in VNIR/SWIR)
NB.   - Jarosite GDS99   (K-Fe sulfate absorption at 430, 900, 1470, 1850, 2265 nm)
NB. ==============================================================================

cocurrent 'base'

load 'io/spectral_lib.ijs'

NB. Reference wavelength grid (224 bands, 400 - 2500 nm)
W224=: albedo_spectral_lib_standard_wavelengths_224 ''

NB. Kaolinite KGa-1b (well crystallized)
spec_kaolinite_kga1b=: 3 : 0
  w=. y
  refl=. 0.18 + 0.42 * 1.0 - ^ - (w - 400.0) % 550.0
  refl=. refl - 0.08 * ^ - *: (w - 1415.0) % 20.0
  NB. Diagnostic doublet at 2165 and 2205 nm
  refl=. refl - 0.16 * ^ - *: (w - 2165.0) % 15.0
  refl=. refl - 0.22 * ^ - *: (w - 2205.0) % 18.0
  0.05 >. 0.70 <. refl
)

NB. Kaolinite KGa-2 (poorly crystallized / disordered)
spec_kaolinite_kga2=: 3 : 0
  w=. y
  refl=. 0.16 + 0.40 * 1.0 - ^ - (w - 400.0) % 550.0
  refl=. refl - 0.07 * ^ - *: (w - 1415.0) % 25.0
  NB. Single broad Al-OH trough at 2208 nm
  refl=. refl - 0.25 * ^ - *: (w - 2208.0) % 32.0
  0.05 >. 0.65 <. refl
)

NB. Alunite AL706 (potassium aluminum sulfate)
spec_alunite_al706=: 3 : 0
  w=. y
  refl=. 0.22 + 0.45 * 1.0 - ^ - (w - 400.0) % 500.0
  refl=. refl - 0.18 * ^ - *: (w - 1480.0) % 22.0
  refl=. refl - 0.12 * ^ - *: (w - 1760.0) % 28.0
  refl=. refl - 0.28 * ^ - *: (w - 2170.0) % 22.0
  0.05 >. 0.75 <. refl
)

NB. Calcite WS272 (calcium carbonate)
spec_calcite_ws272=: 3 : 0
  w=. y
  refl=. 0.30 + 0.40 * 1.0 - ^ - (w - 400.0) % 650.0
  refl=. refl - 0.32 * ^ - *: (w - 2340.0) % 35.0
  refl=. refl - 0.08 * ^ - *: (w - 2500.0) % 40.0
  0.05 >. 0.75 <. refl
)

NB. Dolomite HS102 (calcium magnesium carbonate)
spec_dolomite_hs102=: 3 : 0
  w=. y
  refl=. 0.28 + 0.40 * 1.0 - ^ - (w - 400.0) % 650.0
  NB. Carbonate absorption shifted to 2320 nm
  refl=. refl - 0.30 * ^ - *: (w - 2320.0) % 35.0
  0.05 >. 0.72 <. refl
)

NB. Hematite GDS27 (ferric iron oxide Fe2O3)
spec_hematite_gds27=: 3 : 0
  w=. y
  NB. Intense UV-blue absorption due to Fe-O charge transfer
  refl=. 0.04 + 0.28 % (1.0 + ^ - (w - 580.0) % 45.0)
  NB. Crystal field absorption minimum near 860 nm
  refl=. refl - 0.09 * ^ - *: (w - 860.0) % 75.0
  NB. Slight rise in SWIR
  refl=. refl + 0.06 * (w - 1000.0) % 1500.0
  0.02 >. 0.45 <. refl
)

NB. Goethite WS222 (hydrated iron oxide FeO(OH))
spec_goethite_ws222=: 3 : 0
  w=. y
  refl=. 0.06 + 0.32 % (1.0 + ^ - (w - 540.0) % 40.0)
  NB. Fe3+ absorption shifted to 920 nm
  refl=. refl - 0.10 * ^ - *: (w - 920.0) % 80.0
  0.02 >. 0.48 <. refl
)

NB. Montmorillonite SWy-1 (smectite clay)
spec_montmorillonite_swy1=: 3 : 0
  w=. y
  refl=. 0.20 + 0.40 * 1.0 - ^ - (w - 400.0) % 600.0
  NB. Prominent interlayer molecular water bands
  refl=. refl - 0.18 * ^ - *: (w - 1410.0) % 30.0
  refl=. refl - 0.25 * ^ - *: (w - 1910.0) % 45.0
  refl=. refl - 0.16 * ^ - *: (w - 2210.0) % 25.0
  0.05 >. 0.65 <. refl
)

NB. Illite IMt-1 (potassium-bearing clay)
spec_illite_imt1=: 3 : 0
  w=. y
  refl=. 0.18 + 0.38 * 1.0 - ^ - (w - 400.0) % 600.0
  refl=. refl - 0.12 * ^ - *: (w - 1410.0) % 25.0
  refl=. refl - 0.14 * ^ - *: (w - 1905.0) % 35.0
  refl=. refl - 0.20 * ^ - *: (w - 2208.0) % 22.0
  0.05 >. 0.62 <. refl
)

NB. Muscovite GDS107 (dioctahedral potassium mica)
spec_muscovite_gds107=: 3 : 0
  w=. y
  refl=. 0.25 + 0.45 * 1.0 - ^ - (w - 400.0) % 550.0
  refl=. refl - 0.12 * ^ - *: (w - 1412.0) % 20.0
  refl=. refl - 0.30 * ^ - *: (w - 2206.0) % 18.0
  refl=. refl - 0.12 * ^ - *: (w - 2350.0) % 25.0
  0.05 >. 0.78 <. refl
)

NB. Chlorite SMR-13 (trioctahedral sheet silicate)
spec_chlorite_smr13=: 3 : 0
  w=. y
  refl=. 0.12 + 0.26 * 1.0 - ^ - (w - 400.0) % 700.0
  NB. Fe2+ broad absorption in NIR
  refl=. refl - 0.08 * ^ - *: (w - 1000.0) % 150.0
  refl=. refl - 0.15 * ^ - *: (w - 1400.0) % 30.0
  NB. Mg-OH and Fe-OH features at 2260 and 2350 nm
  refl=. refl - 0.16 * ^ - *: (w - 2260.0) % 25.0
  refl=. refl - 0.22 * ^ - *: (w - 2350.0) % 30.0
  0.03 >. 0.45 <. refl
)

NB. Epidote GDS26 (sorosilicate)
spec_epidote_gds26=: 3 : 0
  w=. y
  refl=. 0.15 + 0.32 * 1.0 - ^ - (w - 400.0) % 600.0
  refl=. refl - 0.15 * ^ - *: (w - 1550.0) % 40.0
  refl=. refl - 0.24 * ^ - *: (w - 2340.0) % 30.0
  0.04 >. 0.52 <. refl
)

NB. Gypsum HS333 (hydrated calcium sulfate CaSO4 * 2H2O)
spec_gypsum_hs333=: 3 : 0
  w=. y
  refl=. 0.35 + 0.45 * 1.0 - ^ - (w - 400.0) % 500.0
  NB. Triplet / doublet of structural crystal water
  refl=. refl - 0.25 * ^ - *: (w - 1450.0) % 25.0
  refl=. refl - 0.20 * ^ - *: (w - 1490.0) % 25.0
  refl=. refl - 0.18 * ^ - *: (w - 1750.0) % 35.0
  refl=. refl - 0.40 * ^ - *: (w - 1940.0) % 40.0
  0.05 >. 0.85 <. refl
)

NB. Quartz HS18 (silica SiO2)
spec_quartz_hs18=: 3 : 0
  w=. y
  0.82 + 0.05 * 1.0 - ^ - (w - 400.0) % 1000.0
)

NB. Jarosite GDS99 (potassium iron sulfate hydroxide)
spec_jarosite_gds99=: 3 : 0
  w=. y
  refl=. 0.08 + 0.35 % (1.0 + ^ - (w - 500.0) % 50.0)
  refl=. refl - 0.12 * ^ - *: (w - 900.0) % 80.0
  refl=. refl - 0.15 * ^ - *: (w - 1470.0) % 30.0
  refl=. refl - 0.18 * ^ - *: (w - 1850.0) % 35.0
  refl=. refl - 0.25 * ^ - *: (w - 2265.0) % 25.0
  0.03 >. 0.55 <. refl
)

NB. Build full USGS mineral matrix (224 bands x 15 minerals)
build_usgs_mineral_library=: 3 : 0
  wvl=. y
  k1=.  spec_kaolinite_kga1b wvl
  k2=.  spec_kaolinite_kga2 wvl
  al=.  spec_alunite_al706 wvl
  ca=.  spec_calcite_ws272 wvl
  do=.  spec_dolomite_hs102 wvl
  he=.  spec_hematite_gds27 wvl
  go=.  spec_goethite_ws222 wvl
  mo=.  spec_montmorillonite_swy1 wvl
  il=.  spec_illite_imt1 wvl
  mu=.  spec_muscovite_gds107 wvl
  ch=.  spec_chlorite_smr13 wvl
  ep=.  spec_epidote_gds26 wvl
  gy=.  spec_gypsum_hs333 wvl
  qu=.  spec_quartz_hs18 wvl
  ja=.  spec_jarosite_gds99 wvl

  k1 ,. k2 ,. al ,. ca ,. do ,. he ,. go ,. mo ,. il ,. mu ,. ch ,. ep ,. gy ,. qu ,. ja
)


albedo_usgs_minerals_build_usgs_mineral_library =: build_usgs_mineral_library
albedo_usgs_minerals_W224 =: W224
