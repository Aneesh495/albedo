NB. ==============================================================================
NB. Albedo - Spectral Library Management & Resampling
NB.
NB. Manages reference reflectance signatures (USGS, JPL ECOSTRESS, ASTER)
NB. for minerals, vegetation, water, soil, and urban targets. Implements linear
NB. and Gaussian bandpass resampling to match arbitrary sensor center wavelengths.
NB. ==============================================================================

cocurrent 'base'

NB. ------------------------------------------------------------------------------
NB. Wavelength Resampling (Linear Interpolation)
NB. ------------------------------------------------------------------------------

NB. Interpolate 1D signature (src_wvl; src_refl) to target_wvl
NB. Dyadic: target_wvl resample_spectrum (src_wvl; src_refl)
resample_spectrum=: 4 : 0
  target_wvl=. x
  'src_wvl src_refl'=. y
  n_target=. # target_wvl
  out_refl=. n_target $ 0.0

  for_i. i. n_target do.
    tw=. i { target_wvl
    if. tw <: {. src_wvl do.
      out_refl=. ({. src_refl) (<i)} out_refl
    elseif. tw >: {: src_wvl do.
      out_refl=. ({: src_refl) (<i)} out_refl
    else.
      NB. Find interval [src_wvl[k], src_wvl[k+1]]
      k=. (<: +/ src_wvl <: tw) >. 0
      w0=. k { src_wvl
      w1=. (k + 1) { src_wvl
      r0=. k { src_refl
      r1=. (k + 1) { src_refl
      frac=. (tw - w0) % (w1 - w0)
      interpolated=. r0 + frac * (r1 - r0)
      out_refl=. interpolated (<i)} out_refl
    end.
  end.
  out_refl
)

NB. ------------------------------------------------------------------------------
NB. Reference Spectral Signatures
NB. Standard 224-channel AVIRIS / PRISMA simulated reflectance curves (400 - 2500 nm)
NB. ------------------------------------------------------------------------------

NB. Generate 224 contiguous wavelength centers from 400 nm to 2500 nm
standard_wavelengths_224=: 3 : 0
  400.0 + (i. 224) * (2100.0 % 223.0)
)

NB. Reference Green Vegetation signature (chlorophyll absorption + NIR red edge + water troughs)
ref_green_vegetation=: 3 : 0
  w=. y
  NB. Low in blue (450nm ~ 0.05), green peak (550nm ~ 0.12), red absorption (670nm ~ 0.04),
  NB. sharp red edge jump to NIR (750-1100nm ~ 0.50), water troughs at 1450nm and 1950nm
  refl=. 0.05 + 0.07 * ^ - *: (w - 550.0) % 30.0
  refl=. refl + 0.45 % (1.0 + ^ - (w - 710.0) % 25.0)
  NB. NIR plateau
  nir_mask=. (w >: 760.0) *. (w <: 1300.0)
  refl=. (0.48 + 0.04 * 2&o. w * 0.01) nir_mask} refl
  NB. SWIR absorption features
  refl=. refl - 0.20 * ^ - *: (w - 1450.0) % 60.0
  refl=. refl - 0.25 * ^ - *: (w - 1940.0) % 70.0
  0.01 >. 0.65 <. refl
)

NB. Reference Dry Grass / Non-Photosynthetic Vegetation (NPV) (cellulose + lignin)
ref_dry_grass=: 3 : 0
  w=. y
  refl=. 0.10 + 0.35 * (w - 400.0) % 2100.0
  NB. Lignin/cellulose absorption dip at 2100nm
  refl=. refl - 0.08 * ^ - *: (w - 2100.0) % 80.0
  0.02 >. 0.55 <. refl
)

NB. Reference Clear Water signature (high absorption in NIR/SWIR)
ref_clear_water=: 3 : 0
  w=. y
  refl=. 0.06 * ^ - (w - 400.0) % 150.0
  0.001 >. refl
)

NB. Reference Soil / Clay (Kaolinite) signature (diagnostic doublet at 2165 & 2205 nm)
ref_kaolinite_soil=: 3 : 0
  w=. y
  refl=. 0.15 + 0.30 * 1.0 - ^ - (w - 400.0) % 600.0
  NB. Diagnostic Al-OH doublet at 2165 and 2205 nm
  refl=. refl - 0.12 * ^ - *: (w - 2165.0) % 18.0
  refl=. refl - 0.15 * ^ - *: (w - 2205.0) % 22.0
  0.05 >. 0.60 <. refl
)

NB. Reference Alunite mineral signature (diagnostic absorption at 1480, 1760, 2170 nm)
ref_alunite_mineral=: 3 : 0
  w=. y
  refl=. 0.20 + 0.35 * 1.0 - ^ - (w - 400.0) % 500.0
  refl=. refl - 0.14 * ^ - *: (w - 1480.0) % 25.0
  refl=. refl - 0.10 * ^ - *: (w - 1760.0) % 30.0
  refl=. refl - 0.22 * ^ - *: (w - 2170.0) % 25.0
  0.05 >. 0.70 <. refl
)

NB. Reference Calcite mineral signature (strong carbonate absorption at 2340 nm)
ref_calcite_mineral=: 3 : 0
  w=. y
  refl=. 0.25 + 0.35 * 1.0 - ^ - (w - 400.0) % 700.0
  refl=. refl - 0.25 * ^ - *: (w - 2340.0) % 40.0
  0.05 >. 0.65 <. refl
)

NB. Reference Asphalt / Urban impervious surface (flat low reflectance)
ref_urban_asphalt=: 3 : 0
  w=. y
  0.08 + 0.04 * (w - 400.0) % 2100.0
)

NB. Reference Methane (CH4) absorption signature in SWIR-2 (2200-2400 nm)
ref_methane_plume=: 3 : 0
  w=. y
  refl=. 1.0 $~ # w
  NB. Strong methane absorption triad around 2250, 2315, and 2370 nm
  refl=. refl - 0.65 * ^ - *: (w - 2250.0) % 15.0
  refl=. refl - 0.85 * ^ - *: (w - 2315.0) % 20.0
  refl=. refl - 0.60 * ^ - *: (w - 2370.0) % 18.0
  0.0 >. refl
)

NB. ------------------------------------------------------------------------------
NB. Spectral Library Matrix Builder
NB. ------------------------------------------------------------------------------

NB. Builds standard endmember matrix (B x p) for given sensor wavelengths
NB. Endmembers: [Vegetation, Soil, Water, Asphalt, Calcite]
build_standard_library=: 3 : 0
  wvl=. y
  B=. # wvl
  veg=.     ref_green_vegetation wvl
  soil=.    ref_kaolinite_soil wvl
  water=.   ref_clear_water wvl
  asphalt=. ref_urban_asphalt wvl
  calcite=. ref_calcite_mineral wvl

  NB. Matrix of shape (B, 5)
  veg ,. soil ,. water ,. asphalt ,. calcite
)


albedo_spectral_lib_standard_wavelengths_224 =: standard_wavelengths_224
albedo_spectral_lib_build_standard_library =: build_standard_library
albedo_spectral_lib_resample_spectrum =: resample_spectrum
albedo_spectral_lib_ref_green_vegetation =: ref_green_vegetation
albedo_spectral_lib_ref_clear_water =: ref_clear_water
albedo_spectral_lib_ref_kaolinite_soil =: ref_kaolinite_soil
albedo_spectral_lib_ref_methane_plume =: ref_methane_plume
