NB. ==============================================================================
NB. Albedo - Physically Realistic Synthetic Hyperspectral Scene Generator
NB.
NB. Synthesizes full 224-band hyperspectral scenes with realistic terrain:
NB.   - River channel (clear water)
NB.   - Agricultural cropland (green vegetation with high NIR reflectance)
NB.   - Hydrothermal mineral outcrop (kaolinite clay Al-OH doublet at 2205 nm)
NB.   - Carbonate formation (calcite absorption at 2340 nm)
NB.   - Urban roadway (impervious asphalt)
NB.   - Subpixel methane (CH4) gas plume leak
NB. Simulates atmospheric solar radiative transfer (solar irradiance, path radiance,
NB. and sensor Gaussian noise).
NB. ==============================================================================

cocurrent 'base'

load 'io/spectral_lib.ijs'
load 'io/cube.ijs'
load 'atmos/solar.ijs'
load 'core/linalg.ijs'

NB. ------------------------------------------------------------------------------
NB. Scene Landscape Synthesis
NB. ------------------------------------------------------------------------------

NB. Generate 224-band synthetic scene of spatial dimensions (H, W)
NB. Dyadic: (H; W; snr_db) generate_scene ''
generate_scene=: 3 : 0
  'H W snr_db'=. y
  wvl=. albedo_spectral_lib_standard_wavelengths_224 ''
  B=. # wvl

  NB. Endmember signatures
  sig_veg=.     albedo_spectral_lib_ref_green_vegetation wvl
  sig_water=.   albedo_spectral_lib_ref_clear_water wvl
  sig_soil=.    albedo_spectral_lib_ref_kaolinite_soil wvl
  sig_calcite=. albedo_spectral_lib_ref_calcite_mineral wvl
  sig_asphalt=. albedo_spectral_lib_ref_urban_asphalt wvl
  sig_methane=. albedo_spectral_lib_ref_methane_plume wvl

  NB. Initialize ground-truth abundance cubes: 5 primary endmembers
  NB. [0: Veg, 1: Water, 2: Soil/Kaolinite, 3: Calcite, 4: Asphalt]
  A_gt=. (5, H, W) $ 0.0

  NB. Terrain zoning
  for_h. i. H do.
    for_w. i. W do.
      NB. River channel: diagonal strip
      is_river=. (| h - w) <: (<. 0.08 * H)
      NB. Roadway: vertical strip
      is_road=. (| w - (<. 0.8 * W)) <: 1
      NB. Calcite quarry: circle near bottom right
      dist_quarry=. %: (*: h - (<. 0.75 * H)) + (*: w - (<. 0.65 * W))
      is_calcite=. dist_quarry <: (<. 0.15 * H)
      NB. Soil outcrop: quadrant top right
      is_soil=. (h < (<. 0.5 * H)) *. (w > (<. 0.5 * W))

      if. is_river do.
        A_gt=. 0.90 (<0, h, w)} A_gt
        A_gt=. 0.10 (<1, h, w)} A_gt
      elseif. is_road do.
        A_gt=. 0.85 (<4, h, w)} A_gt
        A_gt=. 0.15 (<2, h, w)} A_gt
      elseif. is_calcite do.
        A_gt=. 0.80 (<3, h, w)} A_gt
        A_gt=. 0.20 (<2, h, w)} A_gt
      elseif. is_soil do.
        A_gt=. 0.70 (<2, h, w)} A_gt
        A_gt=. 0.30 (<0, h, w)} A_gt
      else.
        NB. Agricultural cropland
        A_gt=. 0.85 (<0, h, w)} A_gt
        A_gt=. 0.15 (<2, h, w)} A_gt
      end.
    end.
  end.

  NB. Endmember matrix M of shape (B, 5)
  M=. |: sig_veg ,: sig_water ,: sig_soil ,: sig_calcite ,: sig_asphalt

  NB. Linear spectral mixing: Reflectance = M * A_gt
  Refl_mat=. M albedo_linalg_matmul ((5, H * W) $ A_gt)
  Refl_cube=. (B, H, W) $ Refl_mat

  NB. Inject subpixel methane plume at (h=H/2, w=W/4) with 2x2 dispersion
  plume_h=. <. 0.5 * H
  plume_w=. <. 0.25 * W
  for_dh. 0 1 do.
    for_dw. 0 1 do.
      orig_spec=. (<a:; (plume_h + dh); (plume_w + dw)) { Refl_cube
      NB. Plume transmission: 85% background + 15% methane gas absorption
      gas_spec=. (0.85 * orig_spec) + (0.15 * sig_methane)
      Refl_cube=. gas_spec (<a:; (plume_h + dh); (plume_w + dw))} Refl_cube
    end.
  end.

  NB. Add Gaussian instrument noise
  noise_sigma=. 10.0 ^ (_0.05 * snr_db)
  noise=. (B, H, W) $ 2.0 * ((? (B * H * W) $ 1000) % 1000.0) - 0.5
  Refl_noisy=. 0.0 >. 1.0 <. (Refl_cube + (noise_sigma * noise))

  NB. Forward radiative transfer to Radiance (W / (m^2 * sr * nm))
  doy=. 120
  sza_deg=. 35.0
  d2=. albedo_solar_earth_sun_dist_sq doy
  cos_sza=. albedo_solar_cos_solar_zenith sza_deg
  E0=. albedo_solar_solar_irradiance_e0 wvl

  Rad_cube=. (B, H, W) $ 0.0
  for_b. i. B do.
    scale=. ((b { E0) * cos_sza) % (3.141592653589793 * d2)
    NB. Add atmospheric Rayleigh path radiance haze
    haze=. 0.05 * (b { E0) * (400.0 % (b { wvl)) ^ 4
    rad_band=. (scale * b { Refl_noisy) + haze
    Rad_cube=. rad_band (<b)} Rad_cube
  end.

  NB. Returns (Rad_cube; Refl_noisy; A_gt; M; wvl; (plume_h, plume_w))
  Rad_cube; Refl_noisy; A_gt; M; wvl; (plume_h, plume_w)
)


albedo_sim_generate_scene =: generate_scene
