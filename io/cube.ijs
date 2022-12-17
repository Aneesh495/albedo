NB. ==============================================================================
NB. Albedo - Hyperspectral Data Cube Container & Slicing
NB.
NB. High-efficiency container for 3D radiance and reflectance cubes (B x H x W).
NB. Provides spatial/spectral subcube slicing, spatial masking, wavelength matching,
NB. and false-color RGB composite generation (e.g. Color Infrared CIR, True Color).
NB. ==============================================================================

cocurrent 'base'

load 'core/tensor.ijs'

NB. ------------------------------------------------------------------------------
NB. Cube Structure & Construction
NB. ------------------------------------------------------------------------------

NB. Create cube box container: (data_3d; wavelengths; metadata_dict)
create_cube=: 3 : 0
  data_3d=. 0 {:: y
  wvl=. 1 {:: y
  meta=. 2 {:: y
  'B H W'=. $ data_3d
  assert. B = # wvl [ 'Band count must match wavelength vector length'
  data_3d; wvl; < meta
)

NB. Accessors
cube_data=: 0&{::
cube_wvl=: 1&{::
cube_meta=: 2&{::
cube_shape=: $@cube_data

NB. ------------------------------------------------------------------------------
NB. Wavelength & Spectral Band Lookup
NB. ------------------------------------------------------------------------------

NB. Find band index nearest to a target wavelength (nm)
NB. Dyadic: target_nm find_nearest_band wvl_vector
find_nearest_band=: 4 : 0
  target_nm=. x
  wvl=. y
  diffs=. | wvl - target_nm
  (<: +/ diffs = <./ diffs) { (I. diffs = <./ diffs)
)

NB. ------------------------------------------------------------------------------
NB. Subcube Slicing
NB. ------------------------------------------------------------------------------

NB. Spatial slice: [row_start, row_end, col_start, col_end]
NB. Dyadic: [r0, r1, c0, c1] slice_spatial cube
slice_spatial=: 4 : 0
  'r0 r1 c0 c1'=. x
  data=. 0 {:: y
  wvl=. 1 {:: y
  meta=. 2 {:: y
  'B H W'=. $ data
  rows=. r0 + i. (r1 - r0)
  cols=. c0 + i. (c1 - c0)
  sub_data=. (<a:; rows; cols) { data
  sub_data; wvl; < meta
)

NB. Spectral slice: [wvl_min, wvl_max]
NB. Dyadic: [wmin, wmax] slice_spectral cube
slice_spectral=: 4 : 0
  'wmin wmax'=. x
  data=. 0 {:: y
  wvl=. 1 {:: y
  meta=. 2 {:: y
  mask=. (wvl >: wmin) *. (wvl <: wmax)
  band_indices=. I. mask
  assert. 0 < # band_indices [ 'No bands found in specified wavelength range'
  sub_data=. band_indices { data
  sub_wvl=. band_indices { wvl
  sub_data; sub_wvl; < meta
)

NB. ------------------------------------------------------------------------------
NB. False-Color RGB Composite Generation
NB. ------------------------------------------------------------------------------

NB. Linear 2% to 98% percentile contrast stretch of a 2D band
contrast_stretch_2d=: 3 : 0
  band=. y
  vals=. , band
  s_vals=. /:~ vals
  n=. # s_vals
  p02=. (<. 0.02 * n) { s_vals
  p98=. (<. 0.98 * n) { s_vals
  denom=. (p98 - p02) >. 1e_6
  stretched=. (band - p02) % denom
  0.0 >. 1.0 <. stretched
)

NB. Generate 3-channel RGB image (3 x H x W) given red, green, blue target wavelengths
NB. Dyadic: (red_nm; green_nm; blue_nm) make_rgb_composite cube
make_rgb_composite=: 4 : 0
  'r_nm g_nm b_nm'=. x
  data=. 0 {:: y
  wvl=. 1 {:: y
  meta=. 2 {:: y

  idx_r=. r_nm find_nearest_band wvl
  idx_g=. g_nm find_nearest_band wvl
  idx_b=. b_nm find_nearest_band wvl

  band_r=. contrast_stretch_2d (idx_r { data)
  band_g=. contrast_stretch_2d (idx_g { data)
  band_b=. contrast_stretch_2d (idx_b { data)

  (,: band_r) , (,: band_g) , (,: band_b)
)

NB. Standard Natural Color (True Color) RGB: 650 nm, 550 nm, 470 nm
true_color_rgb=: 3 : 0
  (650.0; 550.0; 470.0) make_rgb_composite y
)

NB. Standard Color Infrared (CIR) Composite: 850 nm (NIR), 650 nm (Red), 550 nm (Green)
NB. Vigorous vegetation appears bright red.
color_infrared_cir=: 3 : 0
  (850.0; 650.0; 550.0) make_rgb_composite y
)


albedo_cube_create_cube =: create_cube
albedo_cube_cube_data =: cube_data
albedo_cube_cube_wvl =: cube_wvl
albedo_cube_slice_spatial =: slice_spatial
albedo_cube_slice_spectral =: slice_spectral
albedo_cube_color_infrared_cir =: color_infrared_cir
albedo_cube_find_nearest_band =: find_nearest_band
