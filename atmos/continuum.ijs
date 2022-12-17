NB. ==============================================================================
NB. Albedo - Convex Hull Continuum Removal
NB.
NB. Normalizes reflectance spectra by computing the upper convex envelope (continuum)
NB. connecting local reflectance maxima. Isolates diagnostic absorption features,
NB. band center wavelength, maximum band depth, full width at half depth (FWHD),
NB. and integrated absorption band area.
NB. ==============================================================================

cocurrent 'base'

NB. ------------------------------------------------------------------------------
NB. Upper Convex Hull Calculation for 1D Spectral Curve
NB. ------------------------------------------------------------------------------

NB. Computes upper convex hull of points (wvl, refl) using Graham scan variant
upper_convex_hull=: 4 : 0
  wvl=. x
  refl=. y
  n=. # wvl
  assert. n >: 2 [ 'At least two points required for convex hull'

  hull=. n $ 0.0
  hull=. ({. refl) (0)} hull
  hull=. ({: refl) (<: n)} hull

  NB. Identify upper convex envelope vertices
  NB. Start with first point
  hull_indices=. ,: 0
  curr=. 0

  while. curr < (n - 1) do.
    max_slope=. __
    best_next=. n - 1

    for_next. (curr + 1) + i. (n - curr - 1) do.
      dx=. (next { wvl) - (curr { wvl)
      dy=. (next { refl) - (curr { refl)
      slope=. dy % dx
      if. slope > max_slope do.
        max_slope=. slope
        best_next=. next
      end.
    end.

    hull_indices=. hull_indices , best_next
    curr=. best_next
  end.

  NB. Interpolate hull between selected vertex indices
  for_k. i. (<: # hull_indices) do.
    idx0=. k { hull_indices
    idx1=. (k + 1) { hull_indices
    w0=. idx0 { wvl
    w1=. idx1 { wvl
    r0=. idx0 { refl
    r1=. idx1 { refl
    slope=. (r1 - r0) % (w1 - w0)

    for_j. idx0 + i. (1 + idx1 - idx0) do.
      w_j=. j { wvl
      hull=. (r0 + slope * (w_j - w0)) (<j)} hull
    end.
  end.

  NB. Ensure hull is never strictly below original reflectance
  hull >. refl
)

NB. ------------------------------------------------------------------------------
NB. Continuum Removal & Feature Characterization
NB. ------------------------------------------------------------------------------

NB. Removes continuum from 1D reflectance curve: rho_cr = rho / hull
NB. Values are in [0, 1] where 1 is continuum and < 1 indicates absorption
remove_continuum=: 4 : 0
  wvl=. x
  refl=. y
  hull=. wvl upper_convex_hull refl
  refl % (hull >. 1e_6)
)

NB. Characterizes single absorption feature in interval [w_start, w_end]
NB. Returns: (center_wvl; max_depth; band_area; fwhd)
characterize_absorption=: 4 : 0
  'w_start w_end'=. x
  'wvl refl'=. y

  mask=. (wvl >: w_start) *. (wvl <: w_end)
  sub_indices=. I. mask
  if. 0 = # sub_indices do. 0.0; 0.0; 0.0; 0.0 return. end.

  sub_wvl=. sub_indices { wvl
  sub_refl=. sub_indices { refl

  cr=. sub_wvl remove_continuum sub_refl
  depths=. 1.0 - cr

  max_depth=. >./ depths
  min_idx=. (<: +/ depths = max_depth) { (I. depths = max_depth)
  center_wvl=. min_idx { sub_wvl

  NB. Band area: trapezoidal integration of (1 - cr) d(lambda)
  dw=. 1 }. sub_wvl - (_1 }. sub_wvl)
  avg_depth=. -: (1 }. depths) + (_1 }. depths)
  band_area=. +/ avg_depth * dw

  NB. Full Width at Half Depth (FWHD)
  half_depth=. 0.5 * max_depth
  half_mask=. depths >: half_depth
  if. (+/ half_mask) >: 2 do.
    hw_indices=. I. half_mask
    fwhd=. ({: hw_indices { sub_wvl) - ({. hw_indices { sub_wvl)
  else.
    fwhd=. 0.0
  end.

  center_wvl; max_depth; band_area; fwhd
)

NB. Batch continuum removal on 3D cube (B x H x W)
NB. Computes continuum removed cube (B x H x W)
batch_continuum_cube=: 4 : 0
  wvl=. x
  cube=. y
  'B H W'=. $ cube
  out_cube=. (B, H, W) $ 0.0

  for_h. i. H do.
    for_w. i. W do.
      spec=. (<a:; h; w) { cube
      cr_spec=. wvl remove_continuum spec
      out_cube=. cr_spec (<a:; h; w)} out_cube
    end.
  end.
  out_cube
)


albedo_continuum_upper_convex_hull =: upper_convex_hull
albedo_continuum_remove_continuum =: remove_continuum
albedo_continuum_characterize_absorption =: characterize_absorption
