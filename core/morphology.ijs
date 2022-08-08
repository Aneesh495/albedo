NB. ==============================================================================
NB. Albedo - Mathematical Morphology & Morphological Profiles
NB.
NB. Implements 2D grayscale mathematical morphology operators:
NB.   - Dilation & Erosion with flat structuring elements
NB.   - Opening & Closing
NB.   - White Top-Hat & Black Top-Hat transforms
NB.   - Multiscale Morphological Profiles (MP) for spatial-spectral feature extraction
NB. ==============================================================================

cocurrent 'albedo_morphology'

NB. ------------------------------------------------------------------------------
NB. Grayscale Dilation & Erosion with Flat Box Structuring Element (size x size)
NB. ------------------------------------------------------------------------------

NB. Flat box dilation of 2D matrix M
NB. Dyadic: radius dilate_2d M
dilate_2d=: 4 : 0
  r=. x
  M=. y
  'H W'=. $ M
  out=. M

  for_h. i. H do.
    for_w. i. W do.
      r_min=. 0 >. h - r
      r_max=. (<: H) <. h + r
      c_min=. 0 >. w - r
      c_max=. (<: W) <. w + r

      rows=. r_min + i. 1 + r_max - r_min
      cols=. c_min + i. 1 + c_max - c_min
      window=. , M < rows; cols
      out=. (>./ window) (<h, w)} out
    end.
  end.
  out
)

NB. Flat box erosion of 2D matrix M
NB. Dyadic: radius erode_2d M
erode_2d=: 4 : 0
  r=. x
  M=. y
  'H W'=. $ M
  out=. M

  for_h. i. H do.
    for_w. i. W do.
      r_min=. 0 >. h - r
      r_max=. (<: H) <. h + r
      c_min=. 0 >. w - r
      c_max=. (<: W) <. w + r

      rows=. r_min + i. 1 + r_max - r_min
      cols=. c_min + i. 1 + c_max - c_min
      window=. , M < rows; cols
      out=. (<./ window) (<h, w)} out
    end.
  end.
  out
)

NB. ------------------------------------------------------------------------------
NB. Opening, Closing & Top-Hat Transforms
NB. ------------------------------------------------------------------------------

NB. Morphological Opening: gamma = dilate(erode(M))
NB. Removes bright features smaller than structuring element
open_2d=: 4 : 0
  r=. x
  M=. y
  r dilate_2d (r erode_2d M)
)

NB. Morphological Closing: phi = erode(dilate(M))
NB. Fills dark features/holes smaller than structuring element
close_2d=: 4 : 0
  r=. x
  M=. y
  r erode_2d (r dilate_2d M)
)

NB. White Top-Hat (WTH): M - open(M)
NB. Isolates bright objects/anomalies smaller than structuring element
top_hat_white=: 4 : 0
  r=. x
  M=. y
  0.0 >. M - (r open_2d M)
)

NB. Black Top-Hat (BTH): close(M) - M
NB. Isolates dark objects/troughs smaller than structuring element
top_hat_black=: 4 : 0
  r=. x
  M=. y
  0.0 >. (r close_2d M) - M
)

NB. ------------------------------------------------------------------------------
NB. Multiscale Morphological Profile (MP)
NB. Builds stack of openings and closings across multiple radii (e.g. 1, 2, 3, 4)
NB. ------------------------------------------------------------------------------

NB. Build morphological profile for a 2D band
NB. Dyadic: radii_list build_mp_2d M
build_mp_2d=: 4 : 0
  radii=. x
  M=. y
  'H W'=. $ M
  n_r=. # radii

  NB. Profile contains: openings by reconstruction + original + closings
  NB. Total bands: 2 * n_r + 1
  tot_bands=. (2 * n_r) + 1
  MP=. (tot_bands, H, W) $ 0.0

  NB. Opening profiles
  for_k. i. n_r do.
    r=. k { radii
    op=. r open_2d M
    MP=. op (<k)} MP
  end.

  NB. Original image at center
  MP=. M (<n_r)} MP

  NB. Closing profiles
  for_k. i. n_r do.
    r=. k { radii
    cl=. r close_2d M
    MP=. cl (<(n_r + 1 + k))} MP
  end.

  MP
)
