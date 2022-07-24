NB. ==============================================================================
NB. Albedo - Quick Atmospheric Correction (QUAC)
NB.
NB. Implements in-scene Quick Atmospheric Correction (Bernstein et al., 2005)
NB. based on universal spectral statistics of diverse terrestrial endmembers.
NB. Retrieves surface reflectance directly from scene data without external
NB. MODTRAN radiative transfer models.
NB. ==============================================================================

cocurrent 'albedo_quac'

load 'core/linalg.ijs'
load 'core/stats.ijs'
load 'unmix/vca.ijs'

NB. ------------------------------------------------------------------------------
NB. QUAC Algorithm
NB. ------------------------------------------------------------------------------

NB. Perform QUAC atmospheric correction on 3D radiance cube (B x H x W)
NB. Dyadic: num_endmembers quac_correction rad_cube
quac_correction=: 4 : 0
  p_endmembers=. x
  rad_cube=. y
  'B H W'=. $ rad_cube
  N=. H * W
  X_mat=. (B, N) $ rad_cube

  NB. Step 1: Find diverse scene endmember pixels via VCA
  'M_vca vca_indices'=. p_endmembers albedo_vca_extract_vca X_mat

  NB. Step 2: Compute mean and standard deviation across endmember spectra
  NB. In QUAC, the standard deviation across endmembers is assumed approximately
  NB. proportional to the atmospheric transmission * solar irradiance
  mu_endmembers=. (+/ "1 M_vca) % p_endmembers
  diffs=. M_vca - mu_endmembers
  sigma_endmembers=. %: (+/ "1 *: diffs) % (p_endmembers - 1)
  sigma_endmembers=. sigma_endmembers >. 1e_6

  NB. Reference universal standard deviation across natural materials (~ 0.15)
  sigma_ref=. 0.15

  NB. Step 3: Gain calibration vector g(lambda) = sigma_ref / sigma_endmembers(lambda)
  gain_quac=. sigma_ref % sigma_endmembers

  NB. Step 4: Baseline offset estimation (approximate path radiance)
  NB. Estimated as the 1st percentile of radiance across scene
  offset_quac=. B $ 0.0
  for_b. i. B do.
    band=. b { rad_cube
    vals=. /:~ , band
    p01=. (<. 0.01 * # vals) { vals
    offset_quac=. (0.0 >. p01) (<b)} offset_quac
  end.

  NB. Step 5: Surface reflectance retrieval: rho = gain * (L - offset)
  refl_cube=. (B, H, W) $ 0.0
  for_b. i. B do.
    g=. b { gain_quac
    o=. b { offset_quac
    band_rad=. b { rad_cube
    refl_b=. g * (0.0 >. band_rad - o)
    refl_cube=. (0.0 >. 1.0 <. refl_b) (<b)} refl_cube
  end.

  refl_cube; gain_quac; offset_quac
)
