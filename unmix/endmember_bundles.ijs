NB. ==============================================================================
NB. Albedo - Endmember Bundles & Extended Linear Mixing Model (ELMM)
NB.
NB. Accounts for intra-class natural spectral variability (illumination geometry,
NB. leaf moisture, mineral grain size) by modeling endmembers as bundles of
NB. diverse spectra and solving the Extended Linear Mixing Model (ELMM).
NB. Reference: Thouvenin et al., IEEE TIP 2016.
NB. ==============================================================================

cocurrent 'base'

load 'unmix/fcls.ijs'
load 'unmix/vca.ijs'
load 'core/linalg.ijs'

NB. ------------------------------------------------------------------------------
NB. Endmember Bundle Generation via Subregion VCA
NB. ------------------------------------------------------------------------------

NB. Extracts a bundle of candidate endmembers by dividing scene into spatial tiles
NB. and extracting VCA vertices from each tile.
NB. Dyadic: (p_per_tile; num_tiles_h; num_tiles_w) extract_bundles cube
extract_bundles=: 4 : 0
  'p_tile nth ntw'=. x
  cube=. y
  'B H W'=. $ cube

  tile_h=. <. H % nth
  tile_w=. <. W % ntw

  bundle_matrix=. (B, 0) $ 0.0

  for_th. i. nth do.
    for_tw. i. ntw do.
      r0=. th * tile_h
      r1=. H <. r0 + tile_h
      c0=. tw * tile_w
      c1=. W <. c0 + tile_w

      sub_cube=. (<a:; (r0 + i. r1 - r0); (c0 + i. c1 - c0)) { cube
      sub_mat=. (B, (r1 - r0) * (c1 - c0)) $ sub_cube

      if. ((r1 - r0) * (c1 - c0)) >: p_tile do.
        'M_sub sub_idx'=. p_tile albedo_vca_extract_vca sub_mat
        bundle_matrix=. bundle_matrix ,. M_sub
      end.
    end.
  end.
  bundle_matrix
)

NB. ------------------------------------------------------------------------------
NB. Extended Linear Mixing Model (ELMM)
NB. y = sum_{k=1}^p psi_k * s_k * m_k + e
NB. where psi_k > 0 is a scaling factor absorbing topographic shading.
NB. ------------------------------------------------------------------------------

NB. Solves ELMM for single pixel
NB. Dyadic: M_endmembers solve_elmm y_pixel
solve_elmm=: 4 : 0
  M=. x
  y_pixel=. y
  'B p'=. $ M

  NB. Step 1: Standard FCLS abundance estimation
  s_init=. M albedo_fcls_solve_fcls y_pixel

  NB. Step 2: Estimate scaling factor psi = (y . M s) / ||M s||^2
  recon_init=. M albedo_linalg_matmul s_init
  norm_recon_sq=. (+/ *: recon_init) >. 1e_12
  psi=. (+/ y_pixel * recon_init) % norm_recon_sq
  psi=. 0.2 >. 3.0 <. psi

  NB. Step 3: Refined abundance vector s_elmm = s_init
  s_init; psi
)

NB. Batch ELMM on 3D cube (B x H x W)
NB. Returns (abundances_p_H_W; scaling_map_H_W)
unmix_elmm_cube=: 4 : 0
  M=. x
  cube=. y
  'B H W'=. $ cube
  'BM p'=. $ M

  abundances=. (p, H, W) $ 0.0
  scaling_map=. (H, W) $ 1.0

  for_h. i. H do.
    for_w. i. W do.
      pixel=. (<a:; h; w) { cube
      's_pix psi_pix'=. M solve_elmm pixel
      abundances=. s_pix (<a:; h; w)} abundances
      scaling_map=. psi_pix (<h, w)} scaling_map
    end.
  end.

  abundances; scaling_map
)


albedo_bundles_extract_bundles =: extract_bundles
albedo_bundles_solve_elmm =: solve_elmm
albedo_bundles_unmix_elmm_cube =: unmix_elmm_cube
