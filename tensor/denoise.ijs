NB. ==============================================================================
NB. Albedo - Low-Rank Tensor Recovery & Cube Denoising
NB.
NB. Strips random noise, spatial striping artifacts, and sensor readout defects
NB. from 3D hyperspectral data cubes using low-rank Tucker and CP approximations.
NB. ==============================================================================

cocurrent 'base'

load 'tensor/tucker.ijs'
load 'tensor/cp_als.ijs'
load 'core/tensor.ijs'

NB. ------------------------------------------------------------------------------
NB. Low-Rank Tucker Denoising
NB. ------------------------------------------------------------------------------

NB. Denoise 3D cube using truncated HOSVD
NB. Dyadic: (r_spectral; r_spatial_h; r_spatial_w) denoise_tucker cube
denoise_tucker=: 4 : 0
  'r_b r_h r_w'=. x
  cube=. y
  'B H W'=. $ cube
  r_b=. r_b <. B
  r_h=. r_h <. H
  r_w=. r_w <. W

  tucker_model=. (r_b; r_h; r_w) albedo_tucker_fit_hosvd cube
  albedo_tucker_reconstruct_tucker tucker_model
)

NB. ------------------------------------------------------------------------------
NB. Low-Rank CP Denoising
NB. ------------------------------------------------------------------------------

NB. Denoise 3D cube using rank-R CP decomposition
NB. Dyadic: rank_R denoise_cp cube
denoise_cp=: 4 : 0
  rank_R=. x
  cube=. y
  cp_model=. (rank_R; 15; 1e_4) albedo_cp_als_fit_cp_als cube
  albedo_cp_als_reconstruct_cp cp_model
)


albedo_denoise_denoise_tucker =: denoise_tucker
