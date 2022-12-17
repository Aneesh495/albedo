NB. ==============================================================================
NB. ALBEDO - High-Dimensional Hyperspectral Imaging Engine
NB. data/polymers_plastics.ijs : Synthetic Polymers and Plastics Library (20 Classes)
NB. Spectral reflectance signatures for urban materials, polymers, and marine debris
NB. plastics from 0.400 to 2.500 micrometers with diagnostic C-H vibrational overtones.
NB. ==============================================================================

cocurrent 'base'

NB. 224 Wavelength centers in micrometers (0.400 to 2.500 um)
WAVELENGTHS =: 0.400 + 0.009375 * i. 224

NB. 20 Polymer and Plastic Material Names
PLASTIC_NAMES =: <;._1 ' PET' , ' HDPE' , ' PVC' , ' LDPE' , ' PP' , ' PS' , ' PMMA' , ' PTFE' , ' Polyurethane' , ' Polycarbonate' , ' Nylon6' , ' Nylon66' , ' ABS' , ' PLA' , ' Epoxy' , ' PEEK' , ' PVA' , ' CelluloseAcetate' , ' MaritimeRopePE' , ' EPSFoam'

NB. Helper to construct characteristic polymer reflectance curve:
NB. r0: base reflectance, slope: tilt, dips: matrix of (center_um, depth, width_um)
calc_plastic =: 3 : 0
  'r0 slope dips' =. y
  base =. r0 + slope * (WAVELENGTHS - 0.400)
  for_d. dips do.
    'c0 depth w0' =. d
    gauss =. depth * ^. - 0.5 * *: (WAVELENGTHS - c0) % w0
    base =. base - gauss
  end.
  0.01 >. 0.98 <. base
)

NB. Diagnostic C-H, N-H, O-H, and C=O overtones and combination bands:
NB. Polyethylene (HDPE/LDPE): 1.21 um, 1.41 um, 1.73 um, 2.31 um, 2.35 um
pe_dips   =: 5 3 $ 1.210 0.12 0.025  1.410 0.08 0.030  1.730 0.35 0.035  2.310 0.38 0.035  2.350 0.32 0.030
pp_dips   =: 5 3 $ 1.195 0.10 0.025  1.395 0.07 0.030  1.705 0.32 0.035  2.315 0.35 0.035  2.380 0.28 0.030
pet_dips  =: 4 3 $ 1.130 0.06 0.025  1.660 0.22 0.030  2.130 0.24 0.035  2.260 0.18 0.030
pvc_dips  =: 4 3 $ 1.200 0.08 0.025  1.720 0.28 0.035  2.160 0.15 0.030  2.300 0.22 0.035
ps_dips   =: 5 3 $ 1.145 0.08 0.025  1.680 0.36 0.030  2.160 0.30 0.035  2.300 0.20 0.030  2.410 0.16 0.030
pmma_dips =: 4 3 $ 1.185 0.09 0.025  1.400 0.08 0.030  1.720 0.25 0.035  2.150 0.20 0.035
ptfe_dips =: 3 3 $ 1.250 0.04 0.030  1.750 0.06 0.035  2.200 0.05 0.035
pu_dips   =: 4 3 $ 1.210 0.08 0.025  1.510 0.12 0.035  1.730 0.26 0.035  2.220 0.20 0.035
pc_dips   =: 4 3 $ 1.180 0.07 0.025  1.670 0.28 0.030  2.150 0.22 0.035  2.280 0.16 0.030
ny6_dips  =: 4 3 $ 1.200 0.08 0.025  1.530 0.16 0.035  1.720 0.28 0.035  2.300 0.22 0.035
ny66_dips =: 4 3 $ 1.205 0.08 0.025  1.535 0.17 0.035  1.725 0.29 0.035  2.305 0.23 0.035
abs_dips  =: 4 3 $ 1.150 0.07 0.025  1.685 0.30 0.030  2.165 0.24 0.035  2.310 0.20 0.035
pla_dips  =: 4 3 $ 1.200 0.07 0.025  1.450 0.10 0.030  1.720 0.24 0.035  2.280 0.18 0.035
epox_dips =: 4 3 $ 1.180 0.06 0.025  1.670 0.24 0.030  2.140 0.20 0.035  2.280 0.18 0.035
peek_dips =: 3 3 $ 1.180 0.05 0.025  1.670 0.20 0.030  2.140 0.18 0.035
pva_dips  =: 4 3 $ 1.210 0.07 0.025  1.440 0.18 0.035  1.730 0.22 0.035  2.100 0.16 0.035
ca_dips   =: 4 3 $ 1.200 0.06 0.025  1.440 0.14 0.030  1.720 0.22 0.035  2.280 0.17 0.035
rope_dips =: 5 3 $ 1.210 0.14 0.025  1.410 0.10 0.030  1.730 0.38 0.035  2.310 0.40 0.035  2.350 0.34 0.030
eps_dips  =: 5 3 $ 1.145 0.09 0.025  1.680 0.40 0.030  2.160 0.34 0.035  2.300 0.22 0.030  2.410 0.18 0.030

p_pet  =: calc_plastic 0.70 _0.08 pet_dips
p_hdpe =: calc_plastic 0.75 _0.06 pe_dips
p_pvc  =: calc_plastic 0.65 _0.05 pvc_dips
p_ldpe =: calc_plastic 0.72 _0.06 pe_dips
p_pp   =: calc_plastic 0.74 _0.07 pp_dips
p_ps   =: calc_plastic 0.78 _0.08 ps_dips
p_pmma =: calc_plastic 0.82 _0.09 pmma_dips
p_ptfe =: calc_plastic 0.88 _0.04 ptfe_dips
p_pu   =: calc_plastic 0.60 _0.05 pu_dips
p_pc   =: calc_plastic 0.76 _0.08 pc_dips
p_ny6  =: calc_plastic 0.68 _0.06 ny6_dips
p_ny66 =: calc_plastic 0.70 _0.06 ny66_dips
p_abs  =: calc_plastic 0.64 _0.05 abs_dips
p_pla  =: calc_plastic 0.72 _0.07 pla_dips
p_epox =: calc_plastic 0.58 _0.04 epox_dips
p_peek =: calc_plastic 0.52 _0.04 peek_dips
p_pva  =: calc_plastic 0.66 _0.06 pva_dips
p_ca   =: calc_plastic 0.74 _0.07 ca_dips
p_rope =: calc_plastic 0.70 _0.06 rope_dips
p_eps  =: calc_plastic 0.84 _0.08 eps_dips

NB. Shape: 20 rows by 224 columns
PLASTIC_MATRIX =: > p_pet ; p_hdpe ; p_pvc ; p_ldpe ; p_pp ; p_ps ; p_pmma ; p_ptfe ; p_pu ; p_pc ; p_ny6 ; p_ny66 ; p_abs ; p_pla ; p_epox ; p_peek ; p_pva ; p_ca ; p_rope ; p_eps

NB. Query plastic index by exact name
get_plastic_index =: 3 : 0
  idx =. PLASTIC_NAMES i. < y
  if. idx = # PLASTIC_NAMES do. _1 else. idx end.
)

NB. Get spectral signature by plastic name
get_plastic_spectrum =: 3 : 0
  idx =. get_plastic_index y
  if. idx = _1 do.
    224 # 0
  else.
    idx { PLASTIC_MATRIX
  end.
)

NB. Floating Plastic Index (FPI) or Hydrocarbon Detection Index
NB. Uses 1.21 um and 1.73 um absorption overtones relative to continuum
plastic_hydrocarbon_index =: 3 : 0
  spec =. y
  NB. Baseline continuum near 1.60 um and 1.80 um
  i160 =. 128
  i180 =. 149
  i173 =. 142
  cont =. 0.5 * (i160 { spec) + (i180 { spec)
  depth =. cont - (i173 { spec)
  depth % (cont + 1e_6)
)


NB. Exported Aliases
albedo_plastics_PLASTIC_NAMES =: PLASTIC_NAMES
albedo_plastics_WAVELENGTHS =: WAVELENGTHS
albedo_plastics_PLASTIC_MATRIX =: PLASTIC_MATRIX
albedo_plastics_get_plastic_index =: get_plastic_index
albedo_plastics_get_plastic_spectrum =: get_plastic_spectrum
albedo_plastics_plastic_hydrocarbon_index =: plastic_hydrocarbon_index
