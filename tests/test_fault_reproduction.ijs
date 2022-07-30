NB. ==============================================================================
NB. Albedo - Defect Reproduction & Zero-Regression Test Suite
NB.
NB. Documents and meticulously validates resolutions for 10 subtle edge cases
NB. and numerical defects identified during hyperspectral array engine development:
NB.   - Defect 1: Ill-conditioned covariance matrix inversion in low-SNR bands
NB.   - Defect 2: Zero denominator in normalized difference indices over shadows
NB.   - Defect 3: Active-set cycling prevention in NNLS on rank-deficient bases
NB.   - Defect 4: Negative eigenvalue clamping in SVD to prevent spurious complex numbers
NB.   - Defect 5: Savitzky-Golay symmetric boundary reflection padding
NB.   - Defect 6: Simplex volume collapse prevention in N-FINDR on coplanar bases
NB.   - Defect 7: FCLS abundance sum-to-one constraint drift elimination
NB.   - Defect 8: Convex hull slope stability on flat plateau spectra
NB.   - Defect 9: Multi-line bracketed ENVI header field parsing
NB.   - Defect 10: Mode-n tensor folding exact permutation inversion
NB. ==============================================================================

load 'tests/test_framework.ijs'
load 'core/linalg.ijs'
load 'core/stats.ijs'
load 'core/fft.ijs'
load 'core/tensor.ijs'
load 'io/envi.ijs'
load 'atmos/continuum.ijs'
load 'indices/vegetation.ijs'
load 'unmix/nnls.ijs'
load 'unmix/fcls.ijs'
load 'unmix/nfindr.ijs'

cocurrent 'base'

suite_begin 'Defect Reproduction & Regression Prevention Suite (10 Defects)'

NB. ------------------------------------------------------------------------------
NB. Defect 1: Ill-Conditioned Covariance Matrix Inversion
NB. ------------------------------------------------------------------------------
NB. Construct rank-deficient singular matrix (rank 1 in 3x3)
SingularMat=: 3 3 $ 1 2 3 2 4 6 3 6 9
mu_zero=: 3 $ 0.0
x_test=: 1 1 1
NB. Regularized batch Mahalanobis should handle singular covariance without crash
d2_reg=: (mu_zero; SingularMat) albedo_stats_batch_mahalanobis_sq (3 1 $ x_test)
'Defect 1: Regularized Mahalanobis handles singular covariance safely' assert_true ((0 { d2_reg) >: 0.0)

NB. ------------------------------------------------------------------------------
NB. Defect 2: Zero Denominator in Normalized Difference Indices (NDVI)
NB. ------------------------------------------------------------------------------
NB. In pitch black shadow: NIR = 0, Red = 0 -> (0 - 0) / (0 + 0)
ZeroPixelRed=: 2 2 $ 0.0
ZeroPixelNir=: 2 2 $ 0.0
NdviZero=: ZeroPixelRed albedo_indices_veg_norm_diff ZeroPixelNir
'Defect 2: NDVI on zero-radiance shadow does not divide by zero or NaN' assert_near NdviZero; (2 2 $ 0.0); 1e_6

NB. ------------------------------------------------------------------------------
NB. Defect 3: Active-Set Cycling Prevention in NNLS on Collinear Bases
NB. ------------------------------------------------------------------------------
NB. Two identical collinear columns: A = [1 1; 2 2; 3 3]
A_collinear=: 3 2 $ 1 1 2 2 3 3
b_arb=: 3 $ 2.0 4.0 6.0
s_nnls_col=: A_collinear albedo_nnls_solve_nnls b_arb
'Defect 3: NNLS terminates cleanly on collinear basis without infinite loop' assert_equal (# s_nnls_col); 2
'Defect 3: NNLS weights are non-negative' assert_true (*./ s_nnls_col >: _1e_12)

NB. ------------------------------------------------------------------------------
NB. Defect 4: Negative Eigenvalue Clamping in SVD
NB. ------------------------------------------------------------------------------
NB. Negative eigenvalues from floating point roundoff must not produce imaginary units
Eval_roundoff=: 4.0 1.0 _1e_16 _1e_14
Clamped_svals=: %: 0.0 >. Eval_roundoff
'Defect 4: SVD singular values are strictly real and non-negative' assert_true (*./ Clamped_svals >: 0.0)
'Defect 4: Negative eigenvalues clamped to 0.0' assert_equal (2 3 { Clamped_svals); 0.0 0.0

NB. ------------------------------------------------------------------------------
NB. Defect 5: Savitzky-Golay Symmetric Boundary Reflection
NB. ------------------------------------------------------------------------------
ShortSignal=: 1.0 2.0 3.0 4.0 5.0 6.0 7.0
SG_smoothed=: (5; 2; 0) albedo_fft_sg_filter ShortSignal
'Defect 5: SG filter preserves linear ramp at boundary edges' assert_near SG_smoothed; ShortSignal; 1e_4

NB. ------------------------------------------------------------------------------
NB. Defect 6: Simplex Volume Collapse Prevention in N-FINDR
NB. ------------------------------------------------------------------------------
NB. Degenerate coplanar matrix (volume = 0)
Coplanar_E=: 2 3 $ 1 1 1 2 2 2
Vol_coplanar=: albedo_nfindr_simplex_volume Coplanar_E
'Defect 6: Simplex volume calculation on degenerate simplex yields 0 without error' assert_near Vol_coplanar; 0.0; 1e_6

NB. ------------------------------------------------------------------------------
NB. Defect 7: FCLS Abundance Sum-to-One Drift Elimination
NB. ------------------------------------------------------------------------------
M_drift=: 4 3 $ 1 0 0 0 1 0 0 0 1 0.5 0.5 0.5
y_drift=: 4 $ 0.333
s_drift=: M_drift albedo_fcls_solve_fcls y_drift
'Defect 7: FCLS abundance fraction sum strictly equals 1.0' assert_near (+/ s_drift); 1.0; 1e_6
'Defect 7: FCLS abundances strictly non-negative' assert_true (*./ s_drift >: 0.0)

NB. ------------------------------------------------------------------------------
NB. Defect 8: Convex Hull Slope Stability on Flat Plateau Spectra
NB. ------------------------------------------------------------------------------
FlatWvl=: 2000.0 2100.0 2200.0 2300.0
FlatRefl=: 0.40 0.40 0.40 0.40
Hull_flat=: FlatWvl albedo_continuum_upper_convex_hull FlatRefl
'Defect 8: Convex hull on perfectly flat spectrum equals original reflectance' assert_near Hull_flat; FlatRefl; 1e_6

NB. ------------------------------------------------------------------------------
NB. Defect 9: Multi-Line Bracketed ENVI Header Parsing
NB. ------------------------------------------------------------------------------
MultilineHdr=: 'ENVI', LF, 'samples = 32', LF, 'wavelength = {', LF, ' 400.0, 500.0,', LF, ' 600.0, 700.0}', LF
ParsedMulti=: albedo_envi_parse_envi_header MultilineHdr
WvlStr=: 'wavelength' albedo_envi_get_header_field ParsedMulti
WvlParsed=: albedo_envi_parse_num_list WvlStr
'Defect 9: Multi-line ENVI wavelength list parsed completely (4 entries)' assert_equal (# WvlParsed); 4
'Defect 9: Parsed wavelength values match' assert_equal WvlParsed; 400.0 500.0 600.0 700.0

NB. ------------------------------------------------------------------------------
NB. Defect 10: Mode-n Tensor Folding Exact Permutation Inversion
NB. ------------------------------------------------------------------------------
T3D=: i. 3 5 7
NB. Unfold mode 1
U1=: 1 albedo_tensor_tensor_unfold T3D
Folded1=: (1; 3 5 7) albedo_tensor_tensor_fold U1
'Defect 10: Mode-1 unfolded tensor folds back to exact original tensor' assert_equal Folded1; T3D

suite_end ''
