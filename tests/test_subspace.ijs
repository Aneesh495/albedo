NB. ==============================================================================
NB. Test Suite: subspace/pca.ijs, subspace/mnf.ijs, subspace/hysime.ijs, subspace/segmented_pca.ijs
NB. ==============================================================================

load 'tests/test_framework.ijs'
load 'subspace/pca.ijs'
load 'subspace/mnf.ijs'
load 'subspace/hysime.ijs'
load 'subspace/segmented_pca.ijs'
load 'io/spectral_lib.ijs'

cocurrent 'base'

suite_begin 'Subspace Projection & Dimensionality Reduction'

NB. Construct synthetic low-rank data: 3 pure endmembers mixed across 100 pixels
Wvl20=: 400.0 + (i. 20) * 100.0
B_count=: 20
N_pixels=: 100

NB. 3 random endmembers (B x 3)
E_true=: ? B_count, 3 $ 100
NB. 3 abundance fractions (3 x N) summing to 1
A_raw=: ? 3, N_pixels $ 100
A_true=: A_raw % (+/ "1 A_raw)
NB. Linear mixture X = E A + small noise
X_clean=: E_true albedo_linalg_matmul A_true
X_noisy=: X_clean + (0.01 * ? B_count, N_pixels $ 10)

NB. Test 1: PCA Model Training & Variance
PcaModel=: albedo_pca_train_pca X_noisy
'evals evecs'=: 1 2 { PcaModel
'PCA eigenvalues descending' assert_true (>: / evals)
CumVar=: 4 {:: PcaModel
'Cumulative variance reaches 1.0' assert_near ({: CumVar); 1.0; 1e_5

NB. Test 2: PCA Projection and Reconstruction
NB. Project to k=3 components
Y3=: (3; PcaModel) albedo_pca_project_pca X_noisy
'Y3 shape is 3 100' assert_shape Y3; 3 N_pixels
X_rec3=: (3; PcaModel) albedo_pca_reconstruct_pca Y3
'Reconstructed shape is 20 100' assert_shape X_rec3; B_count N_pixels

Mse3=: X_noisy albedo_pca_reconstruction_mse X_rec3
Mse1=: X_noisy albedo_pca_reconstruction_mse ((1; PcaModel) albedo_pca_reconstruct_pca ((1; PcaModel) albedo_pca_project_pca X_noisy))
'MSE with 3 components is smaller than with 1 component' assert_true (Mse3 < Mse1)

NB. Test 3: MNF Transform on 3D Cube
CubeSynthetic=: 10 10 albedo_tensor_spectral_matrix_to_cube X_noisy
MnfModel=: albedo_mnf_train_mnf CubeSynthetic
'MnfModel elements' assert_equal (# MnfModel); 5
'snr_evals'=: 4 {:: MnfModel
'MNF SNR eigenvalues are descending' assert_true (>: / snr_evals)

NB. Forward MNF project
Y_mnf=: (3; MnfModel) albedo_mnf_project_mnf X_noisy
'MNF Y shape is 3 100' assert_shape Y_mnf; 3 N_pixels

NB. Test 4: HySime Subspace Dimension Estimation
'dim_k costs U_hy'=: albedo_hysime_hysime_subspace_dim X_noisy
'HySime estimates intrinsic dimension near 3' assert_true ((dim_k >: 1) *. (dim_k <: 6))

NB. Test 5: Segmented PCA
W_full=: albedo_spectral_lib_standard_wavelengths_224 ''
'vnir swir1 swir2'=: albedo_segmented_pca_partition_bands W_full
'VNIR bands non-empty' assert_true (0 < # vnir)
'SWIR1 bands non-empty' assert_true (0 < # swir1)
'SWIR2 bands non-empty' assert_true (0 < # swir2)

NB. Train segmented PCA on 224-band matrix with 50 pixels
X224=: ? 224 50 $ 100
K_parts=: 2 2 2
SegModel=: (K_parts; W_full) albedo_segmented_pca_train_segmented_pca X224
ProjSeg=: SegModel albedo_segmented_pca_project_segmented_pca X224
'Segmented projection shape is 6 50 (2+2+2 components)' assert_shape ProjSeg; 6 50

suite_end ''
