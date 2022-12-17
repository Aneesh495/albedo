NB. ==============================================================================
NB. Test Suite: unmix/vca.ijs, unmix/nfindr.ijs, unmix/ppi.ijs, unmix/nnls.ijs, unmix/fcls.ijs, unmix/sunsal.ijs
NB. ==============================================================================

load 'tests/test_framework.ijs'
load 'unmix/vca.ijs'
load 'unmix/nfindr.ijs'
load 'unmix/ppi.ijs'
load 'unmix/nnls.ijs'
load 'unmix/fcls.ijs'
load 'unmix/sunsal.ijs'
load 'core/linalg.ijs'

cocurrent 'base'

suite_begin 'Endmember Extraction & Spectral Unmixing'

NB. Create synthetic convex simplex with 3 known pure endmembers in 5 bands
M_true=: 5 3 $ 1 0 0  0 1 0  0 0 1  0.5 0.5 0  0 0.5 0.5
NB. 50 pixel mixtures
A_true=: 3 50 $ 0.0
NB. Pixel 0 is pure endmember 0: [1, 0, 0]
A_true=: 1.0 (<0 0)} A_true
NB. Pixel 1 is pure endmember 1: [0, 1, 0]
A_true=: 1.0 (<1 1)} A_true
NB. Pixel 2 is pure endmember 2: [0, 0, 1]
A_true=: 1.0 (<2 2)} A_true
NB. Other pixels are random convex combinations
w_rest =: ? 47 3 $ 100
w_rest =: w_rest % (+/"1 w_rest)
A_rest =: |: w_rest
A_pure =: =/ ~ i. 3
A_true =: A_pure ,. A_rest

X_scene=: M_true albedo_linalg_matmul A_true

NB. Test 1: VCA Endmember Extraction
'M_vca vca_idx'=: 3 albedo_vca_extract_vca X_scene
'VCA extracts 3 endmembers' assert_equal (# vca_idx); 3
'VCA endmember matrix shape is 5 3' assert_shape M_vca; 5 3

NB. Test 2: N-FINDR Simplex Volume Maximization
'M_nf nf_idx max_vol'=: 3 albedo_nfindr_extract_nfindr X_scene
'N-FINDR extracts 3 endmembers' assert_equal (# nf_idx); 3
'N-FINDR volume is positive' assert_true (max_vol > 0.0)

NB. Test 3: Pixel Purity Index (PPI)
'M_ppi ppi_idx hits'=: (3; 100; 95) albedo_ppi_extract_ppi X_scene
'PPI extracts 3 endmembers' assert_equal (# ppi_idx); 3
'PPI pure endmember pixels 0, 1, or 2 have hits' assert_true (+/ (0 1 2 { hits) > 0)

NB. Test 4: NNLS Non-Negative Least Squares
y_test=: 5 $ 0.5
s_nnls=: M_true albedo_nnls_solve_nnls y_test
'NNLS solution length is 3' assert_equal (# s_nnls); 3
'NNLS satisfies non-negativity constraint s >= 0' assert_true (*./ s_nnls >: _1e_12)

NB. Test 5: Fully Constrained Least Squares (FCLS)
s_fcls=: M_true albedo_fcls_solve_fcls (0 {"1 X_scene)
'FCLS for pure endmember 0 has s0 near 1.0' assert_near (0 { s_fcls); 1.0; 1e_4
'FCLS satisfies non-negativity (ANC)' assert_true (*./ s_fcls >: _1e_12)
'FCLS satisfies sum-to-one (ASC)' assert_near (+/ s_fcls); 1.0; 1e_5

NB. Test 6: Batch FCLS Cube Unmixing
CubeSynthetic=: 5 4 4 $ , 16 {."1 X_scene
Abundances=: M_true albedo_fcls_unmix_cube CubeSynthetic
'Abundance cube shape is 3 4 4' assert_shape Abundances; 3 4 4
'All pixel abundance fractions sum to 1.0' assert_near (+/ Abundances); (4 4 $ 1.0); 1e_4
'All pixel abundance fractions >= 0' assert_true (*./ (*./ , Abundances >: _1e_12))

Rmse=: (M_true; Abundances) albedo_fcls_unmix_rmse CubeSynthetic
'Unmixing RMSE is near zero for noiseless data' assert_near Rmse; 0.0; 1e_4

NB. Test 7: SUNSAL Sparse Unmixing
s_sunsal=: (0.01; M_true) albedo_sunsal_solve_sunsal (0 {"1 X_scene)
'SUNSAL satisfies non-negativity' assert_true (*./ s_sunsal >: 0.0)
'SUNSAL endmember 0 dominant' assert_true ((0 { s_sunsal) > (1 { s_sunsal))

suite_end ''
