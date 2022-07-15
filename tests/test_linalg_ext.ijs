NB. ==============================================================================
NB. Test Suite: core/linalg_ext.ijs
NB. ==============================================================================

load 'tests/test_framework.ijs'
load 'core/linalg_ext.ijs'
load 'core/linalg.ijs'

cocurrent 'base'

suite_begin 'Advanced Linear Algebra Extensions (PCG, Whitening, GSVD)'

NB. Construct positive definite symmetric matrix A
M=: 3 3 $ 3 1 0 1 4 2 0 2 5
SPD=: (albedo_linalg_tr M) albedo_linalg_matmul M

NB. Test 1: Symmetric Matrix Square Root A^(1/2)
SqrtA=: albedo_linalg_ext_sqrt_mat_sym SPD
SqrtA_sq=: SqrtA albedo_linalg_matmul SqrtA
'Matrix square root squared equals SPD' assert_near SqrtA_sq; SPD; 1e_5

NB. Test 2: Inverse Matrix Square Root A^(-1/2) (Whitening Transform)
InvSqrtA=: albedo_linalg_ext_inv_sqrt_mat_sym SPD
Whitened=: InvSqrtA albedo_linalg_matmul SPD albedo_linalg_matmul InvSqrtA
'A^(-1/2) A A^(-1/2) equals Identity' assert_near Whitened; (albedo_linalg_eye 3); 1e_5

NB. Test 3: Log-Determinant and Determinant
DiagMat=: albedo_linalg_diag_mat 2.0 3.0 4.0
LogDet=: albedo_linalg_ext_log_det_spd DiagMat
'Log det of diag(2,3,4) is log(24)' assert_near LogDet; (^. 24.0); 1e_6
DetVal=: albedo_linalg_ext_det_spd DiagMat
'Determinant of diag(2,3,4) is 24' assert_near DetVal; 24.0; 1e_5

NB. Test 4: Preconditioned Conjugate Gradient (PCG)
b_vec=: 3 $ 7.0 11.0 13.0
x_pcg=: (50; 1e_8) albedo_linalg_ext_pcg_solve (SPD; b_vec)
Ax_pcg=: SPD albedo_linalg_matmul x_pcg
'PCG solution satisfies A x = b' assert_near Ax_pcg; b_vec; 1e_5

NB. Test 5: Generalized SVD (GSVD) Pair
B_mat=: albedo_linalg_eye 3
'g_evals g_evecs'=: SPD albedo_linalg_ext_gsvd_pair B_mat
'GSVD eigenvalues count is 3' assert_equal (# g_evals); 3
'GSVD eigenvalues are descending' assert_true (>: / g_evals)

suite_end ''
