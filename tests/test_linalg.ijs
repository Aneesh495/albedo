NB. ==============================================================================
NB. Test Suite: core/linalg.ijs
NB. ==============================================================================

load 'tests/test_framework.ijs'
load 'core/linalg.ijs'

cocurrent 'base'

suite_begin 'Linear Algebra Primitives & Factorizations'

NB. Test 1: Identity, trace, transpose, matmul
I3=: albedo_linalg_eye 3
'I3 trace should be 3' assert_equal (albedo_linalg_trace I3); 3
'I3 shape should be 3 3' assert_shape I3; 3 3

M=: 2 2 $ 1 2 3 4
Mt=: albedo_linalg_tr M
'Transpose of M' assert_equal Mt; (2 2 $ 1 3 2 4)
MMt=: M albedo_linalg_matmul Mt
'MMt shape should be 2 2' assert_shape MMt; 2 2
'MMt top-left should be 5' assert_equal (MMt < 0 0); 5
'MMt bottom-right should be 25' assert_equal (MMt < 1 1); 25

NB. Test 2: Cholesky decomposition: L L^T = A
NB. Construct SPD matrix A = M M^T + I
SPD=: MMt + albedo_linalg_eye 2
L=: albedo_linalg_cholesky SPD
LLt=: L albedo_linalg_matmul (albedo_linalg_tr L)
'Cholesky reconstruction LLt = SPD' assert_near LLt; SPD; 1e_6
'L is lower triangular' assert_equal (L < 0 1); 0.0

NB. Test 3: QR decomposition: Q R = A and Q^T Q = I
A_qr=: 3 2 $ 1 2 3 5 2 4
'Q R'=: albedo_linalg_qr_mgs A_qr
'Q shape should be 3 2' assert_shape Q; 3 2
'R shape should be 2 2' assert_shape R; 2 2
QR_rec=: Q albedo_linalg_matmul R
'QR reconstruction matches A' assert_near QR_rec; A_qr; 1e_6
QtQ=: (albedo_linalg_tr Q) albedo_linalg_matmul Q
'Q columns are orthonormal' assert_near QtQ; (albedo_linalg_eye 2); 1e_6

NB. Test 4: Symmetric Jacobi Eigenvalues: V diag(evals) V^T = A
SymA=: 3 3 $ 4 1 2 1 5 3 2 3 6
'evals V'=: albedo_linalg_eigen_sym SymA
'Eigenvalues should be descending' assert_true (>: / evals)
V_Lambda_Vt=: V albedo_linalg_matmul (albedo_linalg_diag_mat evals) albedo_linalg_matmul (albedo_linalg_tr V)
'Eigenvalue reconstruction matches SymA' assert_near V_Lambda_Vt; SymA; 1e_5
VtV=: (albedo_linalg_tr V) albedo_linalg_matmul V
'Eigenvectors are orthonormal' assert_near VtV; (albedo_linalg_eye 3); 1e_5

NB. Test 5: SVD: U Sigma V^T = A
A_svd=: 4 3 $ 1 2 3 4 5 6 7 8 10 2 1 4
'svals U V_svd'=: albedo_linalg_svd_power A_svd
'Singular values should be descending' assert_true (>: / svals)
USigmaVt=: U albedo_linalg_matmul (albedo_linalg_diag_mat svals) albedo_linalg_matmul (albedo_linalg_tr V_svd)
'SVD reconstruction matches A_svd' assert_near USigmaVt; A_svd; 1e_4

NB. Test 6: Moore-Penrose Pseudoinverse: A A^+ A = A
A_pinv=: albedo_linalg_pinv_svd A_svd
A_pinv_A=: A_svd albedo_linalg_matmul A_pinv albedo_linalg_matmul A_svd
'Moore-Penrose property A A^+ A = A' assert_near A_pinv_A; A_svd; 1e_4

NB. Test 7: Orthogonal Projectors
P=: albedo_linalg_proj_mat A_svd
P2=: P albedo_linalg_matmul P
'Projector is idempotent P^2 = P' assert_near P2; P; 1e_5
P_perp=: albedo_linalg_proj_perp_mat A_svd
P_perp2=: P_perp albedo_linalg_matmul P_perp
'Nullspace projector is idempotent P_perp^2 = P_perp' assert_near P_perp2; P_perp; 1e_5

suite_end ''
