NB. ==============================================================================
NB. Test Suite: core/tensor.ijs
NB. ==============================================================================

load 'tests/test_framework.ijs'
load 'core/tensor.ijs'

cocurrent 'base'

suite_begin 'Tensor Algebra Primitives'

NB. Test 1: Tensor basic order and size
T1=: i. 2 3 4
'T1 order should be 3' assert_equal (albedo_tensor_tensor_order T1); 3
'T1 size should be 24' assert_equal (albedo_tensor_tensor_size T1); 24

NB. Test 2: Frobenius norm
V=: 3 4
'V norm should be 5' assert_near (albedo_tensor_tensor_fnorm V); 5.0; 1e_6
M=: 2 2 $ 1 2 2 4
'M norm should be 5' assert_near (albedo_tensor_tensor_fnorm M); 5.0; 1e_6

NB. Test 3: Kronecker product
A=: 2 2 $ 1 2 3 4
B=: 2 2 $ 0 5 6 7
K=: A albedo_tensor_kronecker B
'Kronecker shape should be 4 4' assert_shape K; 4 4
'Kronecker entry (0,1) should be 5' assert_equal (K < 0 1); 5
'Kronecker entry (2,2) should be 0' assert_equal (K < 2 2); 0
'Kronecker entry (3,3) should be 28' assert_equal (K < 3 3); 28

NB. Test 4: Khatri-Rao product
A_kr=: 2 2 $ 1 2 3 4
B_kr=: 2 2 $ 5 6 7 8
KR=: A_kr albedo_tensor_khatri_rao B_kr
'Khatri-Rao shape should be 4 2' assert_shape KR; 4 2
NB. Column 0: [1, 3] (x) [5, 7] = [5, 7, 15, 21]
'KR column 0 should match' assert_equal (0 {"1 KR); 5 7 15 21
NB. Column 1: [2, 4] (x) [6, 8] = [12, 16, 24, 32]
'KR column 1 should match' assert_equal (1 {"1 KR); 12 16 24 32

NB. Test 5: Mode-n Unfold and Fold round-trip
X3D=: i. 3 4 5
NB. Mode 0 unfold
U0=: 0 albedo_tensor_tensor_unfold X3D
'Mode 0 unfold shape should be 3 20' assert_shape U0; 3 20
R0=: (0; 3 4 5) albedo_tensor_tensor_fold U0
'Mode 0 fold roundtrip match' assert_equal R0; X3D

NB. Mode 1 unfold
U1=: 1 albedo_tensor_tensor_unfold X3D
'Mode 1 unfold shape should be 4 15' assert_shape U1; 4 15
R1=: (1; 3 4 5) albedo_tensor_tensor_fold U1
'Mode 1 fold roundtrip match' assert_equal R1; X3D

NB. Mode 2 unfold
U2=: 2 albedo_tensor_tensor_unfold X3D
'Mode 2 unfold shape should be 5 12' assert_shape U2; 5 12
R2=: (2; 3 4 5) albedo_tensor_tensor_fold U2
'Mode 2 fold roundtrip match' assert_equal R2; X3D

NB. Test 6: Mode-n Matrix-Tensor product
NB. Multiply X (3 x 4 x 5) along mode 0 by U (2 x 3) -> result (2 x 4 x 5)
U_mat=: 2 3 $ 1 0 0 0 1 0
Y_prod=: (0; U_mat) albedo_tensor_mode_n_prod X3D
'Mode 0 product shape should be 2 4 5' assert_shape Y_prod; 2 4 5
'Mode 0 product slice 0 should equal X3D slice 0' assert_equal (0 { Y_prod); (0 { X3D)
'Mode 0 product slice 1 should equal X3D slice 1' assert_equal (1 { Y_prod); (1 { X3D)

NB. Test 7: Cube to spectral matrix and back
Cube=: ? 10 8 8 $ 100
SpecMat=: albedo_tensor_cube_to_spectral_matrix Cube
'Spectral matrix shape should be 10 64' assert_shape SpecMat; 10 64
Cube_rec=: 8 8 albedo_tensor_spectral_matrix_to_cube SpecMat
'Cube reconstruction should match exact' assert_equal Cube_rec; Cube

suite_end ''
