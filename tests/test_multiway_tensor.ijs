NB. ==============================================================================
NB. Test Suite: tensor/cp_als.ijs, tensor/tucker.ijs, tensor/denoise.ijs
NB. ==============================================================================

load 'tests/test_framework.ijs'
load 'tensor/cp_als.ijs'
load 'tensor/tucker.ijs'
load 'tensor/denoise.ijs'
load 'core/tensor.ijs'

cocurrent 'base'

suite_begin 'Multi-Way Tensor Decompositions (CP-ALS & Tucker3 HOSVD)'

NB. Test 1: Rank-1 CP Tensor Decomposition
a_true=: 1.0 2.0 3.0
b_true=: 4.0 5.0
c_true=: 6.0 7.0 8.0 9.0
T_rank1=: (a_true; b_true; c_true) albedo_tensor_vec_outer3 ''
'T_rank1 shape is 3 2 4' assert_shape T_rank1; 3 2 4

Cp1=: (1; 10; 1e_6) albedo_cp_als_fit_cp_als T_rank1
'Cp1 elements' assert_equal (# Cp1); 4
T_rec1=: albedo_cp_als_reconstruct_cp Cp1
'Rank-1 CP reconstruction matches original' assert_near T_rec1; T_rank1; 1e_4

NB. Test 2: Rank-2 CP Tensor Decomposition
T_rank2=: T_rank1 + ((3 2 4 $ ? 24 $ 100) * 0.1)
Cp2=: (2; 15; 1e_4) albedo_cp_als_fit_cp_als T_rank2
T_rec2=: albedo_cp_als_reconstruct_cp Cp2
'Rank-2 CP reconstruction shape is 3 2 4' assert_shape T_rec2; 3 2 4
Err2=: albedo_tensor_tensor_fnorm (T_rank2 - T_rec2)
'Rank-2 CP reconstruction error is low' assert_true (Err2 < 10.0)

NB. Test 3: Tucker3 / HOSVD Full-Rank Exact Recovery
NB. When ranks equal full dimensions (3, 2, 4), HOSVD is an exact orthogonal basis transformation
Hosvd_full=: (3; 2; 4) albedo_tucker_fit_hosvd T_rank1
'G_core U0 U1 U2'=: Hosvd_full
'Core tensor shape is 3 2 4' assert_shape G_core; 3 2 4
'U0 shape is 3 3' assert_shape U0; 3 3
'U1 shape is 2 2' assert_shape U1; 2 2
'U2 shape is 4 4' assert_shape U2; 4 4

T_hosvd_rec=: albedo_tucker_reconstruct_tucker Hosvd_full
'Full-rank Tucker reconstruction matches T_rank1' assert_near T_hosvd_rec; T_rank1; 1e_5

NB. Test 4: Truncated Tucker Low-Rank Approximation
NB. Truncate to (2, 2, 2)
Hosvd_trunc=: (2; 2; 2) albedo_tucker_fit_hosvd T_rank1
'Truncated core tensor shape is 2 2 2' assert_shape (0 {:: Hosvd_trunc); 2 2 2
T_trunc_rec=: albedo_tucker_reconstruct_tucker Hosvd_trunc
'Truncated reconstruction shape is 3 2 4' assert_shape T_trunc_rec; 3 2 4

NB. Test 5: Low-Rank Tensor Denoising
Noise=: (3 2 4 $ ? 24 $ 10) * 0.05
T_noisy=: T_rank1 + Noise
T_denoised=: (1; 1; 1) albedo_denoise_denoise_tucker T_noisy
'Denoised tensor shape is 3 2 4' assert_shape T_denoised; 3 2 4
Err_noisy=: albedo_tensor_tensor_fnorm (T_noisy - T_rank1)
Err_denoised=: albedo_tensor_tensor_fnorm (T_denoised - T_rank1)
'Denoised error is lower than noisy error' assert_true (Err_denoised < Err_noisy)

suite_end ''
