NB. ==============================================================================
NB. Test Suite: tensor/robust_tensor_pca.ijs
NB. ==============================================================================

load 'tests/test_framework.ijs'
load 'tensor/robust_tensor_pca.ijs'
load 'core/tensor.ijs'

cocurrent 'base'

suite_begin 'Robust Tensor Principal Component Analysis (RTPCA)'

NB. Create low-rank 3D tensor L (rank-1 outer product)
a=: 1.0 2.0 3.0 4.0
b=: 1.0 1.0 1.0
c=: 2.0 2.0 2.0
L_true=: (a; b; c) albedo_tensor_vec_outer3 ''
'L_true shape is 4 3 3' assert_shape L_true; 4 3 3

NB. Add sparse spike anomaly at (1, 1, 1)
S_true=: 4 3 3 $ 0.0
S_true=: 50.0 (<0, 1, 1)} S_true
S_true=: 50.0 (<1, 1, 1)} S_true

X_total=: L_true + S_true

NB. Run RTPCA decomposition
'L_est S_est'=: (0.1; 20; 1e_4) albedo_rtpca_fit_rtpca X_total
'L_est shape is 4 3 3' assert_shape L_est; 4 3 3
'S_est shape is 4 3 3' assert_shape S_est; 4 3 3

NB. Verify sparse tensor captures anomaly spike at (0, 1, 1) and (1, 1, 1)
'Sparse tensor isolates anomaly spike' assert_true ((S_est < 0 1 1) > 20.0)
'Sparse tensor is zero away from anomaly' assert_near (S_est < 0 0 0); 0.0; 1.0

NB. Verify sum L + S reconstructs X_total with low error
RecErr=: albedo_tensor_tensor_fnorm (X_total - (L_est + S_est))
'RTPCA reconstruction error is low' assert_true (RecErr < 10.0)

suite_end ''
