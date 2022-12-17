NB. ==============================================================================
NB. Test Suite: core/stats.ijs
NB. ==============================================================================

load 'tests/test_framework.ijs'
load 'core/stats.ijs'

cocurrent 'base'

suite_begin 'High-Dimensional Spectral Statistics'

NB. Test 1: Mean vector and Centering
X=: 2 4 $ 1 2 3 4 10 20 30 40
mu=: albedo_stats_mean_vector X
'Mean vector shape should be 2' assert_equal (# mu); 2
'Mean of row 0 should be 2.5' assert_near (0 { mu); 2.5; 1e_6
'Mean of row 1 should be 25.0' assert_near (1 { mu); 25.0; 1e_6

Xc=: albedo_stats_center_matrix X
'Centered matrix mean should be 0' assert_near (+/ , Xc); 0.0; 1e_6

NB. Test 2: Covariance Matrix
Cov=: albedo_stats_covariance_matrix X
'Covariance shape should be 2 2' assert_shape Cov; 2 2
'Covariance is symmetric' assert_near Cov; (albedo_linalg_tr Cov); 1e_6
'Variance of row 0 should be 5/3 ~ 1.66667' assert_near ((< 0 0) { Cov); (5.0 % 3.0); 1e_5

NB. Test 3: Correlation Matrix
Corr=: albedo_stats_correlation_matrix Cov
'Correlation diagonal should be 1.0' assert_near (albedo_linalg_diag_vec Corr); 1.0 1.0; 1e_6
'Correlation off-diagonals should be 1.0 (perfect correlation in synthetic X)' assert_near ((< 0 1) { Corr); 1.0; 1e_6

NB. Test 4: Shrinkage Covariance
ShrinkCov=: 0.2 albedo_stats_shrinkage_cov X
'Shrinkage cov shape should be 2 2' assert_shape ShrinkCov; 2 2
'Shrinkage trace matches original trace' assert_near (albedo_linalg_trace ShrinkCov); (albedo_linalg_trace Cov); 1e_6

NB. Test 5: Batch vs Single Mahalanobis Distance
NB. Generate random 3-band matrix with 10 samples
X_rand=: ? 3 10 $ 100
mu_r=: albedo_stats_mean_vector X_rand
cov_r=: albedo_stats_covariance_matrix X_rand
batch_d2=: (mu_r; cov_r) albedo_stats_batch_mahalanobis_sq X_rand
'Batch Mahalanobis count should be 10' assert_equal (# batch_d2); 10
NB. Verify sample 0 matches single vector calculation
reg_cov=: cov_r + (1e_7 * albedo_linalg_eye 3)
inv_cov=: %. reg_cov
single_d0=: (mu_r; inv_cov) albedo_stats_mahalanobis_sq (0 {"1 X_rand)
'Sample 0 Mahalanobis matches batch' assert_near (0 { batch_d2); single_d0; 1e_5

NB. Test 6: Median and MAD
Vec=: 1 2 3 4 5 6 7 8 9
'Median of 1..9 is 5' assert_equal (albedo_stats_median Vec); 5
'MAD of 1..9' assert_true ((albedo_stats_mad Vec) > 0.0)

suite_end ''
