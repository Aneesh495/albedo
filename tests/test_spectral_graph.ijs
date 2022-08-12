NB. ==============================================================================
NB. Test Suite: classify/spectral_graph.ijs
NB. ==============================================================================

load 'tests/test_framework.ijs'
load 'classify/spectral_graph.ijs'
load 'core/linalg.ijs'

cocurrent 'base'

suite_begin 'Spectral Graph Laplacian & Manifold Clustering'

NB. Create two well-separated spectral clusters (4 bands, 20 samples total)
Cluster1=: 4 10 $ ? 40 $ 10
Cluster2=: 4 10 $ 100 + ? 40 $ 10
X_two_clusters=: Cluster1 ,. Cluster2

NB. Test 1: Affinity Matrix Construction
W=: 10.0 albedo_spectral_graph_compute_rbf_affinity X_two_clusters
'Affinity matrix shape is 20 20' assert_shape W; 20 20
'Affinity matrix is symmetric' assert_near W; (albedo_linalg_tr W); 1e_6
'Affinity diagonal is zero' assert_near (albedo_linalg_diag_vec W); (20 $ 0.0); 1e_6

NB. Test 2: Normalized Laplacian Properties
'L_sym inv_sqrt_d'=: albedo_spectral_graph_normalized_laplacian W
'Laplacian shape is 20 20' assert_shape L_sym; 20 20
'Laplacian is symmetric' assert_near L_sym; (albedo_linalg_tr L_sym); 1e_6

NB. Test 3: Spectral Embedding
Y_embed=: 2 albedo_spectral_graph_spectral_embedding X_two_clusters
'Embedding shape is 2 20' assert_shape Y_embed; 2 20

NB. Test 4: Spectral Graph Clustering
'Labels Y_norm'=: (2; 10.0) albedo_spectral_graph_cluster_spectral X_two_clusters
'Labels count is 20' assert_equal (# Labels); 20
NB. Verify cluster 1 points are assigned same label, cluster 2 points assigned opposite label
'Cluster 1 points have uniform label' assert_equal (# ~. 10 {. Labels); 1
'Cluster 2 points have uniform label' assert_equal (# ~. 10 }. Labels); 1
'Cluster 1 and Cluster 2 labels differ' assert_true ((0 { Labels) ~: (10 { Labels))

suite_end ''
