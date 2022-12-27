NB. ==============================================================================
NB. ALBEDO - High-Dimensional Hyperspectral Imaging & Tensor Decomposition Engine
NB. albedo.ijs : Master Entrypoint and Interactive Environment
NB.
NB. Loads all core subsystems:
NB.   - core/     : Dense linear algebra, tensors, FFT, wavelets, splines
NB.   - io/       : ENVI format parser/serializer, spectral library, 3D cube container
NB.   - atmos/    : Radiative transfer, TOA reflectance, DOS, QUAC, continuum removal
NB.   - subspace/ : PCA, MNF, HySime, FastICA, segmented PCA
NB.   - unmix/    : VCA, N-FINDR, PPI, NNLS, FCLS, SUNSAL, NMF, Bundles (ELMM)
NB.   - detect/   : Reed-Xiaoli (RX), Dual-Window Local RX, Kernel RX, CEM, MF, ACE, AMSD
NB.   - tensor/   : Tucker3 HOSVD, CP-ALS, RTPCA (ADMM), low-rank denoising
NB.   - indices/  : Narrowband vegetation, water, and mineral diagnostic indices
NB.   - classify/ : SAM classifier, K-Means++, spatial-spectral voting, graph Laplacians
NB.   - signal/   : 1D/2D DWT, Mallat pyramid, spatial filters, Grunwald-Letnikov
NB.   - data/     : Calibrated USGS minerals, vegetation, urban, ECOSTRESS, HITRAN gases
NB.   - sim/      : Radiative transfer synthetic scene generator with trace gas plumes
NB.   - pipeline/ : Automated unmixing and anomaly detection pipelines
NB. ==============================================================================

cocurrent 'base'

albedo_version =: '1.0.0'
albedo_release_date =: '2026-11-22'

NB. Load Core Math Primitives
load 'core/linalg.ijs'
load 'core/linalg_ext.ijs'
load 'core/stats.ijs'
load 'core/tensor.ijs'
load 'core/fft.ijs'
load 'core/morphology.ijs'
load 'core/spline.ijs'
load 'core/quadrature.ijs'

NB. Load I/O & Cube Container
load 'io/envi.ijs'
load 'io/spectral_lib.ijs'
load 'io/cube.ijs'

NB. Load Atmospheric Calibration & Corrections
load 'atmos/solar.ijs'
load 'atmos/radiance.ijs'
load 'atmos/dos.ijs'
load 'atmos/quac.ijs'
load 'atmos/continuum.ijs'

NB. Load Subspace Dimensionality Reduction
load 'subspace/pca.ijs'
load 'subspace/mnf.ijs'
load 'subspace/hysime.ijs'
load 'subspace/segmented_pca.ijs'
load 'subspace/ica.ijs'

NB. Load Spectral Unmixing Algorithms
load 'unmix/vca.ijs'
load 'unmix/nfindr.ijs'
load 'unmix/ppi.ijs'
load 'unmix/nnls.ijs'
load 'unmix/fcls.ijs'
load 'unmix/sunsal.ijs'
load 'unmix/endmember_bundles.ijs'
load 'unmix/nmf.ijs'

NB. Load Target & Anomaly Detection Kernels
load 'detect/rx.ijs'
load 'detect/local_rx.ijs'
load 'detect/matched_filter.ijs'
load 'detect/sam.ijs'
load 'detect/kernel_rx.ijs'
load 'detect/matched_subspace.ijs'

NB. Load Multi-way Tensor Decompositions
load 'tensor/cp_als.ijs'
load 'tensor/tucker.ijs'
load 'tensor/denoise.ijs'
load 'tensor/robust_tensor_pca.ijs'

NB. Load Indices & Classification
load 'indices/vegetation.ijs'
load 'indices/minerals_water.ijs'
load 'classify/sam_classifier.ijs'
load 'classify/kmeans.ijs'
load 'classify/spatial_spectral.ijs'
load 'classify/spectral_graph.ijs'

NB. Load Signal & Spatial Transforms
load 'signal/wavelet.ijs'
load 'signal/spatial_filter.ijs'
load 'signal/fractional_derivative.ijs'

NB. Load Reference Spectral Endmember Libraries
load 'data/usgs_minerals.ijs'
load 'data/usgs_vegetation.ijs'
load 'data/usgs_urban.ijs'
load 'data/atmospheric_lut.ijs'
load 'data/ecostress_rocks.ijs'
load 'data/ecostress_soils.ijs'
load 'data/polymers_plastics.ijs'
load 'data/trace_gases.ijs'
load 'data/aquatic_optics.ijs'

NB. Load Simulation & Pipelines
load 'sim/synthetic_scene.ijs'
load 'pipeline/pipeline_unmix.ijs'
load 'pipeline/pipeline_detect.ijs'

NB. Banner
echo '========================================================================'
echo '        ALBEDO: Hyperspectral Imaging & Tensor Decomposition Engine    '
echo '        Version: ', albedo_version, ' | Pure J Vectorized Numerics      '
echo '========================================================================'
echo 'Modules initialized: core, io, atmos, subspace, unmix, detect, tensor,'
echo '                     indices, classify, signal, data, sim, pipeline.'
echo 'Type: albedo_help '''' for interactive reference.'

albedo_help =: 3 : 0
  echo '------------------------------------------------------------------------'
  echo 'ALBEDO QUICK REFERENCE'
  echo '------------------------------------------------------------------------'
  echo 'Simulate Scene:    scene =: albedo_sim_generate_synthetic_scene 32 32'
  echo 'Extract Endmembers: E =: 4 albedo_vca_extract_vca (albedo_cube_cube_data scene)'
  echo 'Estimate Abundances: A =: (E; 100) albedo_fcls_fcls_cube (albedo_cube_cube_data scene)'
  echo 'Global RX Anomaly: scores =: albedo_rx_detect_global_rx (albedo_cube_cube_data scene)'
  echo 'Run Test Suite:    (load ''scripts/run_tests.py'' via shell)'
  echo 'Run Benchmarks:    (load ''scripts/benchmarks.ijs'')'
  echo '------------------------------------------------------------------------'
)
