NB. ==============================================================================
NB. ALBEDO - Tests: detect/matched_subspace.ijs
NB. ==============================================================================

load 'detect/matched_subspace.ijs'

assert =: 3 : 0
  if. -. y do.
    echo 'ASSERTION FAILED'
    exit 1
  end.
  1
)

echo 'Testing Adaptive Matched Subspace Detector (AMSD)...'

NB. Bands: 8, Target subspace dim: 1, Background subspace dim: 2
b =: 8
pt =: 1
pb =: 2

NB. Target subspace (e.g. signature of specific mineral)
T =: (b, pt) $ 0.1 0.2 0.8 0.9 0.4 0.1 0.05 0.02

NB. Background clutter basis (2 eigenvectors of soil/vegetation)
B =: (b, pb) $ 0.5 0.1  0.6 0.2  0.7 0.3  0.6 0.4  0.5 0.5  0.4 0.6  0.3 0.7  0.2 0.8

NB. Generate 20 test pixels:
NB. Pixels 0-14: pure clutter B * alpha + small noise
NB. Pixels 15-19: target T * beta + clutter B * alpha + small noise
n =: 20
alpha =: (pb, n) $ 0.5
clutter =: |: B +/ .* alpha NB. (n x b)

NB. Add target to pixel 17
target_vec =: , T
X =: clutter
X =: ((17 { X) + 2.0 * target_vec) (< 17) } X

NB. Run AMSD detector
scores =: (T ; B) amsd_detector_albedo_amsd_ X
assert n = # scores

NB. Target pixel 17 must yield the maximum detection response
max_idx =: (i. >./) scores
assert 17 = max_idx

echo 'AMSD tests passed successfully.'
exit 0
