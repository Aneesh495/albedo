NB. ==============================================================================
NB. ALBEDO - Tests: detect/kernel_rx.ijs
NB. ==============================================================================

load 'detect/kernel_rx.ijs'

assert =: 3 : 0
  if. -. y do.
    echo 'ASSERTION FAILED'
    exit 1
  end.
  1
)

echo 'Testing Kernel Reed-Xiaoli nonlinear anomaly detector...'

NB. Synthesize 20 background pixels (5 bands) clustered around mean
n =: 20
b =: 5
bg =: (n, b) $ (0.1 * i. n * b) % (n * b)

NB. Inject a distinct nonlinear anomaly at index 7
anomaly =: 2.5 3.0 _1.5 4.0 2.0
X =: anomaly (< 7) } bg
assert (n, b) -: $ X

NB. Test Gram matrix calculation
K =: 0.1 rbf_gram_matrix_albedo_krx_ X
assert (n, n) -: $ K
NB. Diagonal elements must be exp(0) = 1.0
assert *./ 1e_6 > | 1.0 - (<"1 (i. n) ,. (i. n)) { K

NB. Test Gram matrix centering
K_tilde =: center_gram_matrix_albedo_krx_ K
assert 1e_6 > | +/ , K_tilde

NB. Run Kernel RX detection (gamma=0.5, lambda=1e-3)
scores =: (0.5 1e_3) kernel_rx_albedo_krx_ X
assert n = # scores

NB. The injected anomaly at index 7 must produce the highest anomaly score
max_idx =: (i. >./) scores
assert 7 = max_idx

echo 'Kernel RX tests passed successfully.'
exit 0
