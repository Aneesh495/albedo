NB. ==============================================================================
NB. ALBEDO - Tests: subspace/ica.ijs
NB. ==============================================================================

load 'subspace/ica.ijs'

assert =: 3 : 0
  if. -. y do.
    echo 'ASSERTION FAILED'
    exit 1
  end.
  1
)

approx_eq =: 4 : 0
  tol =. 1e_4
  *./ , tol > | x - y
)

echo 'Testing FastICA Independent Component Analysis...'

NB. Create two independent 1D source signals:
NB. S1: sine wave, S2: sawtooth wave (length 100)
n =: 100
t =: (i. n) % 10.0
s1 =: 1 o. t
s2 =: 2.0 * ((t % 5.0) - <. t % 5.0) - 1.0
S =: s1 ,. s2 NB. (100 x 2)

NB. Mixing matrix A (2 x 2)
A =: 2 2 $ 0.8 0.6  0.4 0.9

NB. Mixed observations X = S * A^T
X =: S +/ .* |: A
assert (100 2) -: $ X

NB. Run FastICA for 2 components
'sources W W_white' =: 2 fastica_albedo_ica_ X

assert (100 2) -: $ sources
assert (2 2) -: $ W

NB. Verify reconstructed sources are mutually uncorrelated / orthogonal
c1 =: 0 {"1 sources
c2 =: 1 {"1 sources
cov12 =: (+/ (c1 - (+/ c1)%n) * (c2 - (+/ c2)%n)) % (<: n)
assert (| cov12) < 0.10

echo 'FastICA tests passed successfully.'
exit 0
