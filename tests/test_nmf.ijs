NB. ==============================================================================
NB. ALBEDO - Tests: unmix/nmf.ijs
NB. ==============================================================================

load 'unmix/nmf.ijs'

assert =: 3 : 0
  if. -. y do.
    echo 'ASSERTION FAILED'
    exit 1
  end.
  1
)

echo 'Testing Non-negative Matrix Factorization spectral unmixing...'

NB. Ground truth: 2 endmembers, 10 bands, 20 pixels
b =: 10
p =: 2
n =: 20

M_true =: 10 2 $ 0.8 0.1  0.7 0.2  0.6 0.3  0.5 0.4  0.4 0.5  0.3 0.6  0.2 0.7  0.1 0.8  0.2 0.7  0.3 0.6
NB. True abundances satisfying ANC and ASC
a1 =: (i. n) % (<: n)
a2 =: 1.0 - a1
S_true =: a1 ,: a2

NB. Linear mixture Y = M * S
Y =: M_true +/ .* S_true
assert (b, n) -: $ Y

NB. Run NMF for 2 endmembers, 60 iterations
'M_est S_est loss' =: (p ; 60) albedo_nmf_nmf_unmix Y

assert (b, p) -: $ M_est
assert (p, n) -: $ S_est

NB. Non-negativity constraints strictly satisfied
assert *./ , M_est >: 0.0
assert *./ , S_est >: 0.0

NB. Low reconstruction loss
assert loss < 0.5

echo 'NMF spectral unmixing tests passed successfully.'
exit 0
