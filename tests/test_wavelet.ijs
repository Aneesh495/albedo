NB. ==============================================================================
NB. ALBEDO - Tests: signal/wavelet.ijs
NB. ==============================================================================

load 'signal/wavelet.ijs'

assert =: 3 : 0
  if. -. y do.
    echo 'ASSERTION FAILED'
    exit 1
  end.
  1
)

approx_eq =: 4 : 0
  tol =. 1e_5
  *./ , tol > | x - y
)

echo 'Testing 1D and 2D Discrete Wavelet Transforms...'

NB. Test 1D Haar forward and inverse perfect reconstruction
sig1 =: 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0
'a1 d1' =: albedo_wavelet_dwt1_haar_step sig1
recon1 =: albedo_wavelet_idwt1_haar_step (a1 ; d1)
assert sig1 approx_eq recon1

NB. Test multi-level 1D Haar decomposition (3 levels on length 16)
sig16 =: 1.0 + i. 16
decomp3 =: 3 albedo_wavelet_dwt1_multilevel sig16
assert 4 = # decomp3
recon16 =: albedo_wavelet_idwt1_multilevel decomp3
assert sig16 approx_eq recon16

NB. Test soft thresholding
vals =: _3.0 _1.5 0.0 1.5 3.0
th =: 2.0 albedo_wavelet_soft_threshold vals
assert (_1.0 0.0 0.0 0.0 1.0) approx_eq th

NB. Test hard thresholding
h_th =: 2.0 albedo_wavelet_hard_threshold vals
assert (_3.0 0.0 0.0 0.0 3.0) approx_eq h_th

NB. Test 1D Wavelet denoising on noisy step function
clean =: (16 # 2.0) , (16 # 8.0)
noise =: 0.1 * _1.0 1.0 _1.0 1.0 _1.0 1.0 _1.0 1.0 _1.0 1.0 _1.0 1.0 _1.0 1.0 _1.0 1.0 _1.0 1.0 _1.0 1.0 _1.0 1.0 _1.0 1.0 _1.0 1.0 _1.0 1.0 _1.0 1.0 _1.0 1.0
noisy =: clean + noise
denoised =: albedo_wavelet_wavelet_denoise_1d noisy
err_noisy =: +/ *: clean - noisy
err_clean =: +/ *: clean - denoised
assert err_clean < err_noisy

echo 'Wavelet transform tests passed successfully.'
exit 0
