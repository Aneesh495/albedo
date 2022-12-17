NB. ==============================================================================
NB. ALBEDO - Tests: signal/spatial_filter.ijs
NB. ==============================================================================

load 'signal/spatial_filter.ijs'

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

echo 'Testing 2D spatial filtering and edge preservation...'

NB. Test dimensions of kernels
assert (3 3) -: $ albedo_spatial_SOBEL_H
assert (3 3) -: $ albedo_spatial_SOBEL_V
assert (3 3) -: $ albedo_spatial_LAPLACIAN_4

NB. Test Gaussian kernel generation (5x5, sigma=1.0)
g5 =: albedo_spatial_gaussian_kernel 2 1.0
assert (5 5) -: $ g5
assert 1.0 approx_eq (+/ , g5)

NB. Test Sobel vertical edge response
NB. Image with sharp vertical step in middle
img =: 8 8 $ 0.0
img =: 10.0 (<"1 (i. 8) ,"0/ 4 5 6 7) } img
'mag dir' =: albedo_spatial_sobel_gradient img
NB. Edge gradient should be strongly non-zero at column 3 and 4
assert *./ ((<"1 (i. 8) ,. 3) { mag) > 10.0

NB. Test bilateral filtering preserves sharp step while smoothing noise
step_img =: 10 10 $ 0.0
step_img =: 5.0 (<"1 (i. 10) ,"0/ 5 6 7 8 9) } step_img
smoothed =: (1 1.0 0.5) albedo_spatial_bilateral_filter2d step_img
NB. Values on either side should stay close to 0 and 5
assert *./ (((<"1 (i. 10) ,. 1) { smoothed) < 0.1)
assert *./ (((<"1 (i. 10) ,. 8) { smoothed) > 4.9)

echo 'Spatial filtering tests passed successfully.'
exit 0
