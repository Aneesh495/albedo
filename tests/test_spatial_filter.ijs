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
assert (3 3) -: $ SOBEL_H_albedo_spatial_
assert (3 3) -: $ SOBEL_V_albedo_spatial_
assert (3 3) -: $ LAPLACIAN_4_albedo_spatial_

NB. Test Gaussian kernel generation (5x5, sigma=1.0)
g5 =: gaussian_kernel_albedo_spatial_ 2 1.0
assert (5 5) -: $ g5
assert 1.0 approx_eq (+/ , g5)

NB. Test Sobel vertical edge response
NB. Image with sharp vertical step in middle
img =: 8 8 $ 0.0
img =: 10.0 (<"1 (i. 8) ,"0/ 4 5 6 7) } img
'mag dir' =: sobel_gradient_albedo_spatial_ img
NB. Edge gradient should be strongly non-zero at column 3 and 4
assert *./ (mag {~ (i. 8) ,"0/ 3) > 10.0

NB. Test bilateral filtering preserves sharp step while smoothing noise
step_img =: 10 10 $ 0.0
step_img =: 5.0 (<"1 (i. 10) ,"0/ 5 6 7 8 9) } step_img
smoothed =: (1 1.0 0.5) bilateral_filter2d_albedo_spatial_ step_img
NB. Values on either side should stay close to 0 and 5
assert *./ ((smoothed {~ (i. 10) ,"0/ 1) < 0.1)
assert *./ ((smoothed {~ (i. 10) ,"0/ 8) > 4.9)

echo 'Spatial filtering tests passed successfully.'
exit 0
