NB. ==============================================================================
NB. ALBEDO - High-Dimensional Hyperspectral Imaging Engine
NB. signal/spatial_filter.ijs : 2D Spatial Filtering & Edge-Preserving Smoothing
NB. Implements 2D spatial convolution, Sobel/Prewitt/Laplacian edge detectors,
NB. Gaussian blurring, and bilateral filtering for spatial-spectral preprocessing.
NB. ==============================================================================

cocurrent 'base'

NB. Standard 3x3 Differential Kernels
SOBEL_H =: 3 3 $ _1 0 1  _2 0 2  _1 0 1
SOBEL_V =: 3 3 $ _1 _2 _1  0 0 0  1 2 1

PREWITT_H =: 3 3 $ _1 0 1  _1 0 1  _1 0 1
PREWITT_V =: 3 3 $ _1 _1 _1  0 0 0  1 1 1

LAPLACIAN_4 =: 3 3 $ 0 1 0  1 _4 1  0 1 0
LAPLACIAN_8 =: 3 3 $ 1 1 1  1 _8 1  1 1 1

NB. 2D Convolution with periodic / zero-padded boundaries
NB. x: 2D kernel (kh x kw), y: 2D image (h x w)
convolve2d =: 4 : 0
  kernel =. x
  img =. y
  'kh kw' =. $ kernel
  'h w' =. $ img
  pad_h =. <. -: kh
  pad_w =. <. -: kw
  out =. (h, w) $ 0.0

  for_r. i. h do.
    for_c. i. w do.
      acc =. 0.0
      for_kr. i. kh do.
        rr =. r + kr - pad_h
        if. (rr >: 0) *. (rr < h) do.
          for_kc. i. kw do.
            cc =. c + kc - pad_w
            if. (cc >: 0) *. (cc < w) do.
              acc =. acc + ((< rr, cc) { img) * ((< kr, kc) { kernel)
            end.
          end.
        end.
      end.
      out =. acc (< r , c) } out
    end.
  end.
  out
)

NB. Generate 2D Gaussian Kernel (size 2k+1, standard deviation sigma)
NB. Input: (half_width ; sigma)
gaussian_kernel =: 3 : 0
  'hw sigma' =. y
  dim =. 1 + 2 * hw
  coords =. (i. dim) - hw
  NB. Outer sum of squares
  d2 =. (coords */ coords) + (coords */ coords) NB. will create dim x dim
  NB. More direct:
  r_mat =. (dim, dim) $ coords
  c_mat =. |: r_mat
  r2 =. (*: r_mat) + (*: c_mat)
  k =. ^ - r2 % 2.0 * *: sigma
  k % +/ , k
)

NB. Spatial gradient magnitude and direction (Sobel)
sobel_gradient =: 3 : 0
  gx =. SOBEL_H convolve2d y
  gy =. SOBEL_V convolve2d y
  mag =. %: (*: gx) + (*: gy)
  dir =. _3 o. gy % gx + (gx = 0.0) NB. arctan2(gy, gx) in J
  mag ; dir
)

NB. 2D Bilateral Filter (Edge-Preserving Spatial-Spectral Smoothing)
NB. Preserves sharp geological / anthropogenic boundaries while smoothing sensor noise
NB. Input: img (2D), x: (half_window ; sigma_spatial ; sigma_range)
bilateral_filter2d =: 4 : 0
  'hw sig_s sig_r' =. x
  img =. y
  'h w' =. $ img
  dim =. 1 + 2 * hw
  out =. (h, w) $ 0.0

  two_sig_s2 =. 2.0 * *: sig_s
  two_sig_r2 =. 2.0 * *: sig_r

  for_r. i. h do.
    for_c. i. w do.
      center_val =. (< r, c) { img
      w_sum =. 0.0
      val_sum =. 0.0
      for_dr. (i. dim) - hw do.
        rr =. r + dr
        if. (rr >: 0) *. (rr < h) do.
          for_dc. (i. dim) - hw do.
            cc =. c + dc
            if. (cc >: 0) *. (cc < w) do.
              neighbor_val =. (< rr, cc) { img
              spatial_dist2 =. (*: dr) + (*: dc)
              range_dist2 =. *: neighbor_val - center_val
              spatial_weight =. ^ - spatial_dist2 % two_sig_s2
              range_weight   =. ^ - range_dist2 % two_sig_r2
              weight =. spatial_weight * range_weight
              w_sum =. w_sum + weight
              val_sum =. val_sum + weight * neighbor_val
            end.
          end.
        end.
      end.
      out =. (val_sum % (w_sum + 1e_8)) (< r , c) } out
    end.
  end.
  out
)


NB. Exported Aliases
albedo_spatial_SOBEL_H =: SOBEL_H
albedo_spatial_SOBEL_V =: SOBEL_V
albedo_spatial_LAPLACIAN_4 =: LAPLACIAN_4
albedo_spatial_gaussian_kernel =: gaussian_kernel
albedo_spatial_sobel_gradient =: sobel_gradient
albedo_spatial_bilateral_filter2d =: bilateral_filter2d
