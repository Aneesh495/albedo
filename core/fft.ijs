NB. ==============================================================================
NB. Albedo - Spectral Signal Processing & Fourier Analysis
NB.
NB. Implements Cooley-Tukey Radix-2 FFT, 2D Spatial FFT, power spectral density,
NB. spectral windowing (Hann, Hamming, Blackman), and Savitzky-Golay polynomial
NB. smoothing and derivative filters for reflectance spectroscopy.
NB. ==============================================================================

cocurrent 'base'

NB. Load linear algebra dependencies
load 'core/linalg.ijs'

NB. Complex imaginary unit constant
PI=: 3.14159265358979323846
TWO_PI=: 2.0 * PI

NB. ------------------------------------------------------------------------------
NB. 1D Radix-2 Cooley-Tukey Fast Fourier Transform
NB. ------------------------------------------------------------------------------

NB. Forward Discrete Fourier Transform (DFT) matrix or recursive radix-2
NB. Computes FFT of complex vector y of length N (where N is power of 2)
fft1d=: 3 : 0
  x=. y
  N=. # x
  if. N <: 1 do. x return. end.

  assert. 0 = 2 | N [ 'FFT length must be even (power of 2)'

  NB. Even and odd elements
  even=. fft1d (0 = 2 | i. N) # x
  odd=.  fft1d (1 = 2 | i. N) # x

  NB. Twiddle factors W_N^k = exp(-2 pi i k / N)
  k=. i. -: N
  angles=. (_1.0 * TWO_PI * k) % N
  NB. In J, complex number from polar (magnitude, angle) is r * (cos theta + i sin theta)
  NB. or using j. primitive
  twiddles=. (2&o. angles) + (0j1 * 1&o. angles)

  t=. twiddles * odd
  (even + t) , (even - t)
)

NB. Inverse 1D FFT: IFFT(X) = (1/N) * conj(FFT(conj(X)))
ifft1d=: 3 : 0
  X=. y
  N=. # X
  conj_X=. + X
  transformed=. fft1d conj_X
  NB. Complex conjugate and scale by 1/N
  (+ transformed) % N
)

NB. Power Spectral Density (PSD): |FFT(x)|^2 / N
psd1d=: 3 : 0
  x=. y
  N=. # x
  fx=. fft1d x
  (*: | fx) % N
)

NB. ------------------------------------------------------------------------------
NB. 2D Spatial FFT & IFFT
NB. ------------------------------------------------------------------------------

NB. 2D FFT on 2D matrix (separable row-wise then column-wise)
fft2d=: 3 : 0
  M=. y
  'H W'=. $ M
  NB. Row-wise 1D FFT
  M_rows=. fft1d "1 M
  NB. Column-wise 1D FFT (transpose -> row FFT -> transpose)
  |: fft1d "1 |: M_rows
)

NB. 2D Inverse FFT
ifft2d=: 3 : 0
  M=. y
  'H W'=. $ M
  M_rows=. ifft1d "1 M
  |: ifft1d "1 |: M_rows
)

NB. ------------------------------------------------------------------------------
NB. Spectral Windows
NB. ------------------------------------------------------------------------------

NB. Hann window of length N: 0.5 * (1 - cos(2 pi n / (N - 1)))
hann_window=: 3 : 0
  N=. y
  if. N <: 1 do. 1 $ 1.0 return. end.
  n=. i. N
  0.5 * 1.0 - 2&o. (TWO_PI * n) % (N - 1)
)

NB. Hamming window of length N: 0.54 - 0.46 * cos(2 pi n / (N - 1))
hamming_window=: 3 : 0
  N=. y
  if. N <: 1 do. 1 $ 1.0 return. end.
  n=. i. N
  0.54 - 0.46 * 2&o. (TWO_PI * n) % (N - 1)
)

NB. Blackman window of length N
blackman_window=: 3 : 0
  N=. y
  if. N <: 1 do. 1 $ 1.0 return. end.
  n=. i. N
  a0=. 0.42
  a1=. 0.5
  a2=. 0.08
  (a0 - a1 * 2&o. (TWO_PI * n) % (N - 1)) + (a2 * 2&o. (4.0 * PI * n) % (N - 1))
)

NB. ------------------------------------------------------------------------------
NB. Savitzky-Golay Polynomial Smoothing & Derivative Filter
NB. ------------------------------------------------------------------------------

NB. Computes Savitzky-Golay convolution coefficients
NB. window_size must be odd. poly_order is degree of polynomial (typically 2, 3, or 4).
NB. deriv_order is 0 for smoothing, 1 for first derivative, 2 for second derivative.
NB. Dyadic: (window_size; poly_order; deriv_order) savitzky_golay_coeffs ''
savitzky_golay_coeffs=: 3 : 0
  'window_size poly_order deriv_order'=. y
  assert. 1 = 2 | window_size [ 'Savitzky-Golay window size must be odd'
  assert. window_size > poly_order [ 'Window size must exceed polynomial order'

  half_w=. <. window_size % 2
  z=. (- half_w) + i. window_size

  NB. Vandermonde matrix J of shape (window_size, poly_order + 1): J_ij = z_i^j
  J=. z */ ^ (i. poly_order + 1)

  NB. Gram matrix J^T J
  JtJ=. (albedo_linalg_tr J) albedo_linalg_matmul J
  inv_JtJ=. %. JtJ

  NB. Filter matrix C = (J^T J)^(-1) J^T
  C=. inv_JtJ albedo_linalg_matmul (albedo_linalg_tr J)

  NB. Select row corresponding to deriv_order, scaled by factorial(deriv_order)
  fact=. ! deriv_order
  fact * (deriv_order { C)
)

NB. Apply Savitzky-Golay filter to a 1D spectral curve
NB. Reflects boundaries to prevent edge distortion.
NB. Dyadic: (window_size; poly_order; deriv_order) sg_filter spectral_curve
sg_filter=: 4 : 0
  'window_size poly_order deriv_order'=. x
  curve=. y
  n=. # curve
  half_w=. <. window_size % 2
  coeffs=. savitzky_golay_coeffs (window_size; poly_order; deriv_order)

  NB. Pad curve with symmetric boundary reflection
  left_pad=. |. (1 + i. half_w) { curve
  right_pad=. |. (_2 - i. half_w) { curve
  padded=. left_pad , curve , right_pad

  NB. Convolution: compute dot product of window with coefficients
  res=. n $ 0.0
  for_i. i. n do.
    sub=. (i + i. window_size) { padded
    val=. +/ sub * coeffs
    res=. val (<i)} res
  end.
  res
)


albedo_fft_fft1d =: fft1d
albedo_fft_ifft1d =: ifft1d
albedo_fft_fft2d =: fft2d
albedo_fft_ifft2d =: ifft2d
albedo_fft_hann_window =: hann_window
albedo_fft_hamming_window =: hamming_window
albedo_fft_sg_filter =: sg_filter
