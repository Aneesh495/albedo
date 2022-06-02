NB. ==============================================================================
NB. Test Suite: core/fft.ijs
NB. ==============================================================================

load 'tests/test_framework.ijs'
load 'core/fft.ijs'

cocurrent 'base'

suite_begin 'Spectral Signal Processing & FFT'

NB. Test 1: 1D FFT and IFFT roundtrip
Signal1D=: 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0
Spec1D=: albedo_fft_fft1d Signal1D
'Spec1D length should be 8' assert_equal (# Spec1D); 8
Rec1D=: albedo_fft_ifft1d Spec1D
'IFFT(FFT(x)) matches x' assert_near Rec1D; Signal1D; 1e_6

NB. Test 2: Parseval's Theorem
NB. Sum of squared magnitudes in time domain equals (1/N) * sum of squared magnitudes in freq domain
EnergyTime=: +/ *: Signal1D
EnergyFreq=: (+/ *: | Spec1D) % 8
'Parsevals theorem holds' assert_near EnergyTime; EnergyFreq; 1e_5

NB. Test 3: 2D Spatial FFT and IFFT roundtrip
Image2D=: 4 4 $ i. 16
Spec2D=: albedo_fft_fft2d Image2D
'Spec2D shape should be 4 4' assert_shape Spec2D; 4 4
Rec2D=: albedo_fft_ifft2d Spec2D
'2D IFFT roundtrip matches Image2D' assert_near Rec2D; Image2D; 1e_6

NB. Test 4: Hann and Hamming windows
Hann8=: albedo_fft_hann_window 8
'Hann window length 8' assert_equal (# Hann8); 8
'Hann window symmetric' assert_near Hann8; (|. Hann8); 1e_6
'Hann endpoints are zero' assert_near (0 { Hann8); 0.0; 1e_6

Ham8=: albedo_fft_hamming_window 8
'Hamming window length 8' assert_equal (# Ham8); 8
'Hamming window symmetric' assert_near Ham8; (|. Ham8); 1e_6

NB. Test 5: Savitzky-Golay Smoothing preserves quadratic polynomial exactly
NB. Quadratic signal y = x^2 + 2x + 1
X_vals=: i. 15
QuadSignal=: 1.0 + (2.0 * X_vals) + (*: X_vals)
Smoothed=: (5; 2; 0) albedo_fft_sg_filter QuadSignal
'SG smoothing preserves exact quadratic polynomial' assert_near Smoothed; QuadSignal; 1e_5

NB. Test 6: Savitzky-Golay 1st derivative of linear ramp y = 3x is constant 3
RampSignal=: 3.0 * X_vals
Deriv1=: (5; 2; 1) albedo_fft_sg_filter RampSignal
'Interior 1st derivative of 3x is 3.0' assert_near (2 }. _2 }. Deriv1); (11 $ 3.0); 1e_5

suite_end ''
