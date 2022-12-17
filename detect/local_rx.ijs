NB. ==============================================================================
NB. Albedo - Dual-Window Local RX Anomaly Detector
NB.
NB. Implements sliding dual-window local RX detection:
NB.   - Inner guard window (W_in x W_in): excludes target pixel and immediate
NB.     halo from contaminating the background model.
NB.   - Outer background window (W_out x W_out): gathers local clutter statistics
NB.     surrounding the guard region.
NB. Isolates subpixel targets in heterogeneous and non-stationary clutter.
NB. ==============================================================================

cocurrent 'base'

load 'core/linalg.ijs'
load 'core/stats.ijs'

NB. ------------------------------------------------------------------------------
NB. Local Dual-Window RX Detection
NB. ------------------------------------------------------------------------------

NB. Computes local RX anomaly score map for 3D cube (B x H x W)
NB. Dyadic: (w_in; w_out) detect_local_rx cube
detect_local_rx=: 4 : 0
  'w_in w_out'=. x
  cube=. y
  'B H W'=. $ cube
  assert. 1 = 2 | w_in [ 'Inner window size must be odd'
  assert. 1 = 2 | w_out [ 'Outer window size must be odd'
  assert. w_out > w_in [ 'Outer window must strictly exceed inner window'

  half_in=. <. w_in % 2
  half_out=. <. w_out % 2

  score_map=. (H, W) $ 0.0

  for_h. half_out + i. (H - 2 * half_out) do.
    for_w. half_out + i. (W - 2 * half_out) do.
      center_pixel=. (<a:; h; w) { cube

      NB. Gather outer window pixels excluding inner window
      bg_pixels=. 0 $ 0.0
      for_dh. (- half_out) + i. w_out do.
        for_dw. (- half_out) + i. w_out do.
          NB. Check if inside guard window
          in_guard=. ((| dh) <: half_in) *. ((| dw) <: half_in)
          if. 0 = in_guard do.
            sample=. (<a:; (h + dh); (w + dw)) { cube
            bg_pixels=. bg_pixels , sample
          end.
        end.
      end.

      N_bg=. (<. (# bg_pixels) % B)
      if. N_bg >: (B + 2) do.
        X_bg=. (B, N_bg) $ bg_pixels
        mu_local=. albedo_stats_mean_vector X_bg
        cov_local=. albedo_stats_covariance_matrix X_bg
        reg_cov=. cov_local + (1e_5 * albedo_linalg_eye B)
        inv_cov=. %. reg_cov

        score=. (mu_local; inv_cov) albedo_stats_mahalanobis_sq center_pixel
        score_map=. score (<h, w)} score_map
      end.
    end.
  end.

  score_map
)


albedo_local_rx_detect_local_rx =: detect_local_rx
