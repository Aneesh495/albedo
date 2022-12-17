NB. ==============================================================================
NB. ALBEDO - High-Dimensional Hyperspectral Imaging Engine
NB. subspace/ica.ijs : FastICA Independent Component Analysis
NB. Implements Hyvärinen fixed-point FastICA with negentropy contrast functions (tanh),
NB. Gram-Schmidt deflationary orthogonalization, and spectral source separation.
NB. ==============================================================================

cocurrent 'base'

load 'core/linalg.ijs'
load 'core/stats.ijs'

NB. Center and whiten data matrix X (samples x bands)
NB. Returns boxed: (X_whitened ; mean_vec ; whitening_mat)
whiten =: 3 : 0
  x =. y
  'n p' =. $ x
  mu =. (mean_vector_albedo_stats_ x)
  NB. Center
  xc =. x -"1 mu
  NB. Covariance
  cov =. ((|: xc) +/ .* xc) % (<: n)
  NB. Eigen-decomposition via SVD of xc: xc = U S V^T
  NB. cov = V (S^2 / (n-1)) V^T
  'u s vt' =. thin_svd_albedo_linalg_ xc
  NB. Singular values of covariance: lambda = (s^2) / (n-1)
  lambdas =. (*: s) % (<: n)
  NB. Keep non-zero components
  inv_sqrt_d =. % %: lambdas + 1e_8
  NB. Whitening matrix: W_white = V * diag(inv_sqrt_d)
  w_white =. (|: vt) * inv_sqrt_d
  x_white =. xc +/ .* w_white
  x_white ; mu ; w_white
)

NB. Hyperbolic tangent contrast function g(u) = tanh(u)
NB. and derivative g'(u) = 1 - tanh(u)^2
NB. J hyperbolic tangent: 7 o. u
contrast_tanh =: 3 : 0
  t =. 7 o. y
  dg =. 1.0 - *: t
  t ; dg
)

NB. FastICA fixed-point iteration for a single unit vector w
NB. x_white: whitened matrix (n x p), w_init: initial unit vector (p)
NB. w_prev: matrix of previously converged components (k x p)
fastica_one_component =: 4 : 0
  'x_white w_prev' =. x
  w =. y
  'n p' =. $ x_white
  max_iter =. 100
  tol =. 1e_5

  for_iter. i. max_iter do.
    NB. Project: u = x_white +/ .* w
    u =. x_white +/ .* w
    'g_u dg_u' =. contrast_tanh u

    NB. Term 1: E[x * g(w^T x)] = (x_white^T g_u) / n
    term1 =. ((|: x_white) +/ .* g_u) % n

    NB. Term 2: E[g'(w^T x)] * w = (sum(dg_u) / n) * w
    term2 =. ((+/ dg_u) % n) * w

    w_new =. term1 - term2

    NB. Deflationary Gram-Schmidt orthogonalization against previous vectors
    if. 0 < # w_prev do.
      for_k. i. # w_prev do.
        wk =. k { w_prev
        proj =. (+/ w_new * wk) * wk
        w_new =. w_new - proj
      end.
    end.

    NB. Normalize
    norm =. %: +/ *: w_new
    if. norm < 1e_10 do.
      w_new =. (? p # 0) - 0.5
      norm =. %: +/ *: w_new
    end.
    w_new =. w_new % norm

    NB. Check convergence: |w_new . w| approx 1
    cos_sim =. | +/ w_new * w
    w =. w_new
    if. (1.0 - cos_sim) < tol do.
      break.
    end.
  end.
  w
)

NB. Run FastICA to extract k independent components
NB. Input: x (n samples x p bands), k: number of components
NB. Returns boxed: (sources_S ; unmixing_W ; whitening_mat)
fastica =: 4 : 0
  k =. x
  mat =. y
  'xw mu w_white' =. whiten mat
  'n p' =. $ xw
  NB. Limit k to min(k, p)
  k =. k <. p

  NB. Demixing matrix W_demix (k x p)
  w_accum =. 0 p $ 0.0

  for_comp. i. k do.
    NB. Deterministic seed or basis vector for initialization
    w_init =. p # 0.0
    w_init =. 1.0 comp} w_init
    w_res =. (xw ; w_accum) fastica_one_component w_init
    w_accum =. w_accum , w_res
  end.

  NB. Source signals: S = xw * W_demix^T (n x k)
  sources =. xw +/ .* |: w_accum
  sources ; w_accum ; w_white
)


NB. Exported Aliases
albedo_ica_fastica =: fastica
