###########################################################################
#                           MONTE CARLO STUDY
#   Estimating the mean function in sparse longitudinal functional data
###########################################################################

# =========================================================================
#  CONFIGURATION
# =========================================================================

## ---- Mean function  mu(u,t) --------------------------------------------
P   <- 8L                       # number of Fourier terms
a   <- 1;  b <- 1;  cc <- pi    # cc avoids clash with base::c()
## The default uses the parameters P, a, b, cc declared just above.
MU_FUN <- function(u, t) {
  nu <- length(u)
  S  <- matrix(0, nu, length(t))                 # start from the zero matrix
  for (p in 1:P) {
    # outer(2*pi*p*u, cc*t, "+")[i,j] = 2*pi*p*u[i] + cc*t[j]   (an nu x nt grid)
    S <- S + ((-1)^p / p^2) * sin(outer(2 * pi * p * u, cc * t, "+"))
  }
  # scale column j by (a*t[j] + b*t[j]^2); 'byrow' repeats that vector down rows
  S * matrix(a * t + b * t^2, nu, length(t), byrow = TRUE)
}

## ---- Second mean of higher boundary smoothness
WIN_B    <- function(t) 16 * t^2 * (1 - t)^2
MU_FUN_B <- function(u, t) MU_FUN(u, t) * matrix(WIN_B(t), length(u), length(t), byrow = TRUE)

## ---- Observation design: density g(t) of the visit times ---------------
## Visit times are i.i.d. draws with density g on [0,1]. A design is a list
## with four entries:
##   g(t)  : density (must be > 0 on [0,1]);
##   G(t)  : its integral G(t) = \int_0^t g  (the c.d.f.);
##   gmax  : any number >= max_t g(t)  (used by the rejection sampler);
##   tex   : how the density is written inside LaTeX captions.
## Default: a non-uniform cosine density.
G_DENSITY <- function(t) 1 + 0.6 * cos(2 * pi * t)
G_CDF     <- function(t) t + (0.6 / (2 * pi)) * sin(2 * pi * t)
G_MAX     <- 1.6
G_TEX     <- "g(t)=1 + 3 \\cos(2\\pi t) / 5"   # how g is written in captions
DES       <- list(g = G_DENSITY, G = G_CDF, gmax = G_MAX, tex = G_TEX)
## Skewed alternative used by the robustness study (same g_min = 2/5):
DES_SKEW  <- list(g    = function(t) 0.4 + 1.2 * t,
                  G    = function(t) 0.4 * t + 0.6 * t^2,
                  gmax = 1.6,
                  tex  = "g_2(t)=(2+6t)/5")
lambda_m <- 5L                  # Poisson mean of the visit counts
m_cap    <- 15L                 # visit counts truncated to [2, m_cap]
##       (P[Poisson(5) > 15] < 1e-4), so the
##       bounded-visits assumption holds exactly

## ---- Subject-specific process  X_i -------------------------------------
## X_i(u,t) = sum_{l,l'} C[l,l'] cos(l*pi*u) cos(l'*pi*t), with
## C[l,l'] ~ N(0,1) / (l^k_decay * l'^k_prime_decay) (decay = smoothness).
L             <- 100L           # number of modes
k_decay       <- 2L             # decay exponent: l^k
k_prime_decay <- 2L             # decay exponent: l'^{k'}

## ---- Noise  eps_ij(u) = tau(T_ij) * eta_ij(u) ---------------------------
tau     <- 0.5                  # baseline noise standard deviation
TAU_FUN <- function(t) rep(tau, length(t))   ## baseline scale function
##  tau(.) == tau; the robustness study (E8)
##  passes genuinely t-dependent versions
L_prime <- 50L                  # number of terms in eta_ij

## ---- Truncation and local-linear bandwidth grids ------------------------
Kmax     <- 20L                 # largest truncation level K ever considered
LL_HGRID <- exp(seq(log(0.03), log(0.35), length.out = 30))
# logarithmic h_t grid (sparse t-direction)
# over which the bivariate local-linear
# benchmark is given its oracle
LL_HGRID_U <- c(0, exp(seq(log(0.03), log(0.30), length.out = 6)))
# h_u grid, across the functional argument
# u; at h_u = 0 the joint fit reduces to
# the slice-wise (no u-smoothing) smoother
# of Yao-Mueller-Wang, nested as a section

## ---- U/T grids for the VISUALISATION panels only -----------------------
## (the per-experiment evaluation grids are in 'cfg' below, keyed Nu_eval/Nt_eval)
Nu      <- 1001L                # fine grid for observed curves
Nu_viz  <- 201L                 # coarser grid for surface plots (u-axis)
Nt_viz  <- 101L                 # grid for surface plots (t-axis)
u_fine <- seq(0, 1, length.out = Nu)
u_viz  <- seq(0, 1, length.out = Nu_viz)
t_viz  <- seq(0, 1, length.out = Nt_viz)

## ---- Output toggles ----------------------------------------------------
WRITE_PDF <- FALSE   # TRUE: also save each figure as fig_*.pdf (besides TikZ .tex)
WRITE_CSV <- FALSE   # TRUE: also write the res_*.csv files with the raw numbers
## (If tikzDevice is NOT installed, figures fall back to fig_*.pdf)

## ---- Run size, output folder, master seed ------------------------------
SETTINGS <- Sys.getenv("SIM_SETTINGS", unset = "FULL")  # "QUICK" | "FULL"
OUTDIR   <- Sys.getenv("SIM_OUTDIR",   unset = "./Téléchargements/Simulation")
BASESEED <- 2026L

# create the output folder if needed
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

## ---- TikZ availability -------------------------------------------------
## TikzDevice needs a working LaTeX install WITH the 'pgf' package (which
## provides TikZ). If 'pgf' is missing you get, we fall
## back to PDF figures so the run still completes.
USE_TIKZ <- requireNamespace("tikzDevice", quietly = TRUE)
if (USE_TIKZ) {
  suppressMessages(requireNamespace("tikzDevice"))
  options(tikzDefaultEngine = "pdftex")
  options(tikzMetricsDictionary = file.path(OUTDIR, "tikzMetrics"))   # cache the metrics
  # append the maths packages to whatever tikzDevice loads by default
  options(tikzLatexPackages = c(getOption("tikzLatexPackages"),
                                "\\usepackage{amsmath}\n", "\\usepackage{amssymb}\n"))
  tikz_works <- function() {                       # try one tiny metric computation
    f <- tempfile(fileext = ".tex")
    ok <- tryCatch({
      tikzDevice::tikz(f, width = 2, height = 2, standAlone = FALSE)
      plot.new(); text(0.5, 0.5, "metric test $x_1$"); TRUE
    }, error = function(e) FALSE)
    while (length(dev.list())) try(dev.off(), silent = TRUE)   # close any device left open
    isTRUE(ok)
  }
  if (!tikz_works()) {
    USE_TIKZ <- FALSE
    message("\n[!] tikzDevice cannot compute LaTeX metrics here, so figures will be\n",
            "    written as PDF instead of TikZ. To enable TikZ output, install the\n",
            "    LaTeX 'pgf' package, then re-run:\n",
            "      TeX Live : tlmgr install pgf preview\n",
            "      tinytex  : tinytex::tlmgr_install(c('pgf','preview'))\n",
            "    (Your final LaTeX document also needs \\usepackage{tikz}, which pgf provides.)\n")
  }
}

# Effort levels. Nu_eval/Nt_eval are the EVALUATION grids used to compute the
# integrated risk (distinct from the visualisation grids Nu/Nu_viz/Nt_viz above).
if (SETTINGS == "FULL") {
  cfg <- list(n_rate=c(100,200,500,1000,2000), R_rate=200L, Nu_eval=50L, Nt_eval=50L,
              n_design=c(125,250,500,1000), R_design=150L, K_design=8L,
              n_riskK=500L, R_riskK=100L, R_case12=80L, n_case12=150L,
              n_plug=c(50,100,200,500), R_plug=100L, R_unbias=2000L,
              n_cv=c(100,500), R_cv=100L,                       
              n_rob=200L, R_rob=150L)                           
} else {  # QUICK
  cfg <- list(n_rate=c(50,100,200), R_rate=15L, Nu_eval=30L, Nt_eval=30L,
              n_design=c(100,200,400), R_design=40L, K_design=6L,
              n_riskK=200L, R_riskK=10L, R_case12=20L, n_case12=120L,
              n_plug=c(80,160), R_plug=15L, R_unbias=300L,
              n_cv=c(80), R_cv=10L,                             
              n_rob=80L, R_rob=12L)                           
}
cat(sprintf("==== SETTINGS=%s | tikzDevice=%s | WRITE_PDF=%s | WRITE_CSV=%s ====\n",
            SETTINGS, USE_TIKZ, WRITE_PDF, WRITE_CSV))


# =========================================================================
#  NUMERICAL CORE  --  the data-generating process and the estimators
# =========================================================================

## make_Xi(): returns A NEW FUNCTION representing one random draw of the
## subject process X. Calling make_Xi() again gives an independent subject.
## The returned function maps (u, t) to the matrix X(u[i], t[j]).
make_Xi <- function(L_modes = L, k = k_decay, kp = k_prime_decay) {
  # matrix of random coefficients, divided entrywise by the decay weights
  C <- matrix(rnorm(L_modes * L_modes), L_modes, L_modes) /
    outer(seq_len(L_modes), seq_len(L_modes), function(l, lp) l^k * lp^kp)
  function(u, t)
    outer(u, seq_len(L_modes), function(u, l)  cos(l  * pi * u)) %*% C %*%
    t(outer(t, seq_len(L_modes), function(t, lp) cos(lp * pi * t)))
}

## gen_noise(): one draw of the noise field eta on the points u, for 'nc'
## curves (columns). Returns a (length(u))-by-nc matrix.
## The multiplication by the scale tau(T_ij) happens in gen_data().
gen_noise <- function(u, nc, Lp = L_prime) {
  B  <- matrix(rnorm(Lp * nc), Lp, nc) / (seq_len(Lp)^4)  # decaying random coefficients
  Cu <- outer(u, seq_len(Lp), function(u, l) cos(l * pi * u))
  v  <- as.vector((Cu^2) %*% (seq_len(Lp)^(-8)))          # pointwise variance v(u) > 0
  (Cu %*% B) / sqrt(v)                                    # unit variance at every u
}

## sample_T(): draw m visit times with density des$g on [0,1] by rejection
## sampling: propose uniform candidates, keep each with probability
## g(candidate)/gmax; the kept points then follow g exactly.
sample_T <- function(m, des = DES) {
  out <- numeric(0)
  while (length(out) < m) {
    x <- runif(2 * m)                       # uniform candidate times
    acc <- runif(2 * m)                     # uniform acceptance variables
    out <- c(out, x[acc <= des$g(x) / des$gmax])
  }
  sort(out[1:m])                            # return the first m, in increasing order
}

## phi_mat(): evaluate the first K Fourier-cosine basis functions at the
## points t (an orthonormal basis on [0,1]): phi_1 = 1 and
## phi_k(t) = sqrt(2) cos((k-1) pi t) for k >= 2. Returns length(t)-by-K.
phi_mat <- function(t, K) {
  Pm <- matrix(1, length(t), K)             # first column is phi_1 == 1
  if (K >= 2) for (k in 2:K) Pm[, k] <- sqrt(2) * cos((k - 1) * pi * t)
  Pm
}

## beta_true(): the true coefficients beta_k(u) = \int_0^1 mu(u,t) phi_k(t) dt.
## Computed by composite Simpson quadrature on a fine grid.
beta_true <- function(u, K) {
  Ns <- 2000L                                    # even: 2*Ns+1 Simpson nodes
  tg <- seq(0, 1, length.out = 2L * Ns + 1L)
  wq <- rep(c(2, 4), length.out = length(tg))    # Simpson pattern 1,4,2,4,...,4,1
  wq[1] <- 1; wq[length(wq)] <- 1
  wq <- wq / (3 * (2 * Ns))                      # times h/3, with h = 1/(2 Ns)
  MU_FUN(u, tg) %*% (phi_mat(tg, K) * wq)
}

## h_of_M(): the spacing window half-width h = ceil(1/2 + 1/2 * ln ln(M+20)).
## It grows so slowly that h = 2 for every realistic sample size.
h_of_M <- function(M) ceiling(0.5 + 0.5 * log(log(M + 20)))

## cn_weights(): the control-neighbors leave-one-out Voronoi weights
## omega_ij^MC = (1 + c_hat_ij - d_hat_ij)/M, computed exactly in 1-D from the
## sorted design ('s' = visit times, 'G' = c.d.f.). d_hat counts how many other
## points have T_ij as nearest neighbour; c_hat is the g-volume they vacate.
cn_weights <- function(s, G) {
  M <- length(s); ord <- order(s); ss <- s[ord]    # ss = sorted times
  mid <- (ss[-M] + ss[-1]) / 2                     # midpoints between neighbours
  Lb <- c(0, mid); Rb <- c(mid, 1)                 # left/right cell boundaries
  Vfull <- G(Rb) - G(Lb)                           # full Voronoi volume of each point
  nn <- integer(M); if (M >= 2) { nn[1] <- 2L; nn[M] <- (M - 1L) }
  if (M >= 3) {
    lg <- ss[2:(M-1)] - ss[1:(M-2)]; rg <- ss[3:M] - ss[2:(M-1)]
    nn[2:(M-1)] <- ifelse(lg <= rg, 1:(M-2), 3:M)  # nearest neighbour of each point
  }
  dh <- tabulate(nn, nbins = M)                    # d_hat: in-degree as nearest neighbour
  newL <- rep(NA_real_, M); if (M >= 2) { newL[2] <- 0; if (M >= 3) newL[3:M] <- (ss[1:(M-2)] + ss[3:M]) / 2 }
  Vl <- ifelse(seq_len(M) > 1, G(Rb) - G(newL), 0)
  newR <- rep(NA_real_, M); if (M >= 2) { newR[M-1] <- 1; if (M >= 3) newR[1:(M-2)] <- (ss[1:(M-2)] + ss[3:M]) / 2 }
  Vr <- ifelse(seq_len(M) < M, G(newR) - G(Lb), 0)
  nna <- (M - 1) - (seq_len(M) > 1) - (seq_len(M) < M)   # # points NOT adjacent to this one
  o <- numeric(M)
  o[ord] <- (1 + (nna * Vfull + Vl + Vr) - dh) / M       # put weights back in original order
  o
}

## beta_weighted(): generic coefficient estimator for a weight vector 'omega':
## beta_hat_k(u) = sum_ij omega_ij Y_ij(u) phi_k(T_ij)/g(T_ij). 'gfun' is the
## density used (true g, or a plug-in estimate). Returns nu-by-K.
beta_weighted <- function(Y, Tp, omega, gfun, K)
  Y %*% (omega * (phi_mat(Tp, K) / gfun(Tp)))

## beta_spacing(): the spacing estimator beta_hat_k(u) = sum_l Y_(l)(u) D_kl,
## with D_kl = (1/2h) \int_{T_(l-h)}^{T_(l+h)} phi_k. It uses NO density g.
## Window endpoints past the ends are obtained by symmetric reflection.
beta_spacing <- function(Y, Tp, K, h) {
  M <- length(Tp); ord <- order(Tp); Ts <- Tp[ord]; Yo <- Y[, ord, drop = FALSE]
  Text <- function(r) {                          # reflected order statistic at rank r
    out <- numeric(length(r)); lo <- r <= 0; hi <- r > M; in_ <- !lo & !hi
    out[in_] <- Ts[r[in_]]
    out[lo]  <- 2 * Ts[1] - Ts[2 - r[lo]]        # reflect at the left endpoint
    out[hi]  <- 2 * Ts[M] - Ts[2 * M - r[hi]]    # reflect at the right endpoint
    out
  }
  a <- Text((1:M) - h); bnd <- Text((1:M) + h)   # window endpoints for every rank
  D <- matrix(0, M, K)
  D[, 1] <- (bnd - a) / (2 * h)                  # phi_1 == 1: the integral is the width
  if (K >= 2) for (k in 2:K)                     # closed form of the cosine integral
    D[, k] <- (1 / (2 * h)) * sqrt(2) / ((k - 1) * pi) *
    (sin((k - 1) * pi * bnd) - sin((k - 1) * pi * a))
  Yo %*% D
}

## kde_g(): a FEASIBLE plug-in estimate of g from the observed times (to show
## the cost of not knowing g). Kernel density with reflection at 0 and 1,
## renormalised to integrate to 1, floored away from 0. Returns a function.
kde_g <- function(Tp) {
  z <- c(-Tp, Tp, 2 - Tp)                         # reflect the sample across both ends
  d <- density(z, from = 0, to = 1, n = 512)      # kernel density estimate on [0,1]
  y <- pmax(d$y, 1e-3); y <- y / mean(y)          # floor and renormalise
  gf <- approxfun(d$x, y, rule = 2)               # turn the estimate into a function
  function(x) pmax(gf(x), 0.05)
}

## mu_hat_grid(): assemble the fitted surface mu_hat(u,t) on a t-grid from the
## coefficient matrix 'beta', using the first K terms. Returns nu-by-(len tgrid).
mu_hat_grid <- function(beta, tgrid, K) beta[, 1:K, drop = FALSE] %*% t(phi_mat(tgrid, K))

## gen_data(): simulate ONE data set of n subjects. Returns a list with the
## pooled visit times T, the (nu) x M matrix Y of observed curves, the subject
## label of each observation, and n. Options:
##   signal_only = TRUE -> drop X and noise (used to isolate the design term);
##   unbalanced  = TRUE -> a fraction p_big of subjects get m_big visits, the
##                         rest get 2 (used to contrast Case 1 vs Case 2);
##   m_fixed = m        -> every subject gets exactly m visits
##                         (balanced arm of E4: Case 1 == Case 2 exactly);
##   des, tau_fun       -> design density and noise-scale function
##                         (the robustness study E8 swaps these in).
gen_data <- function(n, ug, signal_only = FALSE, unbalanced = FALSE,
                     p_big = 0.15, m_big = 20L, m_fixed = NULL,
                     des = DES, tau_fun = TAU_FUN, muf = MU_FUN) {
  Tp <- vector("list", n); Yc <- vector("list", n); subj <- integer(0)
  for (i in seq_len(n)) {
    if (!is.null(m_fixed)) {
      mi <- as.integer(m_fixed)
    } else if (unbalanced) {
      mi <- if (runif(1) < p_big) m_big else 2L
    } else {
      mi <- 0L                                            ## Poisson visits
      while (mi < 2L || mi > m_cap) mi <- rpois(1L, lambda_m)  # truncated to [2, m_cap]
    }
    Ti <- sample_T(mi, des)
    Yi <- muf(ug, Ti)
    if (!signal_only)                                     ## tau(T_ij) scaling
      Yi <- Yi + make_Xi()(ug, Ti) +
      gen_noise(ug, mi) * matrix(tau_fun(Ti), length(ug), mi, byrow = TRUE)
    Tp[[i]] <- Ti; Yc[[i]] <- Yi; subj <- c(subj, rep(i, mi))
  }
  list(T = unlist(Tp), Y = do.call(cbind, Yc), subj = subj, n = n)
}

## make_weights(): the three deterministic/random weight vectors for a data set.
make_weights <- function(d, des = DES) {
  M  <- length(d$T)
  nm <- table(d$subj)[as.character(d$subj)]         # m_i repeated for each observation
  list(unif = rep(1 / M, M),                        # Case 1
       bal  = 1 / (d$n * as.numeric(nm)),           # Case 2
       MC   = cn_weights(d$T, des$G))               # control-neighbors
}

## ise_curve(): integrated squared error of a fit as a function of K = 1..Kup.
## 'MUtrue' is mu on the (ug, tg) grid. Returns a vector of length Kup.
ise_curve <- function(beta, tg, MUtrue, Kup = Kmax)
  sapply(1:Kup, function(K) mean((mu_hat_grid(beta, tg, K) - MUtrue)^2))

## ---- the BIVARIATE LOCAL-LINEAR external benchmark ----------------- 
## ll2_ise(): ISE of the benchmark over the FULL bandwidth grid. Returns a
## length(hgrid_t) x length(hgrid_u) matrix; the oracle pair (h_t*, h_u*)
## is chosen downstream from the Monte Carlo average, exactly as the
## oracle K* of the series estimators.
ll2_ise <- function(Y, Tp, ug, tg, MUtrue,
                    hgrid_t = LL_HGRID, hgrid_u = LL_HGRID_U) {
  M <- length(Tp); Nt <- length(tg); Nu <- length(ug)
  Dt <- outer(Tp, tg, "-")                          # M x Nt:  T_ij - t_s
  Du <- outer(ug, ug, "-")                          # Nu x Nu: u_r  - u_q
  ISE <- matrix(NA_real_, length(hgrid_t), length(hgrid_u))
  for (a in seq_along(hgrid_t)) {
    ht <- hgrid_t[a]
    Kt <- pmax(1 - (Dt / ht)^2, 0)                  # t-kernel weights (M x Nt)
    T0 <- colSums(Kt); T1 <- colSums(Kt * Dt); T2 <- colSums(Kt * Dt * Dt)
    A0 <- Y %*% Kt                                  # response t-moments,
    A1 <- Y %*% (Kt * Dt)                           # per slice (Nu x Nt each)
    for (b in seq_along(hgrid_u)) {
      hu <- hgrid_u[b]
      if (hu <= 0) {
        ## h_u = 0 section: univariate local-linear in t at each u_r,
        ## i.e. the pooled smoother of Yao-Mueller-Wang per slice.
        den1 <- T0 * T2 - T1 * T1
        F <- (A0 * rep(T2, each = Nu) - A1 * rep(T1, each = Nu)) /
          rep(pmax(den1, 1e-300), each = Nu)
        bad <- which(den1 <= 1e-10)                 # degenerate windows: NW
        if (length(bad))
          F[, bad] <- (A0 / rep(pmax(T0, 1e-300), each = Nu))[, bad]
      } else {
        ## exact joint plane fit of eq. (8), via factorised normal equations
        Ku <- pmax(1 - (Du / hu)^2, 0)              # u-kernel weights (Nu x Nu)
        C1 <- Ku * Du; C2 <- Ku * Du * Du
        U0 <- colSums(Ku); U1 <- colSums(C1); U2 <- colSums(C2)
        R00 <- crossprod(Ku, A0)                    # response (u,t)-moments,
        R10 <- crossprod(C1, A0)                    # one per target point
        R01 <- crossprod(Ku, A1)                    # (Nu x Nt each)
        S11 <- outer(U0, T0); S12 <- outer(U1, T0); S13 <- outer(U0, T1)
        S22 <- outer(U2, T0); S23 <- outer(U1, T1); S33 <- outer(U0, T2)
        dS <- S11 * (S22 * S33 - S23^2) -           # 3x3 Gram determinant
          S12 * (S12 * S33 - S23 * S13) +
          S13 * (S12 * S23 - S22 * S13)
        dM <- R00 * (S22 * S33 - S23^2) -           # Cramer: first column
          S12 * (R10 * S33 - S23 * R01) +       # replaced by responses
          S13 * (R10 * S23 - S22 * R01)
        F <- dM / pmax(dS, 1e-300)                  # a0-hat at every target
        bad <- which(dS <= 1e-10 * pmax(S11, 1)^3)  # degenerate windows: NW
        if (length(bad)) F[bad] <- (R00 / pmax(S11, 1e-300))[bad]
      }
      ISE[a, b] <- mean((F - MUtrue)^2)             # ISE on the eval grid
    }
  }
  ISE
}


# =========================================================================
#  OUTPUT HELPERS  --  figures (TikZ/PDF), tables (booktabs), formatting
# =========================================================================

## Plotting dictionaries: point symbol, line type and LaTeX label of each
## estimator (keys: unif, bal, MC, spac, plug, ll).
EST_PCH <- c(unif = 19, bal = 17, MC = 15, spac = 18, plug = 4, ll = 1)
EST_LTY <- c(unif = 1,  bal = 2,  MC = 3,  spac = 5,  plug = 4, ll = 6)
EST_TEX <- c(unif = "$\\widehat\\mu^{(1)}$",
             bal  = "$\\widehat\\mu^{(2)}$",
             MC   = "$\\widehat\\mu^{(\\mathrm{MC})}$",
             spac = "$\\widehat\\mu^{(\\mathrm{sp})}$",
             plug = "$\\widehat\\mu^{(1)}_{\\widehat g}$",
             ll   = "$\\widehat{\\mu}^{(\\mathrm{LL})}$")
PT_CEX <- 1.05; LN_LWD <- 1.5; FIT_LWD <- 0.9

## setup_panel(): open one axis system with light margins and a faint grid.
setup_panel <- function(xlab, ylab, xlim, ylim) {
  par(mar = c(4.0, 4.4, 1.0, 0.9), mgp = c(2.5, 0.8, 0), tcl = -0.3)
  plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
       las = 1, bty = "l", cex.axis = 0.9, cex.lab = 1.0)
  faint_grid()
}
faint_grid <- function() grid(col = "grey85", lty = 3, lwd = 0.6)

## pts_mark(): a larger marker used to flag the oracle minimiser on a curve.
pts_mark <- function(x, y, pch) points(x, y, pch = pch, cex = PT_CEX * 1.45, lwd = 1.4)

## legend_box(): a compact, frameless legend.
legend_box <- function(pos, labels, pch, lty, ...)
  legend(pos, legend = labels, pch = pch, lty = lty, lwd = LN_LWD,
         bty = "n", cex = 0.82, seg.len = 2.8, pt.cex = PT_CEX, ...)

## emit_fig(): draw the figure 'plotfun' into OUTDIR/<name>.tex as TikZ code,
## or -- if TikZ is unavailable -- into <name>.pdf plus a one-line stub
## <name>.tex that \includegraphics's the PDF, so the LaTeX document compiles
## either way. With WRITE_PDF = TRUE a PDF copy is always written too.
emit_fig <- function(name, width, height, plotfun) {
  if (USE_TIKZ) {
    tikzDevice::tikz(file.path(OUTDIR, paste0(name, ".tex")),
                     width = width, height = height, standAlone = FALSE)
    plotfun(); dev.off()
    if (WRITE_PDF) {
      pdf(file.path(OUTDIR, paste0(name, ".pdf")), width = width, height = height)
      plotfun(); dev.off()
    }
  } else {
    pdf(file.path(OUTDIR, paste0(name, ".pdf")), width = width, height = height)
    plotfun(); dev.off()
    writeLines(sprintf("\\includegraphics[width=.92\\linewidth]{%s.pdf}", name),
               file.path(OUTDIR, paste0(name, ".tex")))
  }
  cat(sprintf("  wrote %s.tex (%s)\n", name, if (USE_TIKZ) "TikZ" else "PDF stub"))
}

## save_csv(): raw numbers behind a table/figure (only if WRITE_CSV = TRUE).
save_csv <- function(name, df) if (WRITE_CSV) {
  write.csv(df, file.path(OUTDIR, paste0(name, ".csv")), row.names = FALSE)
  cat(sprintf("  wrote %s.csv\n", name))
}

## latex_table(): one booktabs table float per file. 'headers' is the header
## row (character vector), 'body' a CHARACTER matrix of formatted cells.
latex_table <- function(name, headers, body, caption, label, align) {
  body <- as.matrix(body)
  lines <- c("\\begin{table}[H]", "\\centering",
             sprintf("\\begin{tabular}{%s}", align), "\\toprule",
             paste0(paste(headers, collapse = " & "), " \\\\"), "\\midrule",
             paste0(apply(body, 1, paste, collapse = " & "), " \\\\"),
             "\\bottomrule", "\\end{tabular}",
             sprintf("\\caption{%s}", caption),
             sprintf("\\label{%s}", label), "\\end{table}")
  writeLines(lines, file.path(OUTDIR, paste0(name, ".tex")))
  cat(sprintf("  wrote %s.tex\n", name))
}

## ---- number formatting ---------------------------------------------------
fmt_int <- function(v) sprintf("%d", as.integer(round(v)))
fmt_num <- function(v, d = 3) sprintf(paste0("%.", d, "f"), v)

## fmt_sci(): plain scientific notation $m.mm \times 10^{e}$.
fmt_sci <- function(v) {
  out <- character(length(v))
  for (i in seq_along(v)) {
    if (!is.finite(v[i])) { out[i] <- "--"; next }
    if (v[i] == 0)        { out[i] <- "$0$"; next }
    e <- floor(log10(abs(v[i])))
    out[i] <- sprintf("$%.2f\\times 10^{%d}$", v[i] / 10^e, e)
  }
  out
}

## fmt_val_se(): the TABLE CONVENTION of the paper -- "mean (s.e.) in units
## of the common power of ten": $m.mm\,(s.ss)\times 10^{e}$, or plain
## $m.mm\,(s.ss)$ when no exponent is needed. Decimals are widened (up to 4)
## so the s.e. always shows at least one significant digit.
fmt_val_se <- function(v, se, d = 2) {
  out <- character(length(v))
  for (i in seq_along(v)) {
    if (!is.finite(v[i])) { out[i] <- "--"; next }
    av <- abs(v[i])
    if (av == 0) { out[i] <- "$0$"; next }
    if (av >= 0.1 && av < 1000) {                       # no exponent needed
      dd <- d
      if (is.finite(se[i]) && se[i] > 0)
        dd <- min(max(d, ceiling(-log10(se[i]))), 4L)
      out[i] <- sprintf(paste0("$%.", dd, "f\\,(%.", dd, "f)$"), v[i], se[i])
    } else {                                            # common power of ten
      e <- floor(log10(av)); f <- 10^e
      dd <- d
      if (is.finite(se[i]) && se[i] > 0)
        dd <- min(max(d, ceiling(-log10(se[i] / f))), 4L)
      out[i] <- sprintf(paste0("$%.", dd, "f\\,(%.", dd, "f)\\times 10^{%d}$"),
                        v[i] / f, se[i] / f, e)
    }
  }
  out
}

## slope_se(): least-squares slope of y on x WITH its regression standard
## error -- every fitted log-log slope in the paper is reported this way.
slope_se <- function(x, y) {
  s <- summary(lm(y ~ x))$coefficients
  c(slope = unname(s[2, 1]), se = unname(s[2, 2]))
}
fmt_slope_se <- function(b, s) sprintf("$%.2f\\,(%.2f)$", b, s)

## ratio_se(): mean(x)/mean(y) for PAIRED samples, with the delta-method
## standard error (the covariance term captures the pairing).
ratio_se <- function(x, y) {
  R <- length(x); mx <- mean(x); my <- mean(y)
  va <- (var(x) / my^2 + mx^2 * var(y) / my^4 - 2 * mx * cov(x, y) / my^3) / R
  c(ratio = mx / my, se = sqrt(max(va, 0)))
}


# =========================================================================
#  E1 -- TOTAL RISK vs M: convergence rate, oracle K*, and the
#        local-linear benchmark (Tables tab_rate, tab_Kstar; fig_rate_total)
# =========================================================================
exp_rate <- function() {
  cat("\n[E1] total risk vs M: rate, oracle K*/h*, LL benchmark ...\n")
  set.seed(BASESEED)
  ug <- seq(0, 1, length.out = cfg$Nu_eval)
  tg <- seq(0, 1, length.out = cfg$Nt_eval)
  MU <- MU_FUN(ug, tg)
  keys  <- c("unif", "bal", "MC", "spac")      # the four series estimators
  keys5 <- c(keys, "ll")                       # ... plus the LL benchmark
  res <- NULL; Kres <- NULL
  ISbn <- list(); LLbn <- list()               # per-n ISE blocks (for slope bootstrap)
  for (n in cfg$n_rate) {
    R  <- cfg$R_rate
    IS <- setNames(lapply(keys, function(k) matrix(NA_real_, R, Kmax)), keys)
    LL <- array(NA_real_, c(R, length(LL_HGRID), length(LL_HGRID_U)))
    Mv <- numeric(R)
    for (r in seq_len(R)) {
      d <- gen_data(n, ug); Mv[r] <- length(d$T)
      w <- make_weights(d)
      B <- list(unif = beta_weighted(d$Y, d$T, w$unif, G_DENSITY, Kmax),
                bal  = beta_weighted(d$Y, d$T, w$bal,  G_DENSITY, Kmax),
                MC   = beta_weighted(d$Y, d$T, w$MC,   G_DENSITY, Kmax),
                spac = beta_spacing(d$Y, d$T, Kmax, h_of_M(length(d$T))))
      for (k in keys) IS[[k]][r, ] <- ise_curve(B[[k]], tg, MU)
      LL[r, , ] <- ll2_ise(d$Y, d$T, ug, tg, MU)           ## benchmark
    }
    Mbar <- mean(Mv)
    row <- list(n = n, M = Mbar); krow <- list(n = n, M = Mbar)
    for (k in keys) {                          # oracle K* of each series scheme
      mc <- colMeans(IS[[k]]); Ks <- which.min(mc)
      row[[paste0("IMSE_", k)]] <- mc[Ks]
      row[[paste0("se_",   k)]] <- sd(IS[[k]][, Ks]) / sqrt(R)
      krow[[paste0("Kstar_", k)]] <- Ks
    }
    mcLL <- apply(LL, c(2, 3), mean)                       # mean ISE per (h_t, h_u)
    iLL  <- which(mcLL == min(mcLL), arr.ind = TRUE)[1, ]  # joint oracle pair
    row$IMSE_ll <- mcLL[iLL[1], iLL[2]]
    row$se_ll   <- sd(LL[, iLL[1], iLL[2]]) / sqrt(R)
    krow$ht_ll  <- LL_HGRID[iLL[1]]                        # oracle h_t (sparse dir.)
    krow$hu_ll  <- LL_HGRID_U[iLL[2]]                      # oracle h_u (0 = slice-wise)
    res  <- rbind(res,  as.data.frame(row))
    Kres <- rbind(Kres, as.data.frame(krow))
    ISbn[[length(ISbn) + 1L]] <- IS        # keep blocks for the slope bootstrap
    LLbn[[length(LLbn) + 1L]] <- LL
    cat(sprintf("  n=%5d (M~%5.0f) IMSE %s | LL %.3g  K* %s  h_t*=%.3f h_u*=%.3f\n",
                n, Mbar,
                paste(sprintf("%.3g", unlist(row[paste0("IMSE_", keys)])), collapse = " "),
                row$IMSE_ll,
                paste(unlist(krow[paste0("Kstar_", keys)]), collapse = "/"),
                krow$ht_ll, krow$hu_ll))
  }
  save_csv("res_rate_total", res); save_csv("res_Kstar", Kres)
  
  ## fitted log-log slopes. The POINT estimates come from least squares on the
  ## design-point means; the STANDARD ERRORS are obtained by a nonparametric
  ## bootstrap over the R replications (see below), NOT from the regression. The
  ## regression s.e. would only measure scatter of the means about the line and
  ## badly understates the true uncertainty, which is dominated by Monte Carlo
  ## noise in each mean and by the integer/grid oscillation of the oracle tuning
  ## ($K^\star$ jumps in unit steps, $h_t^\star$ on a grid). The bootstrap
  ## resamples the replications at each design point, recomputes the oracle-tuned
  ## summary, and refits the slope; its s.d. is the honest s.e.
  lM  <- log10(res$M)
  sl  <- sapply(keys5, function(k) slope_se(lM, log10(res[[paste0("IMSE_", k)]])))
  ksl <- sapply(keys,  function(k) slope_se(lM, log10(Kres[[paste0("Kstar_", k)]])))
  hsl <- slope_se(lM, log10(Kres$ht_ll))      # prediction -1/5 concerns h_t;
  # h_u has no rate prediction (it is
  # expected to sit at 0, see the tex)
  ## ---- honest bootstrap standard errors of all fitted slopes --------------
  set.seed(BASESEED + 100L)                   # reproducible, after all data gen
  Bb <- 1000L; nN <- length(ISbn)
  xc <- lM - mean(lM); den <- sum(xc^2)       # closed-form LS slope = (xc . y)/den
  ls_slope <- function(y) sum(xc * y) / den
  nht <- length(LL_HGRID)                     # flatten each LL block to R x (nht*nhu)
  LLflat <- lapply(LLbn, function(A) matrix(A, nrow = dim(A)[1]))
  htidx  <- rep(seq_len(nht), times = length(LL_HGRID_U))   # h_t index of each column
  boot_sl <- matrix(NA_real_, Bb, length(keys5), dimnames = list(NULL, keys5))
  boot_k  <- matrix(NA_real_, Bb, length(keys),  dimnames = list(NULL, keys))
  boot_ht <- numeric(Bb)
  for (b in seq_len(Bb)) {
    yR <- matrix(NA_real_, nN, length(keys5)); yK <- matrix(NA_real_, nN, length(keys))
    yH <- numeric(nN)
    for (i in seq_len(nN)) {
      Ri <- nrow(ISbn[[i]][[1]]); id <- sample.int(Ri, Ri, replace = TRUE)
      for (j in seq_along(keys)) {
        mc <- colMeans(ISbn[[i]][[keys[j]]][id, , drop = FALSE])
        yR[i, j] <- log10(min(mc)); yK[i, j] <- log10(which.min(mc))
      }
      mcLL <- colMeans(LLflat[[i]][id, , drop = FALSE])   # mean ISE per (h_t,h_u) cell
      jmin <- which.min(mcLL)
      yR[i, length(keys5)] <- log10(mcLL[jmin])
      yH[i] <- log10(LL_HGRID[htidx[jmin]])
    }
    for (j in seq_along(keys5)) boot_sl[b, j] <- ls_slope(yR[, j])
    for (j in seq_along(keys))  boot_k[b, j]  <- ls_slope(yK[, j])
    boot_ht[b] <- ls_slope(yH)
  }
  sl["se", ]  <- apply(boot_sl, 2, sd)        # overwrite with bootstrap s.e.
  ksl["se", ] <- apply(boot_k,  2, sd)
  hsl["se"]   <- sd(boot_ht)
  
  ## ---- tab_rate: risks with s.e., slope row at the bottom -----------------
  ## All risks are reported in units of 10^{-3} (a single factor for the whole
  ## table, stated in the caption), so each cell is the compact "value (s.e.)"
  ## without a per-cell power of ten -- this keeps the five-estimator table narrow.
  sc <- 1e3
  fmt_r <- function(v, se) sprintf("$%.2f\\,(%.2f)$", v * sc, se * sc)
  body <- character(0)
  for (i in seq_len(nrow(res)))
    body <- c(body, paste(c(fmt_int(res$n[i]), fmt_int(res$M[i]),
                            fmt_r(unlist(res[i, paste0("IMSE_", keys5)]),
                                  unlist(res[i, paste0("se_",   keys5)]))),
                          collapse = " & "))
  writeLines(c("\\begin{table}[H]", "\\centering",
               "\\begin{tabular}{rrccccc}", "\\toprule",
               paste0("$n$ & $\\overline M$ & ", paste(EST_TEX[keys5], collapse = " & "), " \\\\"),
               "\\midrule", paste0(body, " \\\\"), "\\midrule",
               paste0("\\multicolumn{2}{r}{slope (s.e.)} & ",
                      paste(fmt_slope_se(sl["slope", ], sl["se", ]), collapse = " & "), " \\\\"),
               "\\bottomrule", "\\end{tabular}",
               sprintf(paste0("\\caption{Total integrated risk $\\overline{\\mathrm{IMSE}}$ at the ",
                              "oracle tuning (truncation $K^\\star$ for the series estimators, bandwidth ",
                              "pair $(h_t^\\star,h_u^\\star)$ for the bivariate local-linear benchmark) ",
                              "versus the pooled sample size ($R=%d$ replications ",
                              "per design point). All risks are in units of $10^{-3}$, with MC s.e.\\ in ",
                              "parentheses. Bottom row: fitted slope of $\\log_{10}\\overline{\\mathrm{IMSE}}$ ",
                              "on $\\log_{10}M$ with bootstrap s.e.\\ (resampling the $R$ replications, ",
                              "so it reflects Monte Carlo and oracle-tuning uncertainty).}"),
                       cfg$R_rate),
               "\\label{tab:rate}", "\\end{table}"),
             file.path(OUTDIR, "tab_rate.tex"))
  cat("  wrote tab_rate.tex\n")
  
  ## ---- tab_Kstar: oracle K* (4 schemes) and oracle h_t* (LL) --------------
  kbody <- character(0)
  for (i in seq_len(nrow(Kres)))
    kbody <- c(kbody, paste(c(fmt_int(Kres$n[i]), fmt_int(Kres$M[i]),
                              fmt_int(unlist(Kres[i, paste0("Kstar_", keys)])),
                              fmt_num(Kres$ht_ll[i], 3)), collapse = " & "))
  writeLines(c("\\begin{table}[H]", "\\centering",
               "\\begin{tabular}{rrccccc}", "\\toprule",
               paste0(" & & \\multicolumn{4}{c}{oracle truncation $K^\\star$} & ",
                      "LL \\\\"),
               "\\cmidrule(lr){3-6}\\cmidrule(lr){7-7}",
               paste0("$n$ & $\\overline M$ & ", paste(EST_TEX[keys], collapse = " & "),
                      " & $h_t^\\star$ \\\\"),
               "\\midrule", paste0(kbody, " \\\\"), "\\midrule",
               paste0("\\multicolumn{2}{r}{slope (s.e.)} & ",
                      paste(fmt_slope_se(ksl["slope", ], ksl["se", ]), collapse = " & "),
                      " & ", fmt_slope_se(hsl["slope"], hsl["se"]), " \\\\"),
               "\\bottomrule", "\\end{tabular}",
               sprintf(paste0("\\caption{Oracle truncation $K^\\star$ of the four series ",
                              "estimators and oracle $t$-bandwidth $h_t^\\star$ of the bivariate ",
                              "local-linear benchmark versus $M$ (same runs as Table~\\ref{tab:rate}, ",
                              "$R=%d$). The $u$-bandwidth is omitted: $h_u^\\star=0$ at every $n$ (the ",
                              "$h_u$-grid contains $0$, the slice-wise section), so no $u$-smoothing is ",
                              "selected. Bottom row: fitted slope of $\\log_{10}K^\\star$ (resp.\\ ",
                              "$\\log_{10}h_t^\\star$) on $\\log_{10}M$ with bootstrap s.e.}"), cfg$R_rate),
               "\\label{tab:Kstar}", "\\end{table}"),
             file.path(OUTDIR, "tab_Kstar.tex"))
  cat("  wrote tab_Kstar.tex\n")
  
  ## ---- fig_rate_total: 5 log-log curves, fitted lines, slopes in legend ---
  emit_fig("fig_rate_total", 6.4, 5.8, function() {
    x  <- log10(res$M)
    ys <- lapply(keys5, function(k) log10(res[[paste0("IMSE_", k)]]))
    setup_panel("$\\log_{10} M$", "$\\log_{10}\\overline{\\mathrm{IMSE}}$",
                range(x), range(unlist(ys)) + c(-0.05, 0.45))
    for (j in seq_along(keys5)) {
      k <- keys5[j]
      lines(x, ys[[j]], lty = EST_LTY[k], lwd = LN_LWD)
      points(x, ys[[j]], pch = EST_PCH[k], cex = PT_CEX)
    }
    legend_box("topright", sprintf("%s  ($%.2f$)", EST_TEX[keys5], sl["slope", keys5]),
               EST_PCH[keys5], EST_LTY[keys5])
  })
  invisible(list(res = res, Kres = Kres, slopes = sl, kslopes = ksl, hslope = hsl))
}



# =========================================================================
#  E2 -- MONTE CARLO DIAGNOSTIC for unbiasedness of the coefficient
#        estimators (Table tab_unbias)
# =========================================================================
exp_unbias <- function() {
  cat("\n[E2] unbiasedness of beta_hat_k(u) on signal-only data ...\n")
  set.seed(BASESEED + 2L)
  ug <- seq(0, 1, length.out = 50); K <- 8L; n <- 30L; R <- cfg$R_unbias
  uidx <- c(10L, 25L, 40L)                    # u ~ 0.18, 0.49, 0.80
  keys <- c("unif", "bal", "MC", "spac")
  A  <- setNames(lapply(keys, function(k) matrix(0, length(uidx), K)), keys)
  A2 <- A                                     # running sums of B and B^2
  for (r in seq_len(R)) {
    d <- gen_data(n, ug, signal_only = TRUE)  # mu only: removes X and noise,
    # leaving just the design randomness
    w <- make_weights(d)
    B <- list(unif = beta_weighted(d$Y, d$T, w$unif, G_DENSITY, K),
              bal  = beta_weighted(d$Y, d$T, w$bal,  G_DENSITY, K),
              MC   = beta_weighted(d$Y, d$T, w$MC,   G_DENSITY, K),
              spac = beta_spacing(d$Y, d$T, K, h_of_M(length(d$T))))
    for (k in keys) {
      A[[k]]  <- A[[k]]  + B[[k]][uidx, ]
      A2[[k]] <- A2[[k]] + B[[k]][uidx, ]^2
    }
  }
  bt <- beta_true(ug, K)[uidx, ]
  d_cells <- length(uidx) * K                  # number of (u,k) grid cells
  res <- do.call(rbind, lapply(keys, function(k) {
    m   <- A[[k]] / R                                      # MC mean of beta_hat
    se  <- sqrt(pmax(A2[[k]] / R - m^2, 0) * (R / (R - 1)) / R)  # its MC s.e.
    tval <- (m - bt) / se                                  # STUDENTISED bias:
    ## under exact unbiasedness each t_{u,k} ~ N(0,1), so the panel of |t| over
    ## the (u,k) grid is the clean diagnostic, REGARDLESS of R: max|t| stays at
    ## the N(0,1) extreme of a grid this size, and the proportion of cells with
    ## |t|>2 stays near P(|N(0,1)|>2) ~ 0.0455. A genuine bias instead drives
    ## max|t| (and the proportion) up without bound as R grows and s.e. -> 0.
    data.frame(est = k,
               mean_abs_bias = mean(abs(m - bt)),          # average over cells
               max_abs_bias  = max(abs(m - bt)),           # worst cell (scale)
               max_abs_t     = max(abs(tval)),             # the decisive number
               frac_t_gt2    = mean(abs(tval) > 2))        # d^{-1} #{|t|>2}
  }))
  save_csv("res_unbias", res)
  latex_table("tab_unbias",
              headers = c("estimator",
                          "$\\operatorname{avg}_{u,k}\\bigl|\\widehat{\\mathbb{E}}\\,\\widehat{\\beta}_k(u)-\\beta_k(u)\\bigr|$",
                          "$\\max_{u,k}\\bigl|\\widehat{\\mathbb{E}}\\,\\widehat{\\beta}_k(u)-\\beta_k(u)\\bigr|$",
                          "$\\max_{u,k}|t_{u,k}|$",
                          "$\\widehat p_2$"),
              body = cbind(EST_TEX[res$est], fmt_sci(res$mean_abs_bias),
                           fmt_sci(res$max_abs_bias), fmt_num(res$max_abs_t, 2),
                           fmt_num(res$frac_t_gt2, 3)),
              caption = sprintf(paste0("Diagnostic for the unbiasedness of ",
                                       "the coefficient estimators on signal-only data ",
                                       "($Y_{ij}=\\mu(\\cdot,T_{ij})$), over ",
                                       "$u\\in\\{0.18,\\,0.49,\\,0.80\\}$ and $k\\le%d$ ($d=%d$ cells; $n=%d$, ",
                                       "$R=%d$). Columns: mean and maximum absolute bias; the maximum ",
                                       "studentised bias $\\max_{u,k}|t_{u,k}|$; and the proportion ",
                                       "$\\widehat p_2$ of cells with $|t_{u,k}|>2$."),
                                K, d_cells, n, R),
              label = "tab:unbias", align = "lcccc")
  print(res, digits = 3)
  invisible(res)
}



# =========================================================================
#  E3 -- RISK vs TRUNCATION K, incl. the kernel plug-in
#        (Figure fig_riskK, Table tab_riskK)
# =========================================================================
exp_riskK <- function() {
  cat("\n[E3] risk vs truncation K (incl. kernel plug-in) ...\n")
  set.seed(BASESEED + 3L)
  n <- cfg$n_riskK; R <- cfg$R_riskK
  ug <- seq(0, 1, length.out = cfg$Nu_eval)
  tg <- seq(0, 1, length.out = cfg$Nt_eval)
  MU <- MU_FUN(ug, tg)
  keys <- c("unif", "bal", "MC", "spac", "plug")
  IS <- setNames(lapply(keys, function(k) matrix(NA_real_, R, Kmax)), keys)
  for (r in seq_len(R)) {
    d <- gen_data(n, ug); w <- make_weights(d)
    gh <- kde_g(d$T)                          # feasible plug-in estimate of g
    B <- list(unif = beta_weighted(d$Y, d$T, w$unif, G_DENSITY, Kmax),
              bal  = beta_weighted(d$Y, d$T, w$bal,  G_DENSITY, Kmax),
              MC   = beta_weighted(d$Y, d$T, w$MC,   G_DENSITY, Kmax),
              spac = beta_spacing(d$Y, d$T, Kmax, h_of_M(length(d$T))),
              plug = beta_weighted(d$Y, d$T, w$unif, gh, Kmax))
    for (k in keys) IS[[k]][r, ] <- ise_curve(B[[k]], tg, MU)
  }
  mc   <- sapply(keys, function(k) colMeans(IS[[k]]))      # Kmax x 5 matrix
  Kst  <- apply(mc, 2, which.min)
  Imin <- mc[cbind(Kst, seq_along(keys))]
  Ise  <- sapply(seq_along(keys), function(j) sd(IS[[keys[j]]][, Kst[j]]) / sqrt(R))
  res  <- data.frame(est = keys, Kstar = Kst, IMSE = Imin, se = Ise)
  save_csv("res_riskK", res)
  latex_table("tab_riskK",
              headers = c("estimator", "$K^\\star$",
                          "$\\overline{\\mathrm{IMSE}}$ at $K^\\star$ (s.e.)"),
              body = cbind(EST_TEX[keys], fmt_int(Kst), fmt_val_se(Imin, Ise)),
              caption = sprintf(paste0("Oracle truncation and minimum integrated risk for ",
                                       "the four estimators and the kernel plug-in $\\widehat\\mu^{(1)}_{\\widehat g}$ ",
                                       ", $n=%d$, $R=%d$. ",
                                       "MC s.e.\\ in parentheses, in units of the common power of ten."),
                                n, R),
              label = "tab:riskK", align = "lcc")
  emit_fig("fig_riskK", 6.4, 5.8, function() {
    Y <- log10(mc)
    setup_panel("truncation $K$", "$\\log_{10}\\overline{\\mathrm{IMSE}}$",
                c(1, Kmax), range(Y) + c(-0.05, 0.40))
    for (j in seq_along(keys)) {
      k <- keys[j]
      lines(1:Kmax, Y[, j], lty = EST_LTY[k], lwd = LN_LWD)
      points(1:Kmax, Y[, j], pch = EST_PCH[k], cex = 0.62)
      pts_mark(Kst[j], Y[Kst[j], j], EST_PCH[k])   # flag the oracle K*
    }
    legend_box("topright", EST_TEX[keys], EST_PCH[keys], EST_LTY[keys])
  })
  print(res, digits = 3)
  invisible(res)
}


# =========================================================================
#  E4 -- CASE 1 vs CASE 2 in balanced / unbalanced designs
#        (Figure fig_case12, Table tab_case12) -- PAIRED comparison
# =========================================================================
exp_case12 <- function() {
  cat("\n[E4] Case 1 vs Case 2: balanced vs unbalanced designs ...\n")
  n <- cfg$n_case12; R <- cfg$R_case12
  ug <- seq(0, 1, length.out = cfg$Nu_eval)
  tg <- seq(0, 1, length.out = cfg$Nt_eval)
  MU <- MU_FUN(ug, tg)
  one <- function(type, seed) {               # run one design arm
    set.seed(seed)
    I1 <- matrix(NA_real_, R, Kmax); I2 <- I1
    mA <- numeric(R); mH <- numeric(R)
    for (r in seq_len(R)) {
      d <- if (type == "balanced") gen_data(n, ug, m_fixed = 5L)   ## 
      else                    gen_data(n, ug, unbalanced = TRUE)
      w <- make_weights(d)
      if (type == "balanced" && r == 1L)      ## internal consistency:
        stopifnot(max(abs(w$unif - w$bal)) < 1e-15)  # 1/M == 1/(n*5) exactly
      B1 <- beta_weighted(d$Y, d$T, w$unif, G_DENSITY, Kmax)
      B2 <- beta_weighted(d$Y, d$T, w$bal,  G_DENSITY, Kmax)
      I1[r, ] <- ise_curve(B1, tg, MU); I2[r, ] <- ise_curve(B2, tg, MU)
      mi <- as.numeric(table(d$subj))
      mA[r] <- mean(mi); mH[r] <- 1 / mean(1 / mi)   # arithmetic/harmonic means
    }
    K1 <- which.min(colMeans(I1)); K2 <- which.min(colMeans(I2))
    x1 <- I1[, K1]; x2 <- I2[, K2]            # paired: same replications
    list(v1 = mean(x1), se1 = sd(x1) / sqrt(R),
         v2 = mean(x2), se2 = sd(x2) / sqrt(R),
         ratio = mean(x1) / mean(x2), se_diff = sd(x1 - x2) / sqrt(R),
         mA = mean(mA), mH = mean(mH))
  }
  bal <- one("balanced",   BASESEED + 4L)
  unb <- one("unbalanced", BASESEED + 5L)
  res <- data.frame(design = c("balanced", "unbalanced"),
                    mA = c(bal$mA, unb$mA), mH = c(bal$mH, unb$mH),
                    v1 = c(bal$v1, unb$v1), se1 = c(bal$se1, unb$se1),
                    v2 = c(bal$v2, unb$v2), se2 = c(bal$se2, unb$se2),
                    ratio = c(bal$ratio, unb$ratio),
                    se_diff = c(bal$se_diff, unb$se_diff))
  save_csv("res_case12", res)
  latex_table("tab_case12",
              headers = c("design", "$\\bar m_A$", "$\\bar m_H$",
                          EST_TEX["unif"], EST_TEX["bal"],
                          "ratio $\\widehat\\mu^{(1)}/\\widehat\\mu^{(2)}$",
                          "s.e.\\ of difference"),
              body = cbind(res$design, fmt_num(res$mA, 2), fmt_num(res$mH, 2),
                           fmt_val_se(res$v1, res$se1), fmt_val_se(res$v2, res$se2),
                           fmt_num(res$ratio, 3), fmt_sci(res$se_diff)),
              caption = sprintf(paste0("Integrated risk at the oracle truncation for ",
                                       "$\\widehat{\\mu}^{(1)}$ and $\\widehat{\\mu}^{(2)}$ under a balanced design ",
                                       "($m_i\\equiv5$) and an unbalanced two-point mixture ($15\\%%$ of subjects ",
                                       "with $m_i=20$, the rest $m_i=2$); $n=%d$, $R=%d$ paired replications. ",
                                       "$\\bar m_A$/$\\bar m_H$: arithmetic/harmonic mean visit counts. MC s.e.\\ ",
                                       "in parentheses; the last column is the MC s.e.\\ of the paired difference."), n, R),
              label = "tab:case12", align = "lcccccc")
  emit_fig("fig_case12", 6.0, 4.7, function() {
    H  <- rbind(c(bal$v1, unb$v1), c(bal$v2, unb$v2))
    SE <- rbind(c(bal$se1, unb$se1), c(bal$se2, unb$se2))
    colnames(H) <- c("balanced", "unbalanced")
    par(mar = c(2.6, 4.4, 1.0, 0.9), mgp = c(3.3, 0.8, 0))
    bp <- barplot(H, beside = TRUE, col = c("grey35", "grey80"),
                  ylab = "$\\overline{\\mathrm{IMSE}}$ at $K^\\star$",
                  ylim = c(0, 1.30 * max(H + 2 * SE)), las = 1)
    suppressWarnings(arrows(bp, H - 2 * SE, bp, H + 2 * SE,
                            angle = 90, code = 3, length = 0.04))
    text(colMeans(bp), apply(H + 2 * SE, 2, max) + 0.10 * max(H),
         sprintf("ratio $=%.3f$", c(bal$ratio, unb$ratio)), cex = 0.88)
    legend("topleft", c(EST_TEX["unif"], EST_TEX["bal"]),
           fill = c("grey35", "grey80"), bty = "n", cex = 0.9)
  })
  print(res[, c("design", "mA", "mH", "v1", "v2", "ratio")], digits = 4)
  invisible(list(res = res))
}


# =========================================================================
#  E5 -- DECAY OF THE DESIGN TERM (signal-only data)
#        (Figure fig_design, Table tab_design)
# =========================================================================
exp_design <- function() {
  cat("\n[E5] decay of the design term (signal-only data) ...\n")
  set.seed(BASESEED + 6L)
  ug <- seq(0, 1, length.out = 40); K <- cfg$K_design; R <- cfg$R_design
  keys <- c("unif", "spac", "MC")
  res <- NULL
  for (n in cfg$n_design) {
    ## On signal-only data beta_hat_k(u) = M_k(u): its variance ACROSS the
    ## replications, summed over k and integrated in u, IS the design term.
    s1 <- setNames(lapply(keys, function(k) matrix(0, length(ug), K)), keys)
    s2 <- s1; Mv <- numeric(R)                # running sums of B and B^2
    for (r in seq_len(R)) {
      d <- gen_data(n, ug, signal_only = TRUE); Mv[r] <- length(d$T)
      w <- make_weights(d)
      B <- list(unif = beta_weighted(d$Y, d$T, w$unif, G_DENSITY, K),
                spac = beta_spacing(d$Y, d$T, K, h_of_M(length(d$T))),
                MC   = beta_weighted(d$Y, d$T, w$MC,   G_DENSITY, K))
      for (k in keys) { s1[[k]] <- s1[[k]] + B[[k]]; s2[[k]] <- s2[[k]] + B[[k]]^2 }
    }
    row <- list(n = n, M = mean(Mv))
    for (k in keys) {
      V <- (s2[[k]] / R - (s1[[k]] / R)^2) * R / (R - 1)   # unbiased variance
      row[[paste0("D_", k)]] <- sum(colMeans(V))           # sum_k int Var(M_k(u)) du
    }
    res <- rbind(res, as.data.frame(row))
    cat(sprintf("  n=%5d (M~%5.0f) D: %s\n", n, row$M,
                paste(sprintf("%.3g", unlist(row[paste0("D_", keys)])), collapse = " ")))
  }
  save_csv("res_design", res)
  lM <- log10(res$M)
  sl <- sapply(keys, function(k) slope_se(lM, log10(res[[paste0("D_", k)]])))
  body <- character(0)
  for (i in seq_len(nrow(res)))
    body <- c(body, paste(c(fmt_int(res$M[i]), fmt_sci(res$D_unif[i]),
                            fmt_sci(res$D_spac[i]), fmt_sci(res$D_MC[i])),
                          collapse = " & "))
  writeLines(c("\\begin{table}[H]", "\\centering",
               "\\begin{tabular}{rccc}", "\\toprule",
               paste0("$M$ & ", EST_TEX["unif"], " (det.) & ", EST_TEX["spac"],
                      " (spacing) & ", EST_TEX["MC"], " (contr.-neighb.) \\\\"),
               "\\midrule", paste0(body, " \\\\"), "\\midrule",
               paste0("slope (s.e.) & ",
                      paste(fmt_slope_se(sl["slope", ], sl["se", ]), collapse = " & "), " \\\\"),
               "\\bottomrule", "\\end{tabular}",
               sprintf(paste0("\\caption{Design term $\\mathcal{R}^{\\mathrm{design}}=",
                              "\\sum_{k\\le%d}\\int_{\\mathcal{U}}\\mathrm{Var}(M_k(u))\\,\\mathrm{d}u$, ",
                              "isolated from signal-only data, versus the pooled size $M$ ($R=%d$ ",
                              "replications per design point). Bottom row: fitted log--log slope with ",
                              "regression s.e.}"),
                       K, R),
               "\\label{tab:design}", "\\end{table}"),
             file.path(OUTDIR, "tab_design.tex"))
  cat("  wrote tab_design.tex\n")
  emit_fig("fig_design", 6.4, 4.7, function() {
    x  <- log10(res$M)
    ys <- lapply(keys, function(k) log10(res[[paste0("D_", k)]]))
    setup_panel("$\\log_{10} M$", "$\\log_{10}\\mathcal{R}^{\\mathrm{design}}$",
                range(x), range(unlist(ys)) + c(-0.10, 0.55))
    for (j in seq_along(keys)) {
      k <- keys[j]
      abline(lm(ys[[j]] ~ x), lty = EST_LTY[k], lwd = FIT_LWD, col = "grey60")
      lines(x, ys[[j]], lty = EST_LTY[k], lwd = LN_LWD)
      points(x, ys[[j]], pch = EST_PCH[k], cex = PT_CEX)
    }
    legend_box("topright", sprintf("%s  ($%.2f$)", EST_TEX[keys], sl["slope", keys]),
               EST_PCH[keys], EST_LTY[keys])
  })
  invisible(list(res = res, slopes = sl))
}


# =========================================================================
#  E6 -- SUBJECT-LEVEL 5-FOLD CROSS-VALIDATED TRUNCATION (Table tab_cv) 
# =========================================================================
exp_cv <- function() {
  cat("\n[E6] subject-level 5-fold cross-validated truncation ...\n")
  set.seed(BASESEED + 8L)
  
  V <- 5L
  des <- DES
  ug <- seq(0, 1, length.out = cfg$Nu_eval)
  tg <- seq(0, 1, length.out = cfg$Nt_eval)
  MU <- MU_FUN(ug, tg)
  
  keys <- c("unif", "bal", "MC", "spac")
  res <- NULL
  
  for (n in cfg$n_cv) {
    R   <- cfg$R_cv
    IS  <- setNames(lapply(keys, function(k) matrix(NA_real_, R, Kmax)), keys)
    KCV <- setNames(lapply(keys, function(k) integer(R)), keys)
    
    for (r in seq_len(R)) {
      d <- gen_data(n, ug, des = des)
      w <- make_weights(d, des)
      
      ## Full-sample ISE curves, to translate K_CV into a realised risk
      B <- list(
        unif = beta_weighted(d$Y, d$T, w$unif, des$g, Kmax),
        bal  = beta_weighted(d$Y, d$T, w$bal,  des$g, Kmax),
        MC   = beta_weighted(d$Y, d$T, w$MC,   des$g, Kmax),
        spac = beta_spacing(d$Y, d$T, Kmax, h_of_M(length(d$T)))
      )
      
      for (k in keys) {
        IS[[k]][r, ] <- ise_curve(B[[k]], tg, MU)
      }
      
      ## CV criterion, accumulated fold by fold
      fold <- sample(rep_len(1:V, n))         # random subject-level folds
      sse  <- setNames(lapply(keys, function(k) numeric(Kmax)), keys)
      
      for (v in 1:V) {
        tr  <- fold[d$subj] != v              # TRUE = training observation
        
        Ttr <- d$T[tr]
        Ytr <- d$Y[, tr,  drop = FALSE]
        
        Tte <- d$T[!tr]
        Yte <- d$Y[, !tr, drop = FALSE]
        
        ntr  <- sum(fold != v)                # number of training subjects
        nmtr <- table(d$subj[tr])[as.character(d$subj[tr])]
        
        wtr <- list(
          unif = rep(1 / length(Ttr), length(Ttr)),
          bal  = 1 / (ntr * as.numeric(nmtr)),
          MC   = cn_weights(Ttr, des$G)
        )
        
        Btr <- list(
          unif = beta_weighted(Ytr, Ttr, wtr$unif, des$g, Kmax),
          bal  = beta_weighted(Ytr, Ttr, wtr$bal,  des$g, Kmax),
          MC   = beta_weighted(Ytr, Ttr, wtr$MC,   des$g, Kmax),
          spac = beta_spacing(Ytr, Ttr, Kmax, h_of_M(length(Ttr)))
        )
        
        Pte <- phi_mat(Tte, Kmax)
        
        for (k in keys) {
          Rk <- Yte
          
          for (K in 1:Kmax) {
            Rk <- Rk - Btr[[k]][, K, drop = FALSE] %*%
              t(Pte[, K, drop = FALSE])
            
            sse[[k]][K] <- sse[[k]][K] + sum(Rk * Rk)
          }
        }
      }
      
      for (k in keys) {
        KCV[[k]][r] <- which.min(sse[[k]])
      }
    }
    
    for (k in keys) {
      mc <- colMeans(IS[[k]])
      Ks <- which.min(mc)                     # oracle K*
      
      infl <- IS[[k]][cbind(seq_len(R), KCV[[k]])] / IS[[k]][, Ks]
      
      res <- rbind(
        res,
        data.frame(
          n = n,
          est = k,
          Kstar = Ks,
          K_med = unname(quantile(KCV[[k]], 0.50)),
          K_q1  = unname(quantile(KCV[[k]], 0.25)),
          K_q3  = unname(quantile(KCV[[k]], 0.75)),
          infl_med = median(infl),
          infl_q90 = unname(quantile(infl, 0.90))
        )
      )
    }
    
    cat(sprintf("  n=%5d done\n", n))
  }
  
  save_csv("res_cv", res)
  
  ## Multirow column for n
  n_multi <- rep("", nrow(res))
  idx_by_n <- split(seq_len(nrow(res)), res$n)
  
  for (idx in idx_by_n) {
    n_multi[idx[1]] <- sprintf(
      "\\multirow{%d}{*}{%s}",
      length(idx),
      fmt_int(res$n[idx[1]])
    )
  }
  
  latex_table(
    "tab_cv",
    headers = c(
      "$n$",
      "estimator",
      "$K^\\star$",
      "$\\mathrm{med}\\,K_{\\mathrm{CV}}\\;[\\mathrm{IQR}]$",
      "med.\\ $\\mathrm{IR}_{\\mathrm{CV}}$",
      "$q_{0.9}(\\mathrm{IR}_{\\mathrm{CV}})$"
    ),
    body = cbind(
      n_multi,
      EST_TEX[res$est],
      fmt_int(res$Kstar),
      sprintf("$%g\\;[%g,\\,%g]$", res$K_med, res$K_q1, res$K_q3),
      fmt_num(res$infl_med, 3),
      fmt_num(res$infl_q90, 3)
    ),
    caption = sprintf(
      paste0(
        "Subject-level five-fold cross-validated truncation against the oracle: ",
        "median selected $K_{\\mathrm{CV}}$ with interquartile range next to $K^\\star$. ",
        "Here $\\mathrm{IR}_{\\mathrm{CV}}=\\mathrm{ISE}(K_{\\mathrm{CV}})/\\mathrm{ISE}(K^\\star)$; ",
        "we report its median and its $0.9$-quantile ",
        "($R=%d$ replications)."
      ),
      cfg$R_cv
    ),
    label = "tab:cv",
    align = "clcccc"
  )
  
  ## Insert a separator line between the n-groups in tab_cv.tex
  tab_file <- list.files(
    ".",
    pattern = "^tab_cv\\.tex$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(tab_file) > 0L) {
    tab_file <- tab_file[which.max(file.info(tab_file)$mtime)]
    
    z <- readLines(tab_file, warn = FALSE)
    
    multirow_lines <- grep("\\\\multirow\\{", z)
    
    if (length(multirow_lines) > 1L) {
      offset <- 0L
      
      for (pos in multirow_lines[-1]) {
        p <- pos + offset
        
        if (p > 1L && !grepl("^\\s*\\\\midrule\\s*$", z[p - 1L])) {
          z <- append(z, "\\midrule", after = p - 1L)
          offset <- offset + 1L
        }
      }
      
      writeLines(z, tab_file)
    }
  } else {
    warning("Could not find tab_cv.tex to insert group separators.")
  }
  
  print(res, digits = 3)
  invisible(res)
}


# =========================================================================
#  E7 -- ORACLE g vs KERNEL PLUG-IN vs SPACING
#        (Figure fig_plugin, Table tab_plugin)
# =========================================================================
exp_plugin <- function() {
  cat("\n[E7] oracle g vs kernel plug-in vs spacing ...\n")
  set.seed(BASESEED + 7L)
  ug <- seq(0, 1, length.out = cfg$Nu_eval)
  tg <- seq(0, 1, length.out = cfg$Nt_eval)
  MU <- MU_FUN(ug, tg)
  keys <- c("oracle", "plug", "spac")
  res <- NULL
  for (n in cfg$n_plug) {
    R  <- cfg$R_plug
    IS <- setNames(lapply(keys, function(k) matrix(NA_real_, R, Kmax)), keys)
    Mv <- numeric(R)
    for (r in seq_len(R)) {
      d <- gen_data(n, ug); Mv[r] <- length(d$T)
      wu <- rep(1 / length(d$T), length(d$T))
      gh <- kde_g(d$T)
      B <- list(oracle = beta_weighted(d$Y, d$T, wu, G_DENSITY, Kmax),
                plug   = beta_weighted(d$Y, d$T, wu, gh,        Kmax),
                spac   = beta_spacing(d$Y, d$T, Kmax, h_of_M(length(d$T))))
      for (k in keys) IS[[k]][r, ] <- ise_curve(B[[k]], tg, MU)
    }
    Ks <- sapply(keys, function(k) which.min(colMeans(IS[[k]])))
    xo <- IS$oracle[, Ks["oracle"]]; xp <- IS$plug[, Ks["plug"]]
    xs <- IS$spac[,  Ks["spac"]]
    rp <- ratio_se(xp, xo); rs <- ratio_se(xs, xo)   # paired delta-method s.e.
    res <- rbind(res, data.frame(n = n, M = mean(Mv),
                                 IMSE_oracle = mean(xo), se_oracle = sd(xo) / sqrt(R),
                                 IMSE_plug   = mean(xp), se_plug   = sd(xp) / sqrt(R),
                                 IMSE_spac   = mean(xs), se_spac   = sd(xs) / sqrt(R),
                                 r_plug = unname(rp["ratio"]), se_rplug = unname(rp["se"]),
                                 r_spac = unname(rs["ratio"]), se_rspac = unname(rs["se"])))
    cat(sprintf("  n=%5d  ratios: plug-in %.3f  spacing %.3f\n",
                n, rp["ratio"], rs["ratio"]))
  }
  save_csv("res_plugin", res)
  latex_table("tab_plugin",
              headers = c("$n$", "$\\overline M$", "oracle (true $g$)",
                          "plug-in $\\widehat g$ $/$ oracle", "spacing $/$ oracle"),
              body = cbind(fmt_int(res$n), fmt_int(res$M),
                           fmt_val_se(res$IMSE_oracle, res$se_oracle),
                           fmt_val_se(res$r_plug, res$se_rplug),
                           fmt_val_se(res$r_spac, res$se_rspac)),
              caption = sprintf(paste0("Cost of (not) knowing the design density ",
                                       "$%s$: integrated risk at the oracle truncation for the infeasible ",
                                       "estimator using the true $g$, and risk ratios of the kernel plug-in ",
                                       "$\\widehat g$ and of the design-free spacing estimator to that oracle ",
                                       "($R=%d$ replications per design point; MC s.e.\\ of each ratio by the ",
                                       "delta method on the paired replications)."), G_TEX, cfg$R_plug),
              label = "tab:plugin", align = "rrccc")
  emit_fig("fig_plugin", 6.4, 4.7, function() {
    x    <- log10(res$M)
    ys   <- list(oracle = log10(res$IMSE_oracle),
                 plug   = log10(res$IMSE_plug),
                 spac   = log10(res$IMSE_spac))
    pchv <- c(oracle = 1, plug = 4, spac = 18)
    ltyv <- c(oracle = 1, plug = 4, spac = 5)
    setup_panel("$\\log_{10} M$", "$\\log_{10}\\overline{\\mathrm{IMSE}}$",
                range(x), range(unlist(ys)) + c(-0.05, 0.35))
    for (k in names(ys)) {
      lines(x, ys[[k]], lty = ltyv[k], lwd = LN_LWD)
      points(x, ys[[k]], pch = pchv[k], cex = PT_CEX)
    }
    legend_box("topright",
               c(paste0(EST_TEX["unif"], " (true $g$)"),
                 EST_TEX["plug"],
                 paste0(EST_TEX["spac"], " (no $g$)")),
               pchv, ltyv, y.intersp = 1.4)
  })
  print(res[, c("n", "M", "IMSE_oracle", "r_plug", "r_spac")], digits = 3)
  invisible(res)
}


# =========================================================================
#  E8 -- ROBUSTNESS: heteroscedastic noise / skewed design / high noise
#        (Table tab_robust)                                          
# =========================================================================
exp_robust <- function() {
  cat("\n[E8] robustness: heteroscedastic / skewed design / high noise ...\n")
  set.seed(BASESEED + 9L)
  
  n <- cfg$n_rob
  R <- cfg$R_rob
  
  ug <- seq(0, 1, length.out = cfg$Nu_eval)
  tg <- seq(0, 1, length.out = cfg$Nt_eval)
  MU <- MU_FUN(ug, tg)
  
  keys <- c("unif", "bal", "MC", "spac")
  
  scens <- list(
    base = list(
      des = DES,
      tau_fun = TAU_FUN,
      lab = "baseline"
    ),
    het = list(
      des = DES,
      tau_fun = function(t) 0.2 + 0.6 * t,
      lab = "heterosc."
    ),
    skew = list(
      des = DES_SKEW,
      tau_fun = TAU_FUN,
      lab = "skewed design"
    ),
    hinoi = list(
      des = DES,
      tau_fun = function(t) rep(1, length(t)),
      lab = "high noise"
    )
  )
  
  res <- NULL
  
  for (s in names(scens)) {
    sc <- scens[[s]]
    
    IS <- setNames(
      lapply(keys, function(k) matrix(NA_real_, R, Kmax)),
      keys
    )
    
    for (r in seq_len(R)) {
      d <- gen_data(n, ug, des = sc$des, tau_fun = sc$tau_fun)
      w <- make_weights(d, sc$des)
      
      B <- list(
        unif = beta_weighted(d$Y, d$T, w$unif, sc$des$g, Kmax),
        bal  = beta_weighted(d$Y, d$T, w$bal,  sc$des$g, Kmax),
        MC   = beta_weighted(d$Y, d$T, w$MC,   sc$des$g, Kmax),
        spac = beta_spacing(d$Y, d$T, Kmax, h_of_M(length(d$T)))
      )
      
      for (k in keys) {
        IS[[k]][r, ] <- ise_curve(B[[k]], tg, MU)
      }
    }
    
    row <- list(
      scen = s,
      lab = sc$lab
    )
    
    for (k in keys) {
      mc <- colMeans(IS[[k]])
      Ks <- which.min(mc)
      
      row[[paste0("IMSE_", k)]] <- mc[Ks]
      row[[paste0("se_",   k)]] <- sd(IS[[k]][, Ks]) / sqrt(R)
    }
    
    res <- rbind(
      res,
      as.data.frame(row, stringsAsFactors = FALSE)
    )
    
    cat(sprintf(
      "  %-6s IMSE: %s\n",
      s,
      paste(
        sprintf("%.3g", unlist(row[paste0("IMSE_", keys)])),
        collapse = " "
      )
    ))
  }
  
  rownames(res) <- res$scen
  
  save_csv("res_robust", res[, setdiff(names(res), "lab")])
  
  latex_table(
    "tab_robust",
    headers = c(
      "scenario",
      EST_TEX[keys]
    ),
    body = cbind(
      res$lab,
      fmt_val_se(res$IMSE_unif, res$se_unif),
      fmt_val_se(res$IMSE_bal,  res$se_bal),
      fmt_val_se(res$IMSE_MC,   res$se_MC),
      fmt_val_se(res$IMSE_spac, res$se_spac)
    ),
    caption = sprintf(
      paste0(
        "Integrated risk at the oracle truncation under four robustness scenarios, ",
        "$n=%d$, $R=%d$; MC s.e.\\ are reported in parentheses, in units of the common ",
        "power of ten. The scenarios are: baseline ",
        "($\\tau=1/2$, design density $g$); heteroscedastic noise ",
        "($\\tau(t)=(1+3t)/5$); skewed design ",
        "($g_2(t)=(2+6t)/5$); and high noise ",
        "($\\tau=1$)."
      ),
      n, R
    ),
    label = "tab:robust",
    align = "lcccc"
  )
  
  print(res[, c("scen", paste0("IMSE_", keys))], digits = 3)
  invisible(res)
}


# =========================================================================
#  E9 -- CONVERGENCE RATE UNDER A SECOND, BOUNDARY-SMOOTH MEAN
#        (Table tab_rate_reg, Figure fig_rate_reg)
# =========================================================================
exp_rate_reg <- function() {
  cat("\n[E9] convergence rate under the boundary-smooth mean mu_B ...\n")
  set.seed(BASESEED + 10L)
  ug <- seq(0, 1, length.out = cfg$Nu_eval)
  tg <- seq(0, 1, length.out = cfg$Nt_eval)
  MU <- MU_FUN_B(ug, tg)
  keys <- c("unif", "bal", "MC", "spac")
  res <- NULL; ISbn <- list()
  for (n in cfg$n_rate) {
    R  <- cfg$R_rate
    IS <- setNames(lapply(keys, function(k) matrix(NA_real_, R, Kmax)), keys)
    Mv <- numeric(R)
    for (r in seq_len(R)) {
      d <- gen_data(n, ug, muf = MU_FUN_B); Mv[r] <- length(d$T)   # second mean
      w <- make_weights(d)
      B <- list(unif = beta_weighted(d$Y, d$T, w$unif, G_DENSITY, Kmax),
                bal  = beta_weighted(d$Y, d$T, w$bal,  G_DENSITY, Kmax),
                MC   = beta_weighted(d$Y, d$T, w$MC,   G_DENSITY, Kmax),
                spac = beta_spacing(d$Y, d$T, Kmax, h_of_M(length(d$T))))
      for (k in keys) IS[[k]][r, ] <- ise_curve(B[[k]], tg, MU)
    }
    Mbar <- mean(Mv); row <- list(n = n, M = Mbar)
    for (k in keys) {
      mc <- colMeans(IS[[k]]); Ks <- which.min(mc)
      row[[paste0("IMSE_", k)]]  <- mc[Ks]
      row[[paste0("se_",   k)]]  <- sd(IS[[k]][, Ks]) / sqrt(R)
      row[[paste0("Kstar_", k)]] <- Ks
    }
    res <- rbind(res, as.data.frame(row)); ISbn[[length(ISbn) + 1L]] <- IS
    cat(sprintf("  n=%5d (M~%5.0f) IMSE %s  K* %s\n", n, Mbar,
                paste(sprintf("%.3g", unlist(row[paste0("IMSE_", keys)])), collapse = " "),
                paste(unlist(row[paste0("Kstar_", keys)]), collapse = "/")))
  }
  save_csv("res_rate_reg", res)
  
  ## fitted slopes with bootstrap s.e. (same convention as tab_rate)
  lM <- log10(res$M)
  sl <- sapply(keys, function(k) slope_se(lM, log10(res[[paste0("IMSE_", k)]])))
  set.seed(BASESEED + 110L)
  Bb <- 1000L; nN <- length(ISbn)
  xc <- lM - mean(lM); den <- sum(xc^2); ls_slope <- function(y) sum(xc * y) / den
  boot_sl <- matrix(NA_real_, Bb, length(keys), dimnames = list(NULL, keys))
  for (b in seq_len(Bb)) {
    yR <- matrix(NA_real_, nN, length(keys))
    for (i in seq_len(nN)) {
      Ri <- nrow(ISbn[[i]][[1]]); id <- sample.int(Ri, Ri, replace = TRUE)
      for (j in seq_along(keys)) yR[i, j] <- log10(min(colMeans(ISbn[[i]][[keys[j]]][id, , drop = FALSE])))
    }
    for (j in seq_along(keys)) boot_sl[b, j] <- ls_slope(yR[, j])
  }
  sl["se", ] <- apply(boot_sl, 2, sd)
  
  ## table (risks in 10^{-3}, slope row at the bottom)
  sc <- 1e3; fmt_r <- function(v, se) sprintf("$%.2f\\,(%.2f)$", v * sc, se * sc)
  body <- character(0)
  for (i in seq_len(nrow(res)))
    body <- c(body, paste(c(fmt_int(res$n[i]), fmt_int(res$M[i]),
                            fmt_r(unlist(res[i, paste0("IMSE_", keys)]),
                                  unlist(res[i, paste0("se_",   keys)]))), collapse = " & "))
  writeLines(c("\\begin{table}[H]", "\\centering", "\\begin{tabular}{rrcccc}", "\\toprule",
               paste0("$n$ & $\\overline M$ & ", paste(EST_TEX[keys], collapse = " & "), " \\\\"),
               "\\midrule", paste0(body, " \\\\"), "\\midrule",
               paste0("\\multicolumn{2}{r}{slope (s.e.)} & ",
                      paste(fmt_slope_se(sl["slope", ], sl["se", ]), collapse = " & "), " \\\\"),
               "\\bottomrule", "\\end{tabular}",
               sprintf(paste0("\\caption{Total integrated risk $\\overline{\\mathrm{IMSE}}$ at the ",
                              "oracle truncation for the boundary-smooth mean $\\mu_B=16\\,t^2(1-t)^2\\,\\mu$ ",
                              "versus the pooled ",
                              "sample size ($R=%d$). Risks in units of $10^{-3}$, MC s.e.\\ in parentheses. Bottom ",
                              "row: fitted slope of $\\log_{10}\\overline{\\mathrm{IMSE}}$ on $\\log_{10}M$ with ",
                              "bootstrap s.e.}"), cfg$R_rate),
               "\\label{tab:rate_reg}", "\\end{table}"),
             file.path(OUTDIR, "tab_rate_reg.tex"))
  cat("  wrote tab_rate_reg.tex\n")
  
  emit_fig("fig_rate_reg", 6.4, 5.0, function() {
    x <- log10(res$M); ys <- lapply(keys, function(k) log10(res[[paste0("IMSE_", k)]]))
    setup_panel("$\\log_{10} M$", "$\\log_{10}\\overline{\\mathrm{IMSE}}$",
                range(x), range(unlist(ys)) + c(-0.05, 0.45))
    for (j in seq_along(keys)) { k <- keys[j]
    lines(x, ys[[j]], lty = EST_LTY[k], lwd = LN_LWD); points(x, ys[[j]], pch = EST_PCH[k], cex = PT_CEX) }
    legend_box("topright", sprintf("%s  ($%.2f$)", EST_TEX[keys], sl["slope", keys]), EST_PCH[keys], EST_LTY[keys])
  })
  invisible(list(res = res, slopes = sl))
}


# =========================================================================
#  VISUALISATION  --  one realisation of the data-generating mechanism
#  (Figures fig_mean, fig_subjects, fig_obs)
# =========================================================================
make_viz <- function() {
  cat("\n[viz] data-mechanism figures ...\n")
  set.seed(BASESEED + 1L)
  MUv <- MU_FUN(u_viz, t_viz)
  X1 <- make_Xi()(u_viz, t_viz)               # two independent subjects
  X2 <- make_Xi()(u_viz, t_viz)
  zl <- range(MUv + X1, MUv + X2)
  ## one further subject, with its observed curves (X + scaled noise):
  Xi <- make_Xi()
  mi <- 0L
  while (mi < 2L || mi > m_cap) mi <- rpois(1L, lambda_m)   ## truncated
  Ti <- sample_T(mi)
  Yi <- MU_FUN(u_fine, Ti) + Xi(u_fine, Ti) +
    gen_noise(u_fine, mi) * matrix(TAU_FUN(Ti), length(u_fine), mi, byrow = TRUE)
  Zi <- Xi(u_viz, t_viz)
  wire <- function(Z, zlim, main = "")        # a grey wireframe surface
    persp(u_viz, t_viz, Z, theta = 35, phi = 25, expand = 0.72, zlim = zlim,
          xlab = "u", ylab = "t", zlab = "", ticktype = "detailed",
          cex.axis = 0.55, cex.lab = 0.8, border = NA, col = "grey85",
          shade = 0.55, main = main, cex.main = 0.95)
  emit_fig("fig_mean", 5.6, 4.4, function() {
    par(mar = c(1.0, 1.0, 1.0, 0.6)); wire(MUv, range(MUv))
  })
  emit_fig("fig_subjects", 6.8, 3.6, function() {
    par(mfrow = c(1, 2), mar = c(1.0, 1.0, 1.6, 0.6))
    wire(MUv + X1, zl, "$\\mu+X_1$")
    wire(MUv + X2, zl, "$\\mu+X_2$")
  })
  emit_fig("fig_obs", 6.8, 3.8, function() {
    par(mfrow = c(1, 2))
    par(mar = c(1.0, 1.0, 1.6, 0.6))
    wire(Zi, range(Zi), "$X_i(u,t)$")
    par(mar = c(4.0, 4.2, 1.6, 0.8), mgp = c(2.4, 0.8, 0))
    cols <- grey.colors(mi, start = 0.10, end = 0.70)
    matplot(u_fine, Yi, type = "l", lty = 1, lwd = 0.9, col = cols,
            xlab = "$u$", ylab = "$Y_{ij}(u)$", las = 1, bty = "l",
            main = "observed curves", cex.main = 0.95)
    legend("topright", sprintf("$T_{ij}=%.2f$", Ti), col = cols,
           lty = 1, lwd = 1.2, bty = "n", cex = 0.60)
  })
}


# =========================================================================
#  DRIVER  --  run everything in order, then print a one-screen summary
# =========================================================================
t0 <- proc.time()[3]
make_viz()
E1 <- exp_rate()     # rate + K*/h* + LL benchmark   (tab_rate, tab_Kstar, fig_rate_total)
E2 <- exp_unbias()   # unbiasedness                  (tab_unbias)
E3 <- exp_riskK()    # risk vs K, plug-in            (tab_riskK, fig_riskK)
E4 <- exp_case12()   # Case 1 vs Case 2              (tab_case12, fig_case12)
E5 <- exp_design()   # design-term decay             (tab_design, fig_design)
E6 <- exp_cv()       # cross-validated truncation    (tab_cv)            
E7 <- exp_plugin()   # oracle/plug-in/spacing        (tab_plugin, fig_plugin)
E8 <- exp_robust()   # robustness scenarios          (tab_robust)
E9 <- exp_rate_reg() # rate under a second mean regularity   (tab_rate_reg, fig_rate_reg)


cat("\n==================== SUMMARY ====================\n")
cat("[E1] total-risk slopes (series -> -0.75; LL: leading -0.80):\n")
print(round(E1$slopes, 3))
cat("[E1] K*-slopes (prediction +0.25) and h_t*-slope (prediction -0.20):\n")
print(round(E1$kslopes, 3)); print(round(E1$hslope, 3))
cat("[E1] selected h_u* per n (0 = no u-smoothing helps):",
    paste(format(E1$Kres$hu_ll, digits = 3), collapse = "  "), "\n")
cat("[E5] design-term slopes (prediction -1 / -2 / -3):\n")
print(round(E5$slopes, 3))
cat(sprintf("[E4] risk ratio mu^(1)/mu^(2): balanced %.6f (exactly 1), unbalanced %.3f\n",
            E4$res$ratio[1], E4$res$ratio[2]))
cat("[E6] cross-validated truncation vs oracle (median inflation):\n")
print(E6[, c("n", "est", "Kstar", "K_med", "infl_med", "infl_q90")], digits = 3)
cat("[E7] plug-in/oracle and spacing/oracle risk ratios:\n")
print(round(E7[, c("n", "r_plug", "r_spac")], 3))
cat("[E8] best scheme per scenario (by IMSE at K*):\n")
imse_cols <- paste0("IMSE_", c("unif", "bal", "MC", "spac"))
print(data.frame(scenario = rownames(E8),
                 best = c("unif", "bal", "MC", "spac")[apply(E8[, imse_cols], 1, which.min)]))
cat("[E9] rate slopes under boundary-smooth mean (prediction -7/8 = -0.875):\n")
print(round(E9$slopes["slope", ], 3))
cat(sprintf("\nTotal wall time: %.1f s\nOutputs written to: %s\n",
            proc.time()[3] - t0, normalizePath(OUTDIR)))