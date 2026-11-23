# Independent reference for the Poisson-lognormal (PLN) and bivariate PLN pmf.
# Never calls poilog. Sourced by make-reference.R and kernel-reference.R in this
# directory.
#
# Univariate: q(n; mu, sig) = int phi(u) dpois(n, exp(mu + sig*u)) du.
# In u the log integrand
#   h(u) = n*(mu+sig*u) - exp(mu+sig*u) - lgamma(n+1) - u^2/2 - log(2*pi)/2
# is strictly concave with h'' = -sig^2*exp(mu+sig*u) - 1 <= -1, so it decays at
# least like a standard normal around its mode. We find the mode u* (root of the
# monotone h'), find L < u* < R where h drops by DROP = 60 nats below h(u*)
# (truncation error < (R-L)*exp(-60) relative), and integrate exp(h - h(u*)) on [L, R]
# by the composite trapezoid rule with 4001 nodes (geometrically convergent for a
# smooth integrand that is ~exp(-60) at both ends). With check = TRUE it must agree
# to REF_AGREE relative with the same rule at 2001 nodes and with stats::integrate
# (adaptive Gauss-Kronrod, rel.tol 1e-11; compared whenever integrate reports OK).
# log q = h(u*) + log(integral).
# Everything is on the log scale, so there is no underflow for tiny q.

DROP <- 60
REF_AGREE <- 1e-9

ref_h <- function(u, n, mu, sig) {
  eta <- mu + sig * u
  n * eta - exp(eta) - lgamma(n + 1) - u^2 / 2 - 0.5 * log(2 * pi)
}

ref_mode <- function(n, mu, sig) {
  if (sig == 0) return(0)
  dh <- function(u) n * sig - sig * exp(mu + sig * u) - u
  # dh is strictly decreasing; dh(u) = 0 lies in [-B, B] with B = n*sig + |...|
  lo <- -1; while (dh(lo) < 0) lo <- lo * 2
  hi <- 1;  while (dh(hi) > 0) hi <- hi * 2
  uniroot(dh, c(lo, hi), tol = 1e-14, maxiter = 1000)$root
}

ref_lq1 <- function(n, mu, sig, check = TRUE) {
  stopifnot(length(n) == 1, n >= 0, sig >= 0)
  if (sig == 0) return(dpois(n, exp(mu), log = TRUE))
  us <- ref_mode(n, mu, sig)
  hs <- ref_h(us, n, mu, sig)
  g <- function(u) ref_h(u, n, mu, sig) - hs + DROP
  # h - hs is <= -(u-us)^2/2, so the drop points are within us -/+ sqrt(2*DROP)
  w <- sqrt(2 * DROP) + 1
  L <- uniroot(g, c(us - w, us), tol = 1e-12)$root
  R <- uniroot(g, c(us, us + w), tol = 1e-12)$root
  f <- function(u) exp(ref_h(u, n, mu, sig) - hs)
  trap <- function(m) {
    u <- seq(L, R, length.out = m); fu <- f(u)
    (u[2] - u[1]) * (sum(fu) - (fu[1] + fu[m]) / 2)
  }
  B <- trap(4001)
  if (check) {
    B2 <- trap(2001)
    if (abs(B - B2) / B > REF_AGREE)
      stop(sprintf("trapezoid not converged: n=%g mu=%g sig=%g", n, mu, sig))
    A <- integrate(f, L, R, rel.tol = 1e-11, abs.tol = 0, subdivisions = 2000L,
                   stop.on.error = FALSE)
    if (A$message == "OK" && abs(A$value - B) / B > REF_AGREE)
      stop(sprintf("reference rules disagree: n=%g mu=%g sig=%g A=%.15g B=%.15g",
                   n, mu, sig, A$value, B))
  }
  hs + log(B)
}

ref_lq <- function(n, mu, sig, check = TRUE)
  vapply(n, ref_lq1, 0, mu = mu, sig = sig, check = check)

# Bivariate: q(n1, n2) = int phi(u) dpois(n1, exp(mu1+sig1*u)) *
#                          q(n2; mu2 + rho*sig2*u, sig2*sqrt(1-rho^2)) du
# (condition the second log-intensity on the first standard normal u).
# The integrand is log-concave in u (PLN pmf is log-concave in its location by
# Prekopa), so the same mode and drop-point scheme applies; the integral over
# [L, R] is computed by stats::integrate (rel.tol 1e-10). The inner q is the
# reference univariate pmf above, without its self-check.
ref_lq2_1 <- function(n1, n2, mu1, mu2, sig1, sig2, rho) {
  s2c <- sig2 * sqrt(max(0, 1 - rho^2))
  lk <- function(u) {
    vapply(u, function(ui) {
      eta <- mu1 + sig1 * ui
      n1 * eta - exp(eta) - lgamma(n1 + 1) - ui^2 / 2 - 0.5 * log(2 * pi) +
        ref_lq1(n2, mu2 + rho * sig2 * ui, s2c, check = FALSE)
    }, 0)
  }
  # mode: lk is concave; bracket it with the analogous bounds
  opt <- optimize(lk, c(-40, 40), maximum = TRUE, tol = 1e-10)
  us <- opt$maximum; hs <- opt$objective
  if (abs(us) > 39) stop("bivariate reference: mode outside [-40, 40]")
  g <- function(u) lk(u) - hs + DROP
  w <- sqrt(2 * DROP) + 1
  L <- uniroot(g, c(us - w, us), tol = 1e-10)$root
  R <- uniroot(g, c(us, us + w), tol = 1e-10)$root
  f <- function(u) exp(lk(u) - hs)
  A <- integrate(f, L, R, rel.tol = 1e-10, abs.tol = 0, subdivisions = 2000L)$value
  hs + log(A)
}
ref_lq2 <- function(n1, n2, mu1, mu2, sig1, sig2, rho)
  mapply(ref_lq2_1, n1, n2, MoreArgs = list(mu1 = mu1, mu2 = mu2, sig1 = sig1,
                                            sig2 = sig2, rho = rho))

# Reference zero-truncated log-likelihood for counts n (univariate)
ref_loglik <- function(n, mu, sig, zTrunc = TRUE, check = FALSE) {
  tb <- table(n); un <- as.numeric(names(tb)); nr <- as.numeric(tb)
  lq <- ref_lq(un, mu, sig, check = check)
  b <- if (zTrunc) ref_lpos1(mu, sig) else 0
  sum((lq - b) * nr)
}

# ---- log P(N > 0), computed directly (not as log(1 - q0), which cancels when q0 ~ 1)
# log1mexp(lambda) = log(1 - exp(-lambda)), lambda = exp(eta), accurate for all eta
ref_log1mexp_exp <- function(eta) {
  lam <- exp(eta)
  ifelse(eta < -30, eta - lam / 2,
         ifelse(lam > log(2), log1p(-exp(-lam)), log(-expm1(-lam))))
}
# generic: log int exp(lf(u)) du for a concave lf with lf'' <= -1, mode in [lo, hi]
ref_logint_concave <- function(lf, lo, hi) {
  o <- optimize(lf, c(lo, hi), maximum = TRUE, tol = 1e-12)
  us <- o$maximum; hs <- o$objective
  if (us - lo < 1e-3 || hi - us < 1e-3) stop("mode at the edge of the search bracket")
  g <- function(u) lf(u) - hs + DROP
  w <- sqrt(2 * DROP) + 1
  L <- uniroot(g, c(us - w, us), tol = 1e-12)$root
  R <- uniroot(g, c(us, us + w), tol = 1e-12)$root
  u <- seq(L, R, length.out = 4001); fu <- exp(lf(u) - hs)
  hs + log((u[2] - u[1]) * (sum(fu) - (fu[1] + fu[4001]) / 2))
}
ref_lpos1 <- function(mu, sig) {
  if (sig == 0) return(ref_log1mexp_exp(mu))
  lf <- function(u) ref_log1mexp_exp(mu + sig * u) - u^2 / 2 - 0.5 * log(2 * pi)
  ref_logint_concave(lf, -50, max(50, 2 * sig))
}
# log P(N1 + N2 > 0) = log[ P(N1 > 0) + int phi(u) exp(-l1(u)) P(N2 > 0 | u) du ]
ref_lpos2 <- function(mu1, mu2, sig1, sig2, rho) {
  sc <- sig2 * sqrt(max(0, 1 - rho^2))
  lf <- function(u) vapply(u, function(ui) -exp(mu1 + sig1 * ui) +
    ref_lpos1(mu2 + rho * sig2 * ui, sc) - ui^2 / 2 - 0.5 * log(2 * pi), 0)
  a <- ref_lpos1(mu1, sig1); b <- ref_logint_concave(lf, -50, 50)
  max(a, b) + log1p(exp(-abs(a - b)))
}

# Reference zero-truncated log-likelihood for pairs of counts (n1, n2) (bivariate).
# The truncation term is log1p(-q0) with q0 = q(0, 0) from the reference pmf; this
# requires q0 <= 1/2, where log1p(-q0) has no cancellation (ref_lpos2 is the direct
# form, for any q0, but much slower).
ref_loglik2 <- function(n1, n2, mu1, mu2, sig1, sig2, rho) {
  key <- paste(n1, n2); un <- !duplicated(key)
  nr <- as.numeric(table(key)[key[un]])
  lq <- ref_lq2(n1[un], n2[un], mu1, mu2, sig1, sig2, rho)
  q0 <- exp(ref_lq2(0, 0, mu1, mu2, sig1, sig2, rho))
  if (q0 > 0.5) stop("ref_loglik2: q0 > 1/2, use ref_lpos2 for the truncation term")
  sum(lq * nr) - length(n1) * log1p(-q0)
}
