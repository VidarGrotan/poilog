# Reference values for tests/testthat/test-kernel.R.
#
# They come from an independent quadrature (tests/reference/ref.R) and never from
# poilog: the univariate pmf and log P(N > 0) by the trapezoid rule on the log scale
# between the points where the integrand has dropped 60 nats, with 4001 nodes; the
# pmf values here are self-checked against 2001 nodes and stats::integrate
# (check = TRUE), log P(N > 0) is not. The bivariate pmf uses stats::integrate over the
# same kind of range, and log P(N1 + N2 > 0) the trapezoid rule (see ref.R).
# To recompute, from the package root:
#   Rscript tests/reference/kernel-reference.R
# and paste the printed vectors into test-kernel.R.
source("tests/reference/ref.R")
show <- function(name, x) cat(sprintf("%s <- c(%s)\n", name, paste(sprintf("%.17g", x), collapse = ", ")))

# log q(n; mu, sig): points where two successive trapezoid levels agree before the step
# resolves the integrand (the kernel needs its resolution test there), and points
# among the least accurate of dense comparisons of the kernel with this reference
u <- data.frame(
  n   = c(0, 1, 0, 0, 0, 1e5),
  mu  = c(-5.631, -9.63, -6.89, -51.96, -48.1, -56.15),
  sig = c(2, 2, 2, 10, 15, 30))
u$lq <- mapply(function(n, mu, sig) ref_lq(n, mu, sig, check = TRUE), u$n, u$mu, u$sig)
show("q_n", u$n); show("q_mu", u$mu); show("q_sig", u$sig); show("q_ref", u$lq)

# log P(N > 0): points near the kink of the integrand where exp(eta) passes 1, at small
# and large sig, among the least accurate of dense comparisons and of a search over mu
# for the largest error
p <- data.frame(mu = c(-18.73, 2.73, -18.09, 3.7072671199841, -2.82612498245936, -2.83),
                sig = c(3, 1, 15, 50, 5, 5))
p$lp <- mapply(ref_lpos1, p$mu, p$sig)
show("pos_mu", p$mu); show("pos_sig", p$sig); show("pos_ref", p$lp)

# bivariate at large sig (the regime where the change of variables matters)
b <- data.frame(n1 = c(0, 0, 3), n2 = c(0, 1, 2), mu1 = 1, mu2 = 0.5, sig = 50, rho = 0.5)
b$lq <- mapply(function(n1, n2, mu1, mu2, sig, rho) ref_lq2(n1, n2, mu1, mu2, sig, sig, rho),
               b$n1, b$n2, b$mu1, b$mu2, b$sig, b$rho)
show("biv_n1", b$n1); show("biv_n2", b$n2); show("biv_ref", b$lq)
show("bipos_ref", ref_lpos2(1, 0.5, 50, 50, 0.5))

# |rho| = 1 with a very negative mu2: q(0, 1) = int phi(u) dpois(0, e^u) dpois(1,
# e^(-800 + u)) du = e^-800 int phi(u) e^u exp(-e^u) du (the factor exp(-e^(-800 + u))
# is 1 to double precision), i.e. -800 + log q(1; 0, 1)
show("rho1_ref", -800 + ref_lq(1, 0, 1, check = TRUE))

# rho = 1 with a large second count, where a step of the mode search overflows exp():
# q(0, 500) at mu1 = 0, mu2 = -20, sig = 2
show("rho1_500_ref", ref_lq2(0, 500, 0, -20, 2, 2, 1))

# the gradient of log q(n; mu, sig), d/dmu = E[n - lambda] and d/dsig = E[(n - lambda) u]
# under the normalised integrand phi(u) dpois(n, lambda), lambda = exp(mu + sig u), at a
# point where the stopping rule needs its level-agreement test (b) for the gradient:
# without (b) the kernel's gradient is off by about 1e-8 there. The trapezoid rule on
# 20001 nodes over 14 units either side of the mode, checked here against 40001 nodes
# (and once against stats::integrate: agreement to 1e-15)
ref_score1 <- function(n, mu, sig, nodes = 20001) {
  m <- uniroot(function(u) n * sig - sig * exp(mu + sig * u) - u, c(-50, 50), tol = 1e-14)$root
  u <- seq(m - 14, m + 14, length.out = nodes)
  lf <- n * (mu + sig * u) - exp(mu + sig * u) - u^2 / 2
  w <- exp(lf - max(lf)); a <- n - exp(mu + sig * u)
  c(sum(w * a), sum(w * a * u)) / sum(w)
}
g <- ref_score1(3, -6.4, 2)
stopifnot(max(abs(g - ref_score1(3, -6.4, 2, 40001))) < 1e-14)
show("score_ref", g)
