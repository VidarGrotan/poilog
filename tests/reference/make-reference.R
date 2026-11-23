# Recompute the reference values hard-coded in tests/testthat/.
#
# The values come from ref.R in this directory, an implementation of the
# Poisson-lognormal pmf that is independent of poilog (it never calls the
# package): trapezoid rule on the log scale around the mode of the integrand,
# checked against stats::integrate; see the header of ref.R.
#
# Run from the package root:   Rscript tests/reference/make-reference.R
# It prints every value with the test file that uses it. It takes about 20 minutes,
# most of them for the bivariate maximum at the end.
# This directory is not part of the built package (see .Rbuildignore).

source("tests/reference/ref.R")
p <- function(...) cat(sprintf(...), "\n")

p("== test-dpoilog.R: log q(n; mu, sig) ==")
u <- data.frame(n   = c(1, 3, 10, 100, 2, 5, 0, 0, 0, 1000, 1e6),
                mu  = c(0, 1, 2, 4, 3, 10, -2, 0, -8, 0, 10),
                sig = c(1, 1.5, 1, 1, 0.5, 2, 5, 8, 3, 0.5, 1))
for (i in seq_len(nrow(u))) p("n = %g, mu = %g, sig = %g: %.15g", u$n[i], u$mu[i], u$sig[i],
                               ref_lq(u$n[i], u$mu[i], u$sig[i]))

p("== test-dpoilog.R: log P(N > 0) ==")
for (x in list(c(-200, 25), c(0, 1), c(-5, 2))) p("mu = %g, sig = %g: %.15g", x[1], x[2], ref_lpos1(x[1], x[2]))

p("== test-dbipoilog.R: log q(n1, n2; mu1, mu2, sig1, sig2, rho) ==")
b <- rbind(c(3, 5, 1, 1, 1, 1, 0.5), c(0, 2, 0.5, 1, 1.5, 1, -0.5), c(0, 0, -2, -2, 5, 5, 0.5),
           c(0, 500, 3, 3, 2, 2, 0.95), c(500, 0, 3, 3, 2, 2, 0.95), c(4, 4, 1, 1, 1, 1, 1))
for (i in seq_len(nrow(b))) p("(%s): %.15g", paste(b[i, ], collapse = ", "), do.call(ref_lq2, as.list(b[i, ])))

p("== test-dbipoilog.R: log P(N1 + N2 > 0) ==")
# q(0, 0) is below 1/2 for the first two, 0.742 for the third, and 0.5047 and 0.4954
# for the last two, on both sides of 1/2
for (x in list(c(-1, -1, 2, 2, 0.5), c(0, 1, 1, 1.5, -0.5), c(-3, -3, 2, 2, 0.5),
               c(-1.62, -1.5692, 2, 2, 0.5), c(-1.52, -1.5692, 2, 2, 0.5)))
  p("(%s): %.15g", paste(x, collapse = ", "), do.call(ref_lpos2, as.list(x)))

p("== test-mle.R: maximum of the zero-truncated log-likelihood of `sparse` ==")
# Profile over sig (mu maximised for each sig, in a bracket that widens like sig^2,
# because along the ridge of this likelihood mu falls like -sig^2), then polish
# from the best profile point with Nelder-Mead.
sparse <- c(rep(1, 18), rep(2, 6), 4, 5, 5, 6, 6, 7, 16, 45)
prof <- t(vapply(exp(seq(log(0.05), log(30), length.out = 60)), function(s) {
  o <- optimize(function(m) ref_loglik(sparse, m, s), c(-5 - 2 * s^2, 5), maximum = TRUE)
  c(s, o$maximum, o$objective)
}, numeric(3)))
best <- prof[which.max(prof[, 3]), ]
nm <- optim(c(best[2], log(best[1])), function(z) -ref_loglik(sparse, z[1], exp(z[2])),
            control = list(reltol = 1e-12, maxit = 2000))
p("mu = %.4f, sig = %.4f, log-likelihood = %.10f", nm$par[1], exp(nm$par[2]), -nm$value)

p("== test-mle.R: the zero-truncated log-likelihood of `sparse` far out on the ridge ==")
p("mu = -250, sig = 25: %.15g", ref_loglik(sparse, -250, 25, check = TRUE))

p("== test-mle.R: maximum of the zero-truncated log-likelihood of `few` ==")
# the same profile search as for `sparse`
few <- c(1, 1, 2, 3, 5, 8, 13)
prof <- t(vapply(exp(seq(log(0.05), log(30), length.out = 60)), function(s) {
  o <- optimize(function(m) ref_loglik(few, m, s), c(-5 - 2 * s^2, 5), maximum = TRUE)
  c(s, o$maximum, o$objective)
}, numeric(3)))
best <- prof[which.max(prof[, 3]), ]
nm <- optim(c(best[2], log(best[1])), function(z) -ref_loglik(few, z[1], exp(z[2])),
            control = list(reltol = 1e-12, maxit = 2000))
p("mu = %.4f, sig = %.4f, log-likelihood = %.10f", nm$par[1], exp(nm$par[2]), -nm$value)

p("== test-dbipoilog.R: log P(N1 + N2 > 0) where q(0, 0) is small and rho is near -1 ==")
x <- c(-0.6286, 8.1117, 0.6698, 2.1814, -0.99988)
lq0 <- do.call(ref_lq2, as.list(c(0, 0, x)))
p("(%s): log q(0, 0) = %.17g, log1p(-q0) = %.17g", paste(x, collapse = ", "), lq0, log1p(-exp(lq0)))

p("== test-mle.R: maximum of the zero-truncated log-likelihood of `xy4` (sig about 4) ==")
# 34 species simulated with mu1 = mu2 = -1, sig1 = sig2 = 4, rho = 0.5 (60 species, of
# which 26 were seen in neither sample; drawn by set.seed(13); rbipoilog(S = 60, mu1 = -1,
# mu2 = -1, sig1 = 4, sig2 = 4, rho = 0.5)). The maximum is found by BFGS with numerical
# gradients on the reference likelihood, in (mu1, mu2, log sig1, log sig2, atanh rho),
# started from the parameters of the simulation. This takes most of the run time
# (about a quarter of an hour).
xy4 <- cbind(c(6, 0, 455, 1, 33, 51, 0, 0, 34, 1, 102, 1, 0, 4, 2, 580, 4, 0, 0, 0, 1, 6,
               237, 1, 231, 0, 86, 8, 1, 0, 0, 9, 123, 6),
             c(0, 3, 5, 0, 123, 0, 21, 11, 1, 1, 2578, 0, 4, 16, 2, 38, 138, 3, 11, 4, 2, 1,
               71, 0, 14, 8, 2059, 438, 0, 1, 2, 23, 31, 2))
nll4 <- function(z) -ref_loglik2(xy4[, 1], xy4[, 2], z[1], z[2], exp(z[3]), exp(z[4]), tanh(z[5]))
o4 <- optim(c(-1, -1, log(4), log(4), atanh(0.5)), nll4, method = "BFGS",
            control = list(reltol = 1e-12, maxit = 200, ndeps = rep(1e-4, 5)))
th4 <- c(o4$par[1:2], exp(o4$par[3:4]), tanh(o4$par[5]))
g4 <- vapply(1:5, function(j) { e <- replace(numeric(5), j, 1e-4)
  (nll4(o4$par + e) - nll4(o4$par - e)) / 2e-4 }, 0)
p("convergence %d after %d evaluations; largest |gradient| in z %.1e",
  o4$convergence, o4$counts[[1]], max(abs(g4)))
p("mu1 = %.6f, mu2 = %.6f, sig1 = %.6f, sig2 = %.6f, rho = %.6f, log-likelihood = %.10f",
  th4[1], th4[2], th4[3], th4[4], th4[5], -o4$value)
