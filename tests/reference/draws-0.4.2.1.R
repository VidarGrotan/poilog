# Seeded random draws from poilog 0.4.2.1, for tests/testthat/test-random.R, which
# checks that the current version gives identical draws for the same seed.
#
# Needs poilog 0.4.2.1 installed in some library (one of its own, if the default
# library holds another version). Run from the package root:
#   Rscript tests/reference/draws-0.4.2.1.R <library with poilog 0.4.2.1>
# and paste the printed vectors into test-random.R.
lib <- commandArgs(TRUE)[1]
library(poilog, lib.loc = lib)
stopifnot(packageVersion("poilog", lib.loc = lib) == "0.4.2.1")
show <- function(name, x) cat(sprintf("%s <- c(%s)\n", name, paste(x, collapse = ", ")))

set.seed(101); show("d_keep0", rpoilog(12, 1, 1.5, keep0 = TRUE))
set.seed(102); show("d_pos", rpoilog(30, 0, 2))
set.seed(103); show("d_condS", rpoilog(10, -1, 2, nu = 0.5, condS = TRUE))
set.seed(104); xy <- rbipoilog(8, 0, 0.5, 1, 1.2, 0.6, keep0 = TRUE)
show("b_keep0_1", xy[, 1]); show("b_keep0_2", xy[, 2])
set.seed(105); xy <- rbipoilog(10, -1, -0.5, 2, 1.5, -0.3, nu1 = 0.7, condS = TRUE)
show("b_condS_1", xy[, 1]); show("b_condS_2", xy[, 2])
