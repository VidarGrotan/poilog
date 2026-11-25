## Submission

This is a bug-fix release with one new argument (0.4.2.1 -> 0.5.0). It
corrects inaccurate probabilities in `dpoilog()` and `dbipoilog()`, and
failures in the maximum-likelihood and bootstrap code, adds a `log` argument
to `dpoilog()` and `dbipoilog()`, and adds a test suite; see NEWS.md, which
also lists the calls that now stop with an error. The C code no longer
allocates memory, which supersedes the `R_Calloc`/`R_Free` change the CRAN
team made in 0.4.2.1.

## Test environments

* local Windows 11, R 4.4.2
* GitHub Actions: macOS (arm64), Windows and Ubuntu with R 4.6.2; Ubuntu with
  R-devel and with R 4.5.3
* win-builder (R-devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* "checking for future file timestamps ... NOTE: unable to verify current
  time" appears only locally, where the time server cannot be reached. The
  checks on GitHub Actions and win-builder give no notes.

## Reverse dependencies

One reverse dependency, sads. R CMD check of sads 0.6.5 gives the same result
with this version as with 0.4.2.1 (Status: OK with both). Its maintainer has
been informed.
