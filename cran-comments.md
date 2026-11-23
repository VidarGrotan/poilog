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
* win-builder (R-devel) <!-- TODO: run on the final build -->

## R CMD check results

0 errors | 0 warnings | 1 note

* "checking for future file timestamps ... NOTE: unable to verify current
  time" appears only locally, where the time server cannot be reached.
<!-- TODO: confirm on the final build and on win-builder -->

## Reverse dependencies

One reverse dependency, sads.
<!-- TODO: re-run R CMD check of sads 0.6.5 against the final build, and
inform the sads maintainers, before submitting -->
R CMD check of sads 0.6.5 gives the same result with this version as with
0.4.2.1.
