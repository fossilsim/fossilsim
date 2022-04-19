## Test environments
* local Ubuntu 20.04 install, R 4.1.3
* R-hub
* win-builder, devel

## R CMD check results

0 errors | 0 warnings | 1 note

There was 1 note on R-hub & win-builder:
 * Found the following (possibly) invalid URLs:
    URL: https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13170
      From: README.md
      Status: 503
      Message: Service Unavailable
URL works fine in local check and in several browsers.

## Reverse dependencies

There are no reverse dependencies, as verified using ``devtools::revdep()``.

