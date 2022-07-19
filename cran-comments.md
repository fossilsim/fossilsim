## Test environments
* local Windows install, R 4.2.0
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
  URL: https://doi.org/10.1073/pnas.1319091111
    From: inst/doc/fossils.html
    Status: 503
    Message: Service Unavailable
URLs work fine in local check and in several browsers.

## Reverse dependencies

The reverse dependency "FossilSimShiny" has been checked with the new changes.

