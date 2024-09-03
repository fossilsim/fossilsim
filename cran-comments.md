## Test environments
* local Ubuntu 20.04 install, R 4.4.1
* win-builder, devel and release

## R CMD check results

0 errors | 0 warnings | 1 note

The checks identified one note concerning the use of ::: to access unexported functions of the ggtree package.
Our code extends ggtree plotting to a new data structure and so is tightly linked to the internals of ggtree.
This dependency cannot be removed other than by copying the functions to our own package, which would duplicate a large part of ggtree.

## Reverse dependencies

The reverse dependency "FossilSimShiny" has been checked with the new changes.

