# QDotsCpp
A remake of QDots Python for more performance.

Uses armadillo-4.650.4 as a Linear Algebra library. (The latest version gave segmentation faults).

Using OpenBLAS as a faster alternative to BLAS.

Plotting will be done using a pipe to Gnuplot using gnuplot-iostream. This pipe supports the armadillo data types.
