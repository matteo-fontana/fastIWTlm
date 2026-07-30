# Standalone C++ command-line IWT (non-Rcpp)

This directory holds the original **standalone** C++ implementation of the
two-sample Interval-Wise Test. It does **not** use Rcpp: it reads/writes
plain text files and is driven from the command line via
[GetPot](https://getpot.sourceforge.net/). It is kept here for reference and
is **not** compiled as part of the R package (only `src/` is).

## Build

Requires [Eigen](https://eigen.tuxfamily.org/) (headers via `$mkEigenInc`)
and a C++11 compiler with OpenMP.

```sh
make            # builds ./main
```

## Run

```sh
./main B=1000 THREADS=4 alt=two.sided n1=50 n2=50 p=100
```

For everyday use from R, prefer the compiled `fastIWT2()` engine or the
pure-R `IWT2()` in the package instead.
