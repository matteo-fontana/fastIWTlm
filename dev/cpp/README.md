# Experimental C++ (not built)

These translation units are earlier/experimental Rcpp variants of the
two-sample engine, kept for reference. They are **excluded from the package
build** (via `.Rbuildignore`) because they define the same free functions
(`print_mex`, `check_data`, `compute_T2`, ...) as `src/fastIWT2.cpp` and
would cause duplicate-symbol link errors if compiled together.

- `IWT2C_exp.cpp` — experimental Rcpp two-sample engine with verbose timing.
- `IWT2OMP.cpp`   — Rcpp engine that reads its inputs from text files.

The production engine is `src/fastIWT2.cpp` (exported as `fastIWT2()`).
