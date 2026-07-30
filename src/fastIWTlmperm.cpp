// -----------------------------------------------------------------------------
// fastIWTlm_perm: RcppEigen + OpenMP acceleration of the permutation loop of
// the functional-regression IWT implemented in R/fastIWTlm.R.
//
// The design matrix at each abscissa point is FIXED across permutations (only
// the response is permuted / reconstructed), so the projection
// P_j = (X_j' X_j)^{-1} X_j' and the diagonal of (X_j' X_j)^{-1} are computed
// ONCE per point and reused for every one of the B permutations. This is the
// speed-up the naive R loop (which re-runs lm.fit each iteration) leaves on the
// table.
//
// Randomness is supplied from R as a B x n integer matrix of 1-based
// permutation indices (drawn with sample()), so results are bit-for-bit
// reproducible against the pure-R engine under a common set.seed().
//
// Semantics reproduced exactly (see R/fastIWTlm.R):
//   calculateTglob(regr) = sum((fitted - mean(fitted))^2) / (nvar * resvar)
//   calculateTpart(regr)[i] = |coef_i| / sqrt(diag((X'X)^{-1})_i * resvar)
// with resvar = sum(residuals^2) / (n - k),  k = nvar + 1.
// -----------------------------------------------------------------------------

#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <vector>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

//[[Rcpp::export]]
List fastIWTlm_perm(List Xlist,
                    NumericMatrix Y,
                    IntegerMatrix PermMat,
                    int method_residuals,
                    List fittedPart,
                    List residPart,
                    int nvar,
                    int threads) {

  const int J = Y.ncol();      // number of abscissa points
  const int n = Y.nrow();      // number of statistical units
  const int B = PermMat.nrow(); // number of permutations
  const int k = nvar + 1;      // number of regression coefficients
  const double df = static_cast<double>(n - k);

  // --- raw pointers extracted once, outside the parallel region -------------
  // (no R API is touched inside the OpenMP loop; only plain-memory reads/writes)
  const double* Yptr = &Y[0];            // column-major, element(i,j) = Yptr[i + j*n]
  const int*    Pptr = &PermMat[0];      // column-major, element(p,i) = Pptr[p + i*B], 1-based

  std::vector<const double*> fitPtr(method_residuals ? k : 0);
  std::vector<const double*> resPtr(method_residuals ? k : 0);
  if (method_residuals) {
    for (int ii = 0; ii < k; ii++) {
      NumericMatrix f = fittedPart[ii];   // n x J
      NumericMatrix r = residPart[ii];    // n x J
      fitPtr[ii] = &f[0];
      resPtr[ii] = &r[0];
    }
  }

  // --- per-point precomputation (serial) -----------------------------------
  std::vector<Eigen::MatrixXd> Xmats(J);   // n x k
  std::vector<Eigen::MatrixXd> Pmats(J);   // k x n  = (X'X)^{-1} X'
  std::vector<Eigen::VectorXd> diagInv(J); // length k = diag((X'X)^{-1})

  for (int j = 0; j < J; j++) {
    NumericMatrix Xj = Xlist[j];           // n x k, column-major
    Eigen::Map<const Eigen::MatrixXd> Xe(&Xj[0], n, k);
    Eigen::MatrixXd X = Xe;                // own copy (thread-safe reads)
    Eigen::MatrixXd XtX = X.transpose() * X;             // k x k
    Eigen::MatrixXd XtXinv = XtX.inverse();              // k x k (full-rank design)
    Xmats[j]   = X;
    Pmats[j]   = XtXinv * X.transpose();                 // k x n
    diagInv[j] = XtXinv.diagonal();                      // k
  }

  // --- outputs (allocated by R; written from threads to disjoint cells) -----
  NumericMatrix Tglob(B, J);               // element(p,j) = Tgptr[p + j*B]
  NumericVector Tpart(B * k * J);          // 3D array (B, k, J)
  Tpart.attr("dim") = IntegerVector::create(B, k, J);
  double* Tgptr = &Tglob[0];
  double* Tpptr = &Tpart[0];

#ifdef _OPENMP
  if (threads < 1) threads = 1;
  #pragma omp parallel for schedule(static) num_threads(threads)
#endif
  for (int p = 0; p < B; p++) {

    Eigen::VectorXd yperm(n), coefv(k), fitted(n);

    for (int j = 0; j < J; j++) {

      const Eigen::MatrixXd& X = Xmats[j];
      const Eigen::MatrixXd& P = Pmats[j];
      const Eigen::VectorXd& d = diagInv[j];
      const double* Ycol = Yptr + static_cast<std::size_t>(j) * n;

      // ---- global test: response = coeff[perm, j] ----
      for (int i = 0; i < n; i++)
        yperm(i) = Ycol[Pptr[p + static_cast<std::size_t>(i) * B] - 1];

      coefv  = P * yperm;
      fitted = X * coefv;
      double rss = (yperm - fitted).squaredNorm();
      double resvar = rss / df;
      double mfit = fitted.mean();
      double sstat = (fitted.array() - mfit).square().sum();
      Tgptr[p + static_cast<std::size_t>(j) * B] = sstat / (nvar * resvar);

      if (method_residuals == 0) {
        // responses method: partial t-stats from the same global-permutation fit
        for (int i = 0; i < k; i++) {
          double se = std::sqrt(d(i) * resvar);
          Tpptr[p + static_cast<std::size_t>(i) * B +
                    static_cast<std::size_t>(j) * B * k] = std::abs(coefv(i)) / se;
        }
      } else {
        // residuals method: reconstruct response per partial test ii
        for (int ii = 0; ii < k; ii++) {
          const double* fcol = fitPtr[ii] + static_cast<std::size_t>(j) * n;
          const double* rcol = resPtr[ii] + static_cast<std::size_t>(j) * n;
          for (int i = 0; i < n; i++)
            yperm(i) = fcol[i] + rcol[Pptr[p + static_cast<std::size_t>(i) * B] - 1];

          coefv  = P * yperm;
          fitted = X * coefv;
          double rss_ii = (yperm - fitted).squaredNorm();
          double resvar_ii = rss_ii / df;
          double se = std::sqrt(d(ii) * resvar_ii);
          Tpptr[p + static_cast<std::size_t>(ii) * B +
                    static_cast<std::size_t>(j) * B * k] = std::abs(coefv(ii)) / se;
        }
      }
    }
  }

  return List::create(Named("T_glob") = Tglob,
                      Named("T_part") = Tpart);
}
