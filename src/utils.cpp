#include <algorithm>
#include <cmath>
#include <RcppArmadillo.h>

constexpr double kClampEps = 1e-10;

inline double clamp_log(double x, double eps = kClampEps) {
  return std::log(std::max(x, eps));
}

inline double clamp_value(double x, double lo = kClampEps, double hi = 1.0) {
  return std::max(lo, std::min(hi, x));
}

inline arma::mat softmax_rows(const arma::mat &x) {
  arma::mat out = x;
  for (arma::uword i = 0; i < out.n_rows; ++i) {
    const double row_max = out.row(i).max();
    out.row(i) -= row_max;
    out.row(i) = arma::exp(out.row(i));
    const double denom = std::max(arma::accu(out.row(i)), 1e-10);
    out.row(i) /= denom;
  }
  return out;
}

inline arma::mat log_clamped(const arma::mat &x, double eps = kClampEps) {
  return arma::log(arma::clamp(x, eps, arma::datum::inf));
}

inline arma::mat clamp_matrix(const arma::mat &x, double lo = kClampEps, double hi = 1-kClampEps) {
  return arma::clamp(x, lo, hi);
}
