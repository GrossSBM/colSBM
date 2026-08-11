#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <string>

inline double clamp_log(double x, double eps = TOL) {
  return std::log(std::max(x, eps));
}

inline double clamp_value(double x, double lo = TOL, double hi = 1.0) {
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

inline arma::mat log_clamped(const arma::mat &x, double eps = TOL) {
  return arma::log(arma::clamp(x, eps, arma::datum::inf));
}
//' @param hi a double specifying the max value to clamp, defaults to 1-TOL
//' for the bernoulli distribution
inline arma::mat clamp_matrix(const arma::mat &x, double lo = TOL,
                              double hi = 1 - TOL) {
  return arma::clamp(x, lo, hi);
}

inline void check_emission_distribution_unipartite(std::string distribution) {
  std::vector<std::string> implemented_distributions{"bernoulli", "poisson"};
  std::string result;

  if (std::find(implemented_distributions.begin(),
                implemented_distributions.end(),
                distribution) == implemented_distributions.end()) {

    for (std::size_t i = 0; i < implemented_distributions.size(); ++i) {
      if (i > 0)
        result += ", ";
      result += implemented_distributions[i];
    }

    Rcpp::stop("Distribution " + distribution +
               " not implemented ! Should be one of " + result + ".");
  }
}
