#include "colsbm_vem.h"
#include "utils.cpp"
#include <algorithm>
#include <cmath>


ColSBM::ColSBM(const std::vector<mat> &A_, int Q_, bool directed_, std::string distribution_, bool free_mixture_, bool free_density_)
  : A(A_), M((int)A_.size()), Q(Q_), iterations(0), vbound(0.0), directed(directed_), free_mixture(free_mixture_), free_density(free_density_), distribution(std::move(distribution_)) {
  initialize_state();
}

ColSBM::~ColSBM() {}

void ColSBM::initialize_state() {
  mask.clear();
  tau.clear();
  pim.clear();
  emqr.clear();
  nmqr.clear();
  alpham.clear();
  delta.clear();
  vloss.clear();
  logfactA.clear();

  for (int m = 0; m < M; ++m) {
    // There we initialize the mask to mask the  diagonal at least
    const int n = A[m].n_rows;
    mat mask_m = arma::zeros<mat>(n, n);
    mask_m.diag().ones();
    mask_m.elem(find_nonfinite(A[m])).ones();
    mask.push_back(mask_m);

    // The tau matrix is initialized with equal probability everywhere
    mat tau_m = arma::ones<mat>(n, Q) / static_cast<double>(Q);
    tau.push_back(tau_m);

    // Equiprobability for pis
    rowvec pi_m(Q, arma::fill::value(1.0 / static_cast<double>(Q)));
    pim.push_back(pi_m);

    emqr.push_back(mat(Q, Q, arma::fill::zeros));
    nmqr.push_back(mat(Q, Q, arma::fill::zeros));

    alpham.push_back(mat(Q, Q, arma::fill::value(0.5)));
    delta.push_back(1.0);

    vloss.push_back(-arma::datum::inf);

    double logfact = 0.0;
    for (arma::uword i = 0; i < A[m].n_rows; ++i) {
      for (arma::uword j = 0; j < A[m].n_cols; ++j) {
        if (mask_m(i, j) < 0.5) {
          logfact += std::lgamma(A[m](i, j) + 1.0);
        }
      }
    }
    logfactA.push_back(logfact);
  }

  vbound = 0.0;
}

void ColSBM::compute_aggregates() {
  for (int m = 0; m < M; ++m) {
    const mat &X = A[m];
    const mat &tau_m = tau[m];
    const mat tau_t = tau_m.t();
    const mat Um = 1.0 - mask[m];

    mat emqr_m = tau_t * (X % Um) * tau_m;
    mat nmqr_m = tau_t * Um * tau_m;

    emqr[m] = emqr_m;
    nmqr[m] = nmqr_m;
  }
}

void ColSBM::update_pi() {
  for (int m = 0; m < M; ++m) {
    const mat &tau_m = tau[m];
    rowvec pi_m(Q, arma::fill::zeros);
    for (int q = 0; q < Q; ++q) {
      pi_m[q] = arma::mean(tau_m.col(q));
    }
    const double pi_sum = arma::accu(pi_m);
    if (pi_sum > 0.0) {
      pim[m] = pi_m / pi_sum;
    } else {
      pim[m] = rowvec(Q, arma::fill::value(1.0 / static_cast<double>(Q)));
    }
  }
}

void ColSBM::update_alpha() {
  for (int m = 0; m < M; ++m) {
    const mat denom = arma::clamp(nmqr[m], TOL, arma::datum::inf);
    mat alpha_m;

    if (distribution == "bernoulli"){
      alpha_m = clamp_matrix(emqr[m] / denom, TOL, 1.0);
    } else {
      alpha_m = clamp_matrix(emqr[m] / denom, TOL, arma::datum::inf);
    }

    if (!std::isfinite(arma::accu(alpha_m))) {
      alpha_m.fill(0.5);
    }
    alpham[m] = alpha_m;
  }
}

mat ColSBM::fixed_point_tau(int m, int max_iter, double tol) {
  mat tau_old = tau[m];
  const mat Um = 1.0 - mask[m];
  const mat Am = A[m] % Um;

  double prev_vloss = compute_network_vloss(m);
  vloss[m] = prev_vloss;

  for (int it = 0; it < max_iter; ++it) {
    mat tau_new(tau_old.n_rows, tau_old.n_cols, arma::fill::zeros);
    const mat log_pi = arma::repmat(log_clamped(pim[m]), tau_old.n_rows, 1);
    const mat log_alpha = log_clamped(alpham[m] * delta[m]);
    const mat log_1m_alpha = log_clamped(1.0 - alpham[m] * delta[m]);

    mat score = log_pi + Am * tau_old * log_alpha.t() + Um * tau_old * log_1m_alpha.t();
    if (directed) {
      score += Am.t() * tau_old * log_alpha + Um.t() * tau_old * log_1m_alpha;
    }

    tau_new = softmax_rows(score);
    tau[m] = tau_new;
    compute_aggregates();
    update_pi();

    const double new_vloss = compute_network_vloss(m);
    vloss[m] = new_vloss;

    const double delta_tau = arma::accu(arma::square(tau_new - tau_old)) / static_cast<double>(tau_new.n_elem);
    const double delta_vloss = std::abs(new_vloss - prev_vloss);
    tau_old = tau_new;
    prev_vloss = new_vloss;

    if (delta_tau <= tol || delta_vloss <= tol) {
      break;
    }
  }

  return tau[m];
}

double ColSBM::compute_network_vloss(int m) const {
  
  const mat &tau_m = tau[m];
  const mat &em = emqr[m];
  const mat &nm = nmqr[m];
  const double delta_m = delta[m];
  const mat alpha_m = alpham[m] * delta_m;

  // Vbound tau and alpha

  // Vbound tau and pi
  
  // Vbound entropy
  double obj = 0.0;
  if (distribution == "bernoulli") {
    obj += arma::accu(arma::clamp(em, 0.0, arma::datum::inf) % log_clamped(alpha_m));
    obj += arma::accu(arma::clamp(nm - em, 0.0, arma::datum::inf) % log_clamped(1.0 - alpha_m));
  } else {
    obj += arma::accu(arma::clamp(em, 0.0, arma::datum::inf) % log_clamped(alpha_m));
    obj -= arma::accu(nm % alpha_m);
    obj -= logfactA[m];
  }

  const rowvec pi_m = pim[m];
  for (arma::uword i = 0; i < tau_m.n_rows; ++i) {
    for (arma::uword q = 0; q < tau_m.n_cols; ++q) {
      const double val = tau_m(i, q);
      if (val > 0.0) {
        obj += val * clamp_log(pi_m(q));
        obj -= val * clamp_log(val);
      }
    }
  }
  return obj;
}

void ColSBM::step() {
  for (int m = 0; m < M; ++m) {
    fixed_point_tau(m, 1, 1e-6);
  }
  compute_aggregates();
  update_pi();
  update_alpha();

  vbound = -arma::datum::inf;
  for (int m = 0; m < M; ++m) {
    vbound += compute_network_vloss(m);
  }
  ++iterations;
}

void ColSBM::optimize(int max_step, double tol) {
  iterations = 0;
  vbound = 0.0;
  for (int it = 0; it < max_step; ++it) {
    step();
    if (std::abs(vbound) < tol) {
      break;
    }
  }
}

double ColSBM::get_vbound() const {
  return vbound;
}

// Convert an R list of matrices to an XPtr<ColSBM>
Rcpp::XPtr<ColSBM> colsbm_xptr_from_list(const Rcpp::List &A, int Q) {
  std::vector<mat> mats;
  mats.reserve(A.size());
  for (Rcpp::List::const_iterator it = A.begin(); it != A.end(); ++it) {
    mat m = Rcpp::as<mat>(*it);
    mats.push_back(m);
  }
  ColSBM *ptr = new ColSBM(mats, Q);
  Rcpp::XPtr<ColSBM> xp(ptr, true);
  return xp;
}
