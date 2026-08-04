#include "colsbm_vem.h"
#include "utils.cpp"
#include <algorithm>
#include <cmath>

ColSBM::ColSBM(const std::vector<mat> &A_, int Q_) : A(A_), M((int)A_.size()), Q(Q_), iterations(0), vbound(0.0) {
  initialize_state();
}

ColSBM::~ColSBM() {}

void ColSBM::initialize_state() {
  tau.clear();
  pi.clear();
  emqr.clear();
  nmqr.clear();
  alpha.clear();
  delta.clear();

  for (int m = 0; m < M; ++m) {
    const int n = A[m].n_rows;
    mat tau_m = arma::ones<mat>(n, Q) / static_cast<double>(Q);
    tau.push_back(tau_m);

    rowvec pi_m(Q, arma::fill::value(1.0 / static_cast<double>(Q)));
    pi.push_back(pi_m);

    emqr.push_back(mat(Q, Q, arma::fill::zeros));
    nmqr.push_back(mat(Q, Q, arma::fill::zeros));

    alpha.push_back(mat(Q, Q, arma::fill::value(0.5)));
    delta.push_back(vec(1, arma::fill::value(1.0)));
  }

  vbound = 0.0;
}

void ColSBM::compute_aggregates() {
  for (int m = 0; m < M; ++m) {
    const mat &X = A[m];
    const mat &tau_m = tau[m];
    const mat tau_t = tau_m.t();

    mat emqr_m = tau_t * X * tau_m;
    mat nmqr_m = tau_t * arma::ones<mat>(X.n_rows, X.n_cols) * tau_m;

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
    pi[m] = pi_m;
  }
}

void ColSBM::step() {
  compute_aggregates();
  update_pi();

  double sum_log = 0.0;
  for (int m = 0; m < M; ++m) {
    const mat &em = emqr[m];
    const mat &nm = nmqr[m];
    for (int q = 0; q < Q; ++q) {
      for (int r = 0; r < Q; ++r) {
        sum_log += em(q, r) * clamp_log(std::max(em(q, r), 1e-10)) - nm(q, r) * 0.01;
      }
    }
  }
  vbound = sum_log;
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
