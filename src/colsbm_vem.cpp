#include "colsbm_vem.h"

ColSBM::ColSBM(const std::vector<mat> &A_, int Q_) : A(A_), M((int)A_.size()), Q(Q_) {
}

ColSBM::~ColSBM() {}

void ColSBM::optimize(int max_step, double tol) {
  // Minimal placeholder: no-op, set vbound-like quantity via dummy computation
  for (int it = 0; it < max_step; ++it) {
    step();
  }
}

void ColSBM::step() {
  // Placeholder: currently does nothing
}

double ColSBM::get_vbound() const {
  // Placeholder value
  return 0.0;
}

// Convert an R list of matrices to an XPtr<ColSBM>
Rcpp::XPtr<ColSBM> colsbm_xptr_from_list(const Rcpp::List &A, int Q) {
  std::vector<mat> mats;
  mats.reserve(A.size());
  for (auto it = A.begin(); it != A.end(); ++it) {
    Rcpp::RObject obj = *it;
    mat m = Rcpp::as<mat>(obj);
    mats.push_back(m);
  }
  ColSBM *ptr = new ColSBM(mats, Q);
  Rcpp::XPtr<ColSBM> xp(ptr, true);
  return xp;
}
