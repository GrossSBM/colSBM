// Minimal header for colSBM C++ VEM scaffold
#ifndef COLSBM_VEM_H
#define COLSBM_VEM_H

#include <RcppArmadillo.h>
#include <vector>

using arma::mat;

class ColSBM {
public:
  // basic fields
  std::vector<mat> A; // adjacency matrices
  int M;
  int Q;

  ColSBM(const std::vector<mat> &A_, int Q_);
  ~ColSBM();

  // placeholder methods
  void optimize(int max_step, double tol);
  void step();
  double get_vbound() const;
};

// Helper to wrap/unwrap from R
Rcpp::XPtr<ColSBM> colsbm_xptr_from_list(const Rcpp::List &A, int Q);

#endif // COLSBM_VEM_H
