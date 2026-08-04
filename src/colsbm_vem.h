// Minimal header for colSBM C++ VEM scaffold
#ifndef COLSBM_VEM_H
#define COLSBM_VEM_H

#include <RcppArmadillo.h>
#include <vector>

using arma::mat;
using arma::vec;
using arma::rowvec;

class ColSBM {
public:
  std::vector<mat> A; // adjacency matrices
  int M;
  int Q;
  int iterations;
  double vbound;
  std::vector<mat> tau;
  std::vector<rowvec> pi;
  std::vector<mat> emqr;
  std::vector<mat> nmqr;
  std::vector<mat> alpha;
  std::vector<vec> delta;

  ColSBM(const std::vector<mat> &A_, int Q_);
  ~ColSBM();

  void initialize_state();
  void compute_aggregates();
  void update_pi();
  void step();
  void optimize(int max_step, double tol);
  double get_vbound() const;
};

// Helper to wrap/unwrap from R
Rcpp::XPtr<ColSBM> colsbm_xptr_from_list(const Rcpp::List &A, int Q);

#endif // COLSBM_VEM_H
