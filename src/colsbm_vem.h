// Minimal header for colSBM C++ VEM scaffold
#ifndef COLSBM_VEM_H
#define COLSBM_VEM_H

#include "shared_conf.h"
#include <RcppArmadillo.h>
#include <vector>

using arma::cube;
using arma::mat;
using arma::rowvec;
using arma::uvec;
using arma::vec;

class ColSBM {
public:
  std::vector<mat> A; // adjacency matrices
  std::vector<mat> mask;
  int M;
  int Q;
  std::vector<mat> tau;
  int iterations;
  double vbound;
  bool directed;
  bool free_mixture;
  bool free_density;
  std::string distribution;
  std::vector<double> logfactA;
  std::vector<rowvec> pim;
  cube emqr;
  cube nmqr;
  mat alpha;
  std::vector<double> delta;
  std::vector<double> vloss;
  // TODO USE Will be used later to implement free mixture computation
  // this is the support of size QxM indicating which block q
  // is populated by network m
  std::vector<mat> Cpi;
  std::vector<mat> Calpha;

  ColSBM(const std::vector<mat> &A_, int Q_, std::vector<mat> &tau_,
         bool directed_ = false, std::string distribution_ = "bernoulli",
         bool free_mixture_ = true, bool free_density_ = true);
  ~ColSBM();

  void initialize_state();
  void compute_aggregates();
  void update_pi();
  void update_alpha();
  mat fixed_point_tau(int m, int max_iter, double tol = TOL);
  double compute_network_vloss(int m) const;
  void step();
  void optimize(int max_step, double tol);
  double get_vbound() const;
};

// Helper to wrap/unwrap from R
Rcpp::XPtr<ColSBM> colsbm_xptr_from_list(const Rcpp::List &A, int Q,
                                         const Rcpp::List &tau, bool directed,
                                         std::string distribution,
                                         bool free_mixture, bool free_density);

#endif // COLSBM_VEM_H
