// Minimal header for colSBM C++ VEM scaffold
#ifndef COLSBM_VEM_H
#define COLSBM_VEM_H

#include <RcppArmadillo.h>
#include <vector>
#include "shared_conf.h"

using arma::mat;
using arma::vec;
using arma::rowvec;
using arma::uvec;

class ColSBM {
public:
  std::vector<mat> A; // adjacency matrices
  std::vector<mat> mask;
  int M;
  int Q;
  int iterations;
  double vbound;
  bool directed;
  bool free_mixture;
  bool free_density;
  std::string distribution;
  std::vector<double> logfactA;
  std::vector<mat> tau;
  std::vector<rowvec> pim;
  std::vector<mat> emqr;
  std::vector<mat> nmqr;
  std::vector<mat> alpham;
  std::vector<double> delta;
  std::vector<double> vloss;
  // TODO USE Will be used later to implement free mixture computation
  // this is the support of size QxM indicating which block q 
  // is populated by network m
  std::vector<mat> Cpi;
  std::vector<mat> Calpha;


  ColSBM(const std::vector<mat> &A_, int Q_, bool directed_ = false, std::string distribution_ = "bernoulli", bool free_mixture_ = true, bool free_density_ = true);
  ~ColSBM();

  void initialize_state();
  void compute_aggregates();
  void update_pi();
  void update_alpha();
  mat fixed_point_tau(int m, int max_iter = 1000, double tol = 1e-9);
  double compute_network_vloss(int m) const;
  void step();
  void optimize(int max_step, double tol);
  double get_vbound() const;
};

// Helper to wrap/unwrap from R
Rcpp::XPtr<ColSBM> colsbm_xptr_from_list(const Rcpp::List &A, int Q);

#endif // COLSBM_VEM_H
