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
  std::vector<double> vbound;
  bool directed;
  bool free_mixture;
  bool free_density;
  std::string distribution;
  std::vector<double> logfactA;
  std::vector<rowvec> pim;
  std::vector<rowvec> pi;
  cube emqr;
  cube nmqr;
  mat alpha;
  std::vector<mat> alpham;
  std::vector<double> delta;
  Rcpp::List vloss;
  // Support masks for free-mixture models (pi, deltapi).
  // Cpi is Q x M: 1 if block q is populated in network m.
  // Calpha is Q x Q: 1 if the pair of blocks is populated in at least one
  // network (the union of the per-network outer products of Cpi columns).
  mat Cpi;
  mat Calpha;

  ColSBM(const std::vector<mat> &A_, int Q_, std::vector<mat> &tau_,
         bool directed_ = false, std::string distribution_ = "bernoulli",
         bool free_mixture_ = true, bool free_density_ = true,
         const mat &Cpi_ = arma::mat(), const mat &Calpha_ = arma::mat());
  ~ColSBM();

  void initialize_state();
  void compute_aggregates();
  void reorder_parameters();
  void compute_vbound(bool store_vloss);
  void update_pi();
  void update_alpha();
  void fixed_point_alpha_delta(double tol = 1e-6, int max_iter = 50);
  mat fixed_point_tau(int m, int max_iter, double tol = TOL);
  double compute_network_vloss(int m) const;
  double compute_network_vloss_params(int m, const mat &tau_m,
                                      const mat &em, const mat &nm) const;
  void step();
  void optimize(int max_step, double tol);
  std::vector<double> get_vbound() const;
};

// Helper to wrap/unwrap from R
Rcpp::XPtr<ColSBM> colsbm_xptr_from_list(const Rcpp::List &A, int Q,
                                         const Rcpp::List &tau, bool directed,
                                         std::string distribution,
                                         bool free_mixture, bool free_density,
                                         const Rcpp::Nullable<Rcpp::LogicalMatrix> &Cpi = R_NilValue,
                                         const Rcpp::Nullable<Rcpp::LogicalMatrix> &Calpha = R_NilValue);

#endif // COLSBM_VEM_H
