// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
#include "colsbm_vem.h"

using namespace Rcpp;

// [[Rcpp::export]]
XPtr<ColSBM> colsbm_create(List A, int Q) {
  return colsbm_xptr_from_list(A, Q);
}

// [[Rcpp::export]]
void colsbm_optimize(XPtr<ColSBM> ptr, int max_step = 100L, double tol = 1e-6) {
  if (!ptr) stop("NULL pointer");
  ptr->optimize(max_step, tol);
}

// [[Rcpp::export]]
double colsbm_vbound(XPtr<ColSBM> ptr) {
  if (!ptr) stop("NULL pointer");
  return ptr->get_vbound();
}

// [[Rcpp::export]]
List colsbm_info(XPtr<ColSBM> ptr) {
  if (!ptr) stop("NULL pointer");
  List out;
  out["M"] = ptr->M;
  out["Q"] = ptr->Q;
  out["vbound"] = ptr->get_vbound();
  out["iterations"] = ptr->iterations;
  out["tau"] = ptr->tau;
  return out;
}
