// [[Rcpp::depends(RcppArmadillo)]]
#include "colsbm_vem.h"
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
XPtr<ColSBM> colsbm_create(List A, int Q) {
  return colsbm_xptr_from_list(A, Q);
}

// [[Rcpp::export]]
void colsbm_optimize(XPtr<ColSBM> ptr, int max_step = 100L, double tol = 1e-6) {
  if (!ptr)
    stop("NULL pointer");
  ptr->optimize(max_step, tol);
}

// [[Rcpp::export]]
double colsbm_vbound(XPtr<ColSBM> ptr) {
  if (!ptr)
    stop("NULL pointer");
  return ptr->get_vbound();
}

// [[Rcpp::export]]
List colsbm_info(XPtr<ColSBM> ptr) {
  if (!ptr)
    stop("NULL pointer");
  List out;
  out["M"] = ptr->M;
  out["Q"] = ptr->Q;
  out["vbound"] = ptr->get_vbound();
  out["iterations"] = ptr->iterations;
  List mask_out(ptr->mask.size());
  for (std::size_t i = 0; i < ptr->mask.size(); ++i) {
    mask_out[i] = wrap(ptr->mask[i]);
  }
  out["mask"] = mask_out;
  List tau_out(ptr->tau.size());
  for (std::size_t i = 0; i < ptr->tau.size(); ++i) {
    tau_out[i] = wrap(ptr->tau[i]);
  }
  out["tau"] = tau_out;

  return out;
}
