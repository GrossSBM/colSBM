// [[Rcpp::depends(RcppArmadillo)]]
#include "colsbm_vem.h"
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
XPtr<ColSBM> colsbm_create(List A, int Q, List tau, bool directed = false,
                           std::string distribution = "bernoulli",
                           bool free_mixture = true, bool free_density = true) {
  return colsbm_xptr_from_list(
      A, Q, tau, directed = directed, distribution = distribution,
      free_mixture = free_mixture, free_density = free_density);
}

// [[Rcpp::export]]
void colsbm_optimize(XPtr<ColSBM> ptr, int max_step = 100, double tol = 1e-6) {
  if (!ptr)
    stop("NULL pointer");
  ptr->optimize(max_step, tol);
}

// [[Rcpp::export]]
std::vector<double> colsbm_vbound(XPtr<ColSBM> ptr) {
  if (!ptr)
    stop("NULL pointer");
  return ptr->get_vbound();
}

// [[Rcpp::export]]
List colsbm_info(XPtr<ColSBM> ptr) {
  if (!ptr)
    stop("NULL pointer");
  List out;
  List A_out(ptr->A.size());
  for (std::size_t i = 0; i < ptr->A.size(); ++i) {
    mat Ai = ptr->A[i];
    // Replacing the masked constant by NaN
    Ai.elem(find(ptr->mask[i] == 1)).fill(NA_REAL);
    A_out[i] = wrap(Ai);
  }
  out["A"] = A_out;
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
  out["emqr"] = ptr->emqr;
  out["nmqr"] = ptr->nmqr;
  out["alpha"] = ptr->alpha;
  out["delta"] = ptr->delta;
  out["pim"] = ptr->pim;
  out["pi"] = ptr->pi;
  out["vloss"] = ptr->vloss;

  return out;
}
