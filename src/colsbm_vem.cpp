#include "colsbm_vem.h"
#include "utils.h"
#include <cmath>
#include <string>

ColSBM::ColSBM(const std::vector<mat> &A_, int Q_, std::vector<mat> &tau_,
               bool directed_, std::string distribution_, bool free_mixture_,
               bool free_density_, const mat &Cpi_, const mat &Calpha_)
    : A(A_), M((int)A_.size()), Q(Q_), tau(tau_), iterations(0), vbound(0.0),
      directed(directed_), free_mixture(free_mixture_),
      free_density(free_density_), distribution(std::move(distribution_)),
      Cpi(Cpi_), Calpha(Calpha_) {

  initialize_state();
}

ColSBM::~ColSBM() {}

void ColSBM::initialize_state() {
  mask.clear();
  pim.clear();
  pi.clear();
  emqr.clear();
  nmqr.clear();
  delta.clear();
  logfactA.clear();
  alpham.clear();

  // Checking if distribution is implemented
  check_emission_distribution_unipartite(distribution);

  // Support masks, mirroring the R constructor: non-free-mixture models or
  // missing masks force full support. When Calpha is not provided it is the
  // union of the per-network outer products of the Cpi columns.
  if (!free_mixture || Cpi.n_elem == 0) {
    Cpi = arma::ones<mat>(Q, M);
    Calpha = arma::ones<mat>(Q, Q);
  } else {
    Cpi = arma::conv_to<mat>::from(Cpi > 0.5);
    if (Calpha.n_elem == 0) {
      Calpha = arma::zeros<mat>(Q, Q);
      for (int m = 0; m < M; ++m) {
        Calpha += Cpi.col(m) * Cpi.col(m).t();
      }
      Calpha = arma::conv_to<mat>::from(Calpha > 0.5);
    } else {
      Calpha = arma::conv_to<mat>::from(Calpha > 0.5);
    }
  }

  // Enforce the support on the initial memberships, mirroring R's
  // init_clust("empty") which zeroes tau[, !Cpi[, m]] and renormalizes the
  // remaining columns over the support before the VEM. For unconstrained
  // models (Cpi all TRUE) this is a no-op.
  for (int m = 0; m < M; ++m) {
    const arma::uvec supp = arma::find(Cpi.col(m) > 0.5);
    const arma::uvec not_supp = arma::find(Cpi.col(m) < 0.5);
    if (not_supp.n_elem > 0) {
      tau[m].cols(not_supp).zeros();
      const mat tau_supp = tau[m].cols(supp);
      tau[m].cols(supp) =
          tau_supp / arma::repmat(arma::sum(tau_supp, 1), 1, tau_supp.n_cols);
    }
  }

  for (int m = 0; m < M; ++m) {
    // There we initialize the mask to mask the  diagonal at least
    const int n = A[m].n_rows;
    mat mask_m = arma::zeros<mat>(n, n);
    mask_m.diag().ones();
    mask_m.elem(find_nonfinite(A[m])).ones();
    mask.push_back(mask_m);

    // And we set the NA in A to a value to allow computations
    A[m].elem(find_nonfinite(A[m])).fill(NA_REPLACE_VALUE);

    // Initialize alpham
    mat alpha_m(Q, Q, arma::fill::value(1.0 / static_cast<double>(Q)));
    alpham.push_back(alpha_m);

    // Equiprobability for pis
    rowvec pi_m(Q, arma::fill::value(1.0 / static_cast<double>(Q)));
    pim.push_back(pi_m);
    pi.push_back(pi_m);

    delta.push_back(1.0);

    vloss.push_back(std::vector<double>{});

    if (distribution == "poisson") {
      double logfact = 0.0;
      for (arma::uword i = 0; i < A[m].n_rows; ++i) {
        for (arma::uword j = 0; j < A[m].n_cols; ++j) {
          if (mask_m(i, j) < 0.5) {
            logfact += std::lgamma(A[m](i, j) + 1.0);
          }
        }
      }
      logfactA.push_back(logfact);
    } else {
      logfactA.push_back(0.0);
    }
  }

  emqr.zeros(Q, Q, M);
  nmqr.zeros(Q, Q, M);
  compute_aggregates();

  alpha.ones(Q, Q);
  update_alpha();

  update_pi();
  reorder_parameters();
  update_alpha();

  compute_vbound(true);
  iterations = 1;
}

void ColSBM::compute_aggregates() {
  for (int m = 0; m < M; ++m) {
    const mat &X = A[m];
    const mat &tau_m = tau[m];
    const mat tau_t = tau_m.t();
    const mat Um = 1.0 - mask[m];

    mat emqr_m = tau_t * (X % Um) * tau_m;
    mat nmqr_m = tau_t * Um * tau_m;

    // Mask out the blocks that are empty in network m, as R's update_mqr.
    const mat mask_m = Cpi.col(m) * Cpi.col(m).t();
    emqr.slice(m) = mask_m % emqr_m;
    nmqr.slice(m) = mask_m % nmqr_m;
  }
}

void ColSBM::update_pi() {
  for (int m = 0; m < M; ++m) {
    // Get tau_m
    const mat &tau_m = tau[m];
    rowvec pi_m(Q, arma::fill::zeros);
    for (int q = 0; q < Q; ++q) {
      pi_m[q] = Cpi(q, m) * arma::mean(tau_m.col(q));
    }
    pim[m] = pi_m;
  }

  if (!free_mixture) {
    // If in iid the pi != pim, but weighted_mean(pim, w = n_m)
    rowvec mean_pi(Q, arma::fill::zeros);
    int total_nb_nodes{0};
    for (int m = 0; m < M; ++m) {
      total_nb_nodes += A[m].n_rows;
    }
    if (DEBUG_VE) {
      Rcpp::Rcout << "total nb nodes : " << total_nb_nodes << "\n";
    }
    for (int m = 0; m < M; ++m) {
      if (DEBUG_VE) {
        Rcpp::Rcout << "pim (" << m << ") : " << pim[m] << "\n";
        Rcpp::Rcout << "weighted pi to add : "
                    << ((float)A[m].n_rows / (float)total_nb_nodes) * pim[m]
                    << "\n";
      }
      mean_pi += ((float)A[m].n_rows / (float)total_nb_nodes) * pim[m];
      if (DEBUG_VE) {
        Rcpp::Rcout << "mean pi : " << mean_pi << "\n";
      }
    }

    mean_pi = mean_pi / arma::sum(mean_pi);

    if (DEBUG_VE) {
      Rcpp::Rcout << "mean pi : " << mean_pi
                  << "and sum : " << arma::sum(mean_pi) << "\n";
    }
    for (int m = 0; m < M; ++m) {
      pi[m] = mean_pi;
    }
  } else {
    // If we are in a case of free mixture the pi = pim
    pi = pim;
  }
}

void ColSBM::update_alpha() {
  if (DEBUG_VE) {
    Rcpp::Rcout << "sum(EMQR) = " << arma::sum(emqr, 2)
                << " | sum(NMQR) = " << arma::sum(nmqr, 2) << "\n";
  }
  if (distribution == "bernoulli") {
    for (int m = 0; m < M; ++m) {
      const mat mask_m = Cpi.col(m) * Cpi.col(m).t();
      mat ratio = emqr.slice(m) / nmqr.slice(m);
      // Clamp, then map NaN (empty blocks) to 0, then mask, mirroring R's
      // update_alpham where empty blocks yield exactly 0.
      ratio = clamp_matrix(ratio, TOL, 1.0);
      ratio = nan_to_zero(ratio);
      alpham[m] = mask_m % ratio;
    }
    if (free_density) {
      // Delta model for Bernoulli: same fixed-point on the rates and network
      // densities as the (approx-Poisson) R path, then clamp the rates back
      // into the probability range.
      fixed_point_alpha_delta();
      alpha = Calpha % clamp_matrix(alpha, TOL, 1.0);
    } else {
      mat ratio = arma::sum(emqr, 2) / arma::sum(nmqr, 2);
      ratio = nan_to_zero(ratio);
      ratio = clamp_matrix(ratio, TOL, 1.0);
      alpha = Calpha % ratio;
    }
  } else {
    // Poisson: per-network rates (used for block reordering) are the
    // empirical rates emqr[m] / nmqr[m] with NaN mapped to 0 as in R.
    for (int m = 0; m < M; ++m) {
      const mat mask_m = Cpi.col(m) * Cpi.col(m).t();
      mat ratio = emqr.slice(m) / nmqr.slice(m);
      ratio = nan_to_zero(ratio);
      alpham[m] = mask_m % ratio;
    }
    if (free_density) {
      // The delta-model M step: alternate updates of the shared rates alpha
      // and the network densities delta until convergence.
      fixed_point_alpha_delta();
      alpha = Calpha % alpha;
    } else {
      // iid: pooled rates over all networks, NaN mapped to 0 as in R.
      alpha = Calpha %
              nan_to_zero(arma::sum(emqr, 2) / arma::sum(nmqr, 2));
    }
  }

  if (!std::isfinite(arma::accu(alpha))) {
    Rcpp::warning("Some values of alpha were not finite !");
    if (DEBUG_VE) {
      Rcpp::Rcout << "Some values of alpha were not finite !\n";
    }
  }
}

// Fixed-point M step for the delta models (Poisson with free_density, and
// Bernoulli with approx_pois in R). Mirrors R's fixed_point_alpha_delta:
//   delta[m] = sum(emqr[m]) / sum(alpha * nmqr[m])   (delta[1] = 1)
//   alpha    = sum_m emqr[m]    / sum_m delta[m] * nmqr[m]
// iterated until mean((a-a_old)^2) + mean((d-d_old)^2) <= 2*tol.
void ColSBM::fixed_point_alpha_delta(double tol, int max_iter) {
  mat a = alpha;
  vec d(M);
  for (int m = 0; m < M; ++m) {
    d[m] = delta[m];
  }

  for (int it = 0; it < max_iter; ++it) {
    const mat a_old = a;
    const vec d_old = d;

    for (int m = 0; m < M; ++m) {
      const double num = arma::accu(emqr.slice(m));
      const double den = arma::accu(a % nmqr.slice(m));
      d[m] = (den > 0.0) ? num / den : 0.0;
    }
    d[0] = 1.0;

    mat num(Q, Q, arma::fill::zeros);
    mat den(Q, Q, arma::fill::zeros);
    for (int m = 0; m < M; ++m) {
      num += emqr.slice(m);
      den += d[m] * nmqr.slice(m);
    }
    a = nan_to_zero(num / den);

    const double diff = arma::mean(arma::mean(arma::square(a - a_old))) +
                        arma::mean(arma::square(d - d_old));
    if (diff <= 2.0 * tol) {
      break;
    }
  }

  alpha = a;
  for (int m = 0; m < M; ++m) {
    delta[m] = d[m];
  }
}

// The VE step
mat ColSBM::fixed_point_tau(int m, int max_iter = 1, double tol) {
  mat tau_old = tau[m];
  const mat Um = 1.0 - mask[m];
  const mat Am = A[m] % Um;

  // Support of network m: only these blocks can receive variational mass.
  const arma::uvec supp = arma::find(Cpi.col(m) > 0.5);
  const arma::uvec not_supp = arma::find(Cpi.col(m) < 0.5);
  // As in R: tau_old <- t(Cpi[, m] * t(tau[[m]]))
  tau_old.cols(not_supp).zeros();

  double prev_vloss = compute_network_vloss(m);

  for (int it = 0; it < max_iter; ++it) {
    mat tau_new(tau_old.n_rows, tau_old.n_cols, arma::fill::zeros);
    const mat log_pi = arma::repmat(log_clamped(pi[m]), tau_old.n_rows, 1);

    mat score;
    if (distribution == "bernoulli") {
      const mat logit_alpha =
          log_clamped(alpha * delta[m] / (1.0 - alpha * delta[m]));
      const mat log_1m_alpha = log_clamped(1.0 - alpha * delta[m]);
      score = log_pi + Am * tau_old * logit_alpha.t() +
              Um * tau_old * log_1m_alpha.t();
      if (directed) {
        score +=
            Am.t() * tau_old * logit_alpha + Um.t() * tau_old * log_1m_alpha;
      }
    } else { // poisson
      const mat lambda = alpha * delta[m];
      const mat log_lambda = log_clamped(lambda, POISSON_EPS);
      score =
          log_pi + Am * tau_old * log_lambda.t() - Um * tau_old * lambda.t();
      if (directed) {
        score += Am.t() * tau_old * log_lambda - Um.t() * tau_old * lambda;
      }
    }

    // Mask then softmax over the support only, as in R's fixed_point_tau:
    // tau_new[, !Cpi] <- 0 ; tau_new[, Cpi] <- softmax(...) ; clamp ;
    // renormalize over the support ; zero the empty blocks again.
    mat tau_supp = score.cols(supp);
    tau_supp = softmax_rows(tau_supp);
    tau_supp = clamp_matrix(tau_supp, TAU_EPS, 1 - TAU_EPS);
    // To renormalize after clamping
    tau_supp =
        tau_supp / arma::repmat(arma::sum(tau_supp, 1), 1, tau_supp.n_cols);
    tau_new.cols(supp) = tau_supp;
    tau_new.cols(not_supp).zeros();

    // Recompute the per-network aggregates from the candidate tau and accept
    // the update only if the network vloss improves, mirroring R's
    // fixed_point_tau. This keeps the VE step monotone and avoids the
    // oscillation observed in the undirected Bernoulli case.
    const mat mask_m = Cpi.col(m) * Cpi.col(m).t();
    mat emqr_cand = mask_m % (tau_new.t() * Am * tau_new);
    mat nmqr_cand = mask_m % (tau_new.t() * Um * tau_new);
    const double new_vloss =
        compute_network_vloss_params(m, tau_new, emqr_cand, nmqr_cand);

    if (new_vloss <= prev_vloss) {
      // No improvement: keep the previous tau and stop. Still enforce the
      // support constraint: forbidden blocks are masked out of the vloss, so
      // zeroing them here does not change the objective (R's fixed_point_tau
      // also always keeps tau_new[, !Cpi] zeroed).
      tau[m].cols(not_supp).zeros();
      break;
    }

    tau[m] = tau_new;
    const double delta_tau = arma::accu(arma::square(tau_new - tau_old)) /
                             static_cast<double>(tau_new.n_elem);
    tau_old = tau_new;
    prev_vloss = new_vloss;

    if (delta_tau <= tol) {
      break;
    }
  }

  return tau[m];
}

void ColSBM::reorder_parameters() {
  uvec ord_m;
  for (int m = 0; m < M; ++m) {
    ord_m = arma::sort_index(alpham[m].diag(), "descend");
    if (DEBUG_VE) {
      Rcpp::Rcout << "Reordering network " << m << " order is ";
      ord_m.print(Rcpp::Rcout);
      Rcpp::Rcout << "\n";
    }

    // Sorting all related objects

    if (DEBUG_VE) {
      Rcpp::Rcout << "Tau will be\n";
      tau[m].cols(ord_m).brief_print(Rcpp::Rcout);
      Rcpp::Rcout << "\n";
    }

    tau[m] = tau[m].cols(ord_m);

    // Retrieving slice m as a matrix
    mat emqr_m = emqr.slice(m);
    emqr_m = emqr_m.rows(ord_m);
    emqr_m = emqr_m.cols(ord_m);
    if (DEBUG_VE) {
      Rcpp::Rcout << "emqr will be\n";
      emqr_m.brief_print(Rcpp::Rcout);
      Rcpp::Rcout << "\n";
    }
    // Restoring the slice from the permutated matrix
    emqr.slice(m) = emqr_m;

    // Retrieving slice m as a matrix
    mat nmqr_m = nmqr.slice(m);
    nmqr_m = nmqr_m.rows(ord_m);
    nmqr_m = nmqr_m.cols(ord_m);
    if (DEBUG_VE) {
      Rcpp::Rcout << "nmqr will be\n";
      nmqr_m.brief_print(Rcpp::Rcout);
      Rcpp::Rcout << "\n";
    }
    // Restoring the slice from the permutated matrix
    nmqr.slice(m) = nmqr_m;

    if (DEBUG_VE) {
      Rcpp::Rcout << "pi is\n";
      pim[m].brief_print(Rcpp::Rcout);
      Rcpp::Rcout << "pi will be\n";
      pim[m](ord_m).t().brief_print(Rcpp::Rcout);
      Rcpp::Rcout << "\n";
    }
    // Turning back pim to rowvec
    pim[m] = pim[m](ord_m).t();

    // Permute the support of network m the same way, as R's
    // reorder_parameters: Cpi[, m] <- Cpi[ord_m, m]
    vec cpi_m = Cpi.col(m);
    Cpi.col(m) = cpi_m(ord_m);
  }

  // Recompute Calpha from the reordered Cpi so it stays consistent with the
  // new block ordering (R permutes Cpi then relies on Calpha; recomputing is
  // the identity on the full-support case and correct otherwise).
  Calpha = arma::zeros<mat>(Q, Q);
  for (int m = 0; m < M; ++m) {
    Calpha += Cpi.col(m) * Cpi.col(m).t();
  }
  Calpha = arma::conv_to<mat>::from(Calpha > 0.5);
}

double ColSBM::compute_network_vloss_params(int m, const mat &tau_m,
                                            const mat &em,
                                            const mat &nm) const {
  const double delta_m = delta[m];
  const rowvec pi_m = pi[m];
  // Support of network m: pi and entropy terms are restricted to it, and the
  // alpha terms only count supported pairs (mirroring R's vb_tau_alpha,
  // vb_tau_pi and entropy_tau).
  const arma::uvec supp = arma::find(Cpi.col(m) > 0.5);
  const mat mask_m = Cpi.col(m) * Cpi.col(m).t();

  double obj = 0.0;
  double dircoef = directed ? 1 : 0.5;

  // Vbound tau and alpha
  if (distribution == "bernoulli") {
    obj += dircoef *
           arma::accu(Calpha %
                      (mask_m % em % log_clamped(delta_m * alpha, POISSON_EPS) +
                       mask_m % (nm - em) %
                           log_clamped(1 - delta_m * alpha, POISSON_EPS)));
  } else if (distribution == "poisson") {
    obj += dircoef *
           (arma::accu(Calpha % (mask_m % em) %
                       log_clamped(delta_m * alpha, POISSON_EPS)) -
            arma::accu(Calpha % (mask_m % nm) % (delta_m * alpha)) -
            logfactA[m]);
  }
  // Vbound tau and pi (plain log as in R's vb_tau_pi)
  obj += arma::accu(tau_m.cols(supp) * arma::log(pi_m.cols(supp).t()));

  // Vbound entropy (R's entropy_tau uses x*log(x), zero when x is null)
  const mat tau_supp = tau_m.cols(supp);
  const mat log_tau =
      arma::log(arma::clamp(tau_supp, 1e-300, arma::datum::inf));
  obj += -arma::accu(tau_supp % log_tau);

  return obj;
}

double ColSBM::compute_network_vloss(int m) const {
  return compute_network_vloss_params(m, tau[m], emqr.slice(m), nmqr.slice(m));
}

void ColSBM::step() {
  // VE step
  for (int m = 0; m < M; ++m) {
    fixed_point_tau(m, 1, TOL);
  }
  compute_aggregates();
  if (DEBUG_VE) {
    Rcpp::Rcout << "Before reordering\n";
    for (int m = 0; m < M; ++m) {
      tau[m].brief_print(Rcpp::Rcout);
      alpham[m].brief_print(Rcpp::Rcout);
      pim[m].brief_print(Rcpp::Rcout);
    }
  }
  reorder_parameters();
  if (DEBUG_VE) {
    Rcpp::Rcout << "After reordering\n";
    for (int m = 0; m < M; ++m) {
      tau[m].brief_print(Rcpp::Rcout);
      alpham[m].brief_print(Rcpp::Rcout);
      pim[m].brief_print(Rcpp::Rcout);
    }
  }
  // M step
  update_pi();
  update_alpha();
}

// A function to loop over all networks and store the vbound
void ColSBM::compute_vbound(bool store_vloss = false) {
  double current_vbound = 0.0;
  for (int m = 0; m < M; ++m) {
    double vloss_m = compute_network_vloss(m);
    if (DEBUG_VE) {
      Rcpp::Rcout << "vloss (" << m << ") : " << vloss_m << "\n";
    }
    if (store_vloss) {
      // Casting the m element of the list as a NumericVector to allow pushing
      // the loss value to it
      // Rcpp::as<Rcpp::NumericVector>(vloss[m]);
      std::vector<double> vloss_m_values =
          Rcpp::as<std::vector<double>>(vloss[m]);
      vloss_m_values.push_back((double)vloss_m);
      vloss[m] = vloss_m_values;
      if (DEBUG_VE) {
        Rcpp::Rcout << "vloss[" << m
                    << "] = " << Rcpp::as<std::vector<double>>(vloss[m]).back()
                    << "\n";
      }
    }
    current_vbound += vloss_m;
  }
  vbound.push_back(current_vbound);
}

void ColSBM::optimize(int max_step, double tol = VBOUND_TOL) {
  if (Q == 1) {
    // In this case we only need to return the computed parameters at this step
    // and fix the tau parameters
    for (int m = 0; m < M; ++m) {
      tau[m] = arma::ones(A[m].n_rows, 1);
    }
    update_alpha();
  } else {
    // Compute the first vbound
    double prev_vbound = vbound.back();

    // max_step - 1 to account for the first M step after provided taus
    for (int it = 0; it < max_step - 1; ++it) {
      if (DEBUG_VE) {
        Rcpp::Rcout << "Iteration " << it << "\n";
      }

      // Make VE and M step
      step();

      // Update the vbound after step
      prev_vbound = vbound.back();
      compute_vbound(true);
      ++iterations;
      if (DEBUG_LOOP) {
        Rcpp::Rcout << "Iteration " << iterations
                    << " |Delta vb| = " << std::abs(vbound.back() - prev_vbound)
                    << " | tol = " << tol << "\n";
      }
      if (std::abs(vbound.back() - prev_vbound) < tol) {
        break;
      }
    }
    if (iterations > max_step) {
      Rcpp::warning("The VEM failed to converge in %i max steps for Q=%i !",
                    max_step, Q);
    }
  }
}

std::vector<double> ColSBM::get_vbound() const { return vbound; }

// Convert an R list of matrices to an XPtr<ColSBM>
Rcpp::XPtr<ColSBM> colsbm_xptr_from_list(
    const Rcpp::List &A, int Q, const Rcpp::List &tau, bool directed,
    std::string distribution, bool free_mixture, bool free_density,
    const Rcpp::Nullable<Rcpp::LogicalMatrix> &Cpi,
    const Rcpp::Nullable<Rcpp::LogicalMatrix> &Calpha) {
  std::vector<mat> mats;
  mats.reserve(A.size());
  for (Rcpp::List::const_iterator it = A.begin(); it != A.end(); ++it) {
    mat m = Rcpp::as<mat>(*it);
    mats.push_back(m);
  }
  std::vector<mat> taus;
  taus.reserve(tau.size());
  for (Rcpp::List::const_iterator it = tau.begin(); it != tau.end(); ++it) {
    mat m = Rcpp::as<mat>(*it);
    m = clamp_matrix(m);
    taus.push_back(m);
  }

  // Support masks: empty matrix marks "not provided" -> full support.
  mat Cpi_m, Calpha_m;
  if (Cpi.isNotNull()) {
    const Rcpp::LogicalMatrix Cp = Rcpp::as<Rcpp::LogicalMatrix>(Cpi);
    Cpi_m = arma::mat(Cp.nrow(), Cp.ncol(), arma::fill::zeros);
    for (int j = 0; j < Cp.ncol(); ++j) {
      for (int i = 0; i < Cp.nrow(); ++i) {
        Cpi_m(i, j) = (Cp(i, j) == TRUE) ? 1.0 : 0.0;
      }
    }
  }
  if (Calpha.isNotNull()) {
    const Rcpp::LogicalMatrix Ca = Rcpp::as<Rcpp::LogicalMatrix>(Calpha);
    Calpha_m = arma::mat(Ca.nrow(), Ca.ncol(), arma::fill::zeros);
    for (int j = 0; j < Ca.ncol(); ++j) {
      for (int i = 0; i < Ca.nrow(); ++i) {
        Calpha_m(i, j) = (Ca(i, j) == TRUE) ? 1.0 : 0.0;
      }
    }
  }

  ColSBM *ptr = new ColSBM(
      mats, Q, taus, directed = directed, distribution = distribution,
      free_mixture = free_mixture, free_density = free_density, Cpi_m,
      Calpha_m);
  Rcpp::XPtr<ColSBM> xp(ptr, true);
  return xp;
}
