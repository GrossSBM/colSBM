test_that("R with C++ yields the same results as the R procedure on modular networks", {
  M <- 4L
  true_pi <- c(0.6, 0.3, 0.1)
  nbNodes <- 1000
  K <- length(true_pi)
  true_alpha <- diag(c(0.9, 0.5, 0.1))
  sbm_sampler <- sbm::sampleSimpleSBM(nbNodes = nbNodes, blockProp = true_pi, connectParam = list(mean = true_alpha), directed = TRUE)
  set.seed(1234)
  sbm_realisations <- lapply(seq(M), sbm_sampler[["rNetwork"]])
  effective_pim <- lapply(sbm_realisations, function(sbm_real) sbm_real[["indMemberships"]] |> colMeans())
  effective_pi <- rowMeans(simplify2array(effective_pim))
  A <- lapply(seq(M), function(m) sbm_realisations[[m]][["networkData"]])

  # A spectral clustering for the start
  spectral_labels <- lapply(A, spectral_clustering, K)

  spectral_taus <- lapply(spectral_labels, .one_hot, K)

  spectral_taus <- lapply(spectral_taus, function(tau) {
    tau[tau < 1e-6] <- 1e-6
    tau[tau > 1 - 1e-6] <- 1 - 1e-6
    tau <- tau / rowSums(tau)
  })


  R_fit <- estimate_colSBM(netlist = A, colsbm_model = "iid", fit_opts = list(use_cpp = FALSE, verbosity = 2L))

  RCpp_fit <- estimate_colSBM(netlist = A, colsbm_model = "iid", fit_opts = list(use_cpp = TRUE, verbosity = 2L), global_opts = list(backend = "no_mc", verbosity = 2L), nb_run = 1L)
})

test_that("R with C++ yields the same results as the R procedure on core-periphery networks", {
  M <- 4
  true_pi <- c(0.6, 0.4)
  nbNodes <- 1000
  K <- length(true_pi)
  true_alpha <- matrix(c(
    0.9, 0.7,
    0.35, 0.05
  ), nrow = K, byrow = TRUE)
  sbm_sampler <- sbm::sampleSimpleSBM(nbNodes = nbNodes, blockProp = true_pi, connectParam = list(mean = true_alpha), directed = TRUE)
  set.seed(1234)
  sbm_realisations <- lapply(seq(M), sbm_sampler[["rNetwork"]])
  effective_pim <- lapply(sbm_realisations, function(sbm_real) sbm_real[["indMemberships"]] |> colMeans())
  effective_pi <- rowMeans(simplify2array(effective_pim))
  A <- lapply(seq(M), function(m) sbm_realisations[[m]][["networkData"]])

  # A spectral clustering for the start
  spectral_labels <- lapply(A, spectral_clustering, K)

  spectral_taus <- lapply(spectral_labels, .one_hot, K)

  spectral_taus <- lapply(spectral_taus, function(tau) {
    tau[tau < 1e-6] <- 1e-6
    tau[tau > 1 - 1e-6] <- 1 - 1e-6
    tau <- tau / rowSums(tau)
  })

  # colSBM R logic
  R_fit <- fitSimpleSBMPop$new(A = A, Q = K, Z = spectral_taus, directed = TRUE, free_mixture = FALSE, free_density = FALSE, init_method = "given", fit_opts = default_fit_opts_unipartite())
  R_fit$optimize()

  # C++ logic
  ptr <- colsbm_create(A = A, Q = K, tau = spectral_taus, directed = TRUE, distribution = "bernoulli", free_mixture = FALSE, free_density = FALSE)
  colsbm_optimize(ptr, max_step = 100L, tol = 1e-3)

  info <- colsbm_info(ptr)
  labelize <- function(tau) {
    sapply(seq_len(nrow(tau)), function(i) which.max(tau[i, ]))
  }
  info$Z <- lapply(info$tau, labelize)

  # No incorrect output from C++
  expect_equal(object = rowSums(info$tau[[1]]), expected = rep(1, nrow(info$tau[[1]])), tolerance = 1e-6)
  expect_true(all(info$tau[[1]] >= 0))
  expect_type(info, type = "list")
  expect_identical(info$M, 4L)
  expect_identical(info$Q, K)
  expect_all_true(is.finite(info$vbound))


  # True quality checks

  expect_gte(tail(info$vbound, 1), tail(R_fit$vbound, 1))
  expect_equal(info$tau, R_fit$tau, tolerance = 1e-1)
  ## Parameters checks

  expect_equal(info$pi, rep(list(true_pi), M), tolerance = 1e-2)
  expect_equal(info$pim, R_fit$pim, tolerance = 1e-2)

  ## Check that clusterings coincides
  expect_identical(sapply(seq_along(info$Z), function(m) aricode::ARI(R_fit[["Z"]][[m]], info$Z[[m]])), rep(1, M))
})
