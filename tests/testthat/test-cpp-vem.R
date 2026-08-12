test_that("cpp runs normally and stops on non implemented emission distribution", {
  A <- list(
    matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE),
    matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
  )

  taus <- list(matrix(rep(c(0, 1), nrow(A[[1]])), nrow = 2, byrow = TRUE), matrix(rep(c(1, 0), nrow(A[[1]])), nrow = 2, byrow = TRUE))

  expect_no_error(colsbm_create(A, Q = 2L, tau = taus, distribution = "bernoulli"))
  expect_no_error(colsbm_create(A, Q = 2L, tau = taus, distribution = "poisson"))

  # currently not implemented
  expect_error(colsbm_create(A, Q = 2L, tau = taus, distribution = "gaussian"))

  # typos
  expect_error(colsbm_create(A, Q = 2L, tau = taus, distribution = "bernouilli"))
  expect_error(colsbm_create(A, Q = 2L, tau = taus, distribution = "poiton"))
  expect_error(colsbm_create(A, Q = 2L, tau = taus, distribution = "Poisson"))
})

test_that("cpp VE yields the same results as the R procedure on modular networks", {
  M <- 4
  true_pi <- c(0.1, 0.3, 0.6)
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
  spectral_labels <- lapply(A, spectral_clustering, K = 3)

  spectral_taus <- lapply(spectral_labels, .one_hot, Q = 3)

  spectral_taus <- lapply(spectral_taus, function(tau) {
    tau[tau < 1e-6] <- 1e-6
    tau[tau > 1 - 1e-6] <- 1 - 1e-6
    tau <- tau / rowSums(tau)
  })

  # colSBM R logic
  R_fit <- fitSimpleSBMPop$new(A = A, Q = 3L, Z = spectral_taus, directed = TRUE, free_mixture = FALSE, free_density = FALSE, init_method = "given", fit_opts = default_fit_opts_unipartite())
  R_fit$optimize()

  # C++ logic
  ptr <- colsbm_create(A = A, Q = 3L, tau = spectral_taus, directed = TRUE, distribution = "bernoulli", free_mixture = FALSE, free_density = FALSE)
  colsbm_optimize(ptr, max_step = 100L, tol = 1e-9)

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
  expect_identical(info$Q, 3L)
  expect_all_true(is.finite(info$vbound))


  # True quality checks
  cpp_order <- order(diag(info$alpha), decreasing = TRUE)
  expect_equal(info$vbound, R_fit$vbound, tolerance = 1e-2)
  cpp_taus <- lapply(info$tau, function(tau) tau[, cpp_order])
  expect_equal(R_fit$tau, cpp_taus, tolerance = 1e-6)
  ## Parameters checks

  cpp_pi <- lapply(info$pi, function(pi) pi[cpp_order])
  cpp_pim <- lapply(info$pim, function(pim) pim[cpp_order])
  expect_equal(cpp_pi, R_fit$pi, tolerance = 1e-2)
  expect_equal(cpp_pim, R_fit$pim, tolerance = 1e-2)

  ## Check that clusterings coincides
  expect_identical(sapply(seq_along(info$Z), function(m) aricode::ARI(R_fit$Z[[m]], info$Z[[m]])), rep(1, M))
})

test_that("cpp VE yields the same results as the R procedure on core-periphery networks", {
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
  colsbm_optimize(ptr, max_step = 100L, tol = 1e-9)

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
  cpp_order <- order(diag(info$alpha), decreasing = TRUE)

  expect_gte(tail(info$vbound, 1), tail(R_fit$vbound, 1))
  cpp_taus <- lapply(info$tau, function(tau) tau[, cpp_order])
  expect_equal(R_fit$tau, cpp_taus, tolerance = 1e-1)
  ## Parameters checks

  cpp_pi <- lapply(info$pi, function(pi) pi[cpp_order])
  cpp_pim <- lapply(info$pim, function(pim) pim[cpp_order])
  expect_equal(cpp_pi, rep(list(true_pi), M), tolerance = 1e-2)
  expect_equal(cpp_pim, R_fit$pim, tolerance = 1e-2)

  ## Check that clusterings coincides
  expect_identical(sapply(seq_along(info$Z), function(m) aricode::ARI(R_fit[["Z"]][[m]], info$Z[[m]])), rep(1, M))
})
