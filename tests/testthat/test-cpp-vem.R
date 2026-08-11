test_that("cpp runs normally and stops on non implemented emission distribution", {
  A <- list(
    matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE),
    matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
  )

  expect_no_error(colsbm_create(A, Q = 2L, distribution = "bernoulli"))
  expect_no_error(colsbm_create(A, Q = 2L, distribution = "poisson"))

  # currently not implemented
  expect_error(colsbm_create(A, Q = 2L, distribution = "gaussian"))

  # typos
  expect_error(colsbm_create(A, Q = 2L, distribution = "bernouilli"))
  expect_error(colsbm_create(A, Q = 2L, distribution = "poiton"))
  expect_error(colsbm_create(A, Q = 2L, distribution = "Poisson"))
})

test_that("cpp VE yields the same results as the R procedure", {
  M <- 4
  true_pi <- c(0.1, 0.3, 0.6)
  nbNodes <- 100
  K <- length(true_pi)
  true_alpha <- matrix(rev(seq(1, 10 * K^2, by = 10)) / 100, nrow = K)
  sbm_sampler <- sbm::sampleSimpleSBM(nbNodes = nbNodes, blockProp = true_pi, connectParam = list(mean = true_alpha), directed = TRUE)
  set.seed(1234)
  sbm_realisations <- lapply(seq(M), sbm_sampler[["rNetwork"]])
  A <- lapply(seq(M), function(m) sbm_realisations[[m]][["networkData"]])

  # A random start of Z before optimization
  fake_Z <- lapply(seq(M), function(m) {
    t(sapply(seq(nbNodes), function(idx) as.integer(sample.int(n = K, size = 1, replace = TRUE, prob = true_pi) == seq(K))))
  })

  # colSBM R logic
  R_fit <- fitSimpleSBMPop$new(A = A, Q = 3L, Z = fake_Z, directed = TRUE, free_mixture = FALSE, free_density = FALSE, init_method = "spectral", fit_opts = default_fit_opts_unipartite())
  R_fit$optimize()

  # C++ logic
  ptr <- colsbm_create(A = A, Q = 3L, tau = fake_Z, directed = TRUE, distribution = "bernoulli", free_mixture = FALSE, free_density = FALSE)
  colsbm_optimize(ptr, max_step = 100L, tol = 1e-9)

  info <- colsbm_info(ptr)

  expect_equal(rowSums(info$tau[[1]]), rep(1, nrow(info$tau[[1]])))
  expect_true(all(info$tau[[1]] >= 0))
  expect_true(is.list(info))
  expect_equal(info$M, 4L)
  expect_equal(info$Q, 3L)
  expect_true(is.finite(info$vbound))

  labelize <- function(tau) {
    sapply(seq_len(nrow(tau)), function(i) which.max(tau[i, ]))
  }

  # True quality checks
  cpp_order <- order(diag(info$alpha), decreasing = TRUE)
  expect_identical(R_fit$compute_vbound(), info$vbound)
  cpp_taus <- lapply(info$tau, function(tau) tau[, cpp_order])
  expect_equal(R_fit$tau, cpp_taus, tolerance = 1e-2)
})
