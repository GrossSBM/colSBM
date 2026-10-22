test_that("cpp runs normally and stops on non implemented emission distribution", {
  A <- list(
    matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE),
    matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
  )

  taus <- list(
    matrix(rep(c(0, 1), nrow(A[[1]])), nrow = 2, byrow = TRUE),
    matrix(rep(c(1, 0), nrow(A[[1]])), nrow = 2, byrow = TRUE)
  )

  expect_no_error(colsbm_create(
    A,
    Q = 2L,
    tau = taus,
    distribution = "bernoulli"
  ))
  expect_no_error(colsbm_create(
    A,
    Q = 2L,
    tau = taus,
    distribution = "poisson"
  ))

  # currently not implemented
  expect_error(colsbm_create(A, Q = 2L, tau = taus, distribution = "gaussian"))

  # typos
  expect_error(colsbm_create(
    A,
    Q = 2L,
    tau = taus,
    distribution = "bernouilli"
  ))
  expect_error(colsbm_create(A, Q = 2L, tau = taus, distribution = "poiton"))
  expect_error(colsbm_create(A, Q = 2L, tau = taus, distribution = "Poisson"))
})

test_that("cpp poisson vbound agrees with the R procedure on identical parameters", {
  set.seed(42)
  M <- 4L
  Q <- 3L
  n <- 120L
  true_pi <- c(0.5, 0.3, 0.2)
  dens <- c(1, 0.5, 1 / 3, 1)
  lambda <- matrix(c(4, 1, 0.5, 1, 3, 0.5, 0.5, 0.5, 2), Q, Q)
  rNetwork <- function(d) {
    z <- sample(Q, n, TRUE, true_pi)
    Am <- matrix(0, n, n)
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        if (i != j) Am[i, j] <- rpois(1, d * lambda[z[i], z[j]])
      }
    }
    Am
  }
  A <- lapply(dens, rNetwork)
  taus <- lapply(seq_len(M), function(m) {
    tau <- matrix(1e-6, n, Q)
    tau[cbind(seq_len(n), (seq_len(n) %% Q) + 1)] <- 1 - 1e-6
    tau / rowSums(tau)
  })

  for (directed in c(TRUE, FALSE)) {
    for (free_density in c(FALSE, TRUE)) {
      ptr <- colsbm_create(
        A = A,
        Q = Q,
        tau = taus,
        directed = directed,
        distribution = "poisson",
        free_mixture = TRUE,
        free_density = free_density
      )
      colsbm_optimize(ptr, max_step = 200L, tol = 1e-3)
      info <- colsbm_info(ptr)
      expect_true(all(is.finite(info$vbound)))
      expect_equal(info$M, M)
      expect_equal(info$Q, Q)

      # Feed the C++ final parameters back into the R object and ask R for its
      # own vbound: the two must agree to machine precision.
      fit <- colSBM:::fitSimpleSBMPop$new(
        A = A,
        Q = Q,
        Z = taus,
        directed = directed,
        distribution = "poisson",
        free_mixture = TRUE,
        free_density = free_density,
        init_method = "given",
        fit_opts = colSBM:::default_fit_opts_unipartite()
      )
      fit$tau <- info$tau
      fit$pi <- info$pi
      fit$pim <- info$pim
      fit$alpha <- info$alpha
      fit$delta <- info$delta
      fit$alpham <- info$alpham
      fit$emqr <- aperm(info$emqr, c(3, 1, 2))
      fit$nmqr <- aperm(info$nmqr, c(3, 1, 2))
      fit$logfactA <- vapply(
        seq_along(A),
        function(m) sum(lfactorial(A[[m]]) * (1 - fit$mask[[m]]), na.rm = TRUE),
        FUN.VALUE = .1
      )
      expect_equal(tail(info$vbound, 1), fit$compute_vbound(), tolerance = 1e-3)

      # The delta models must recover the network densities
      if (free_density) {
        expect_equal(info$delta[2], 1 / 2, tolerance = 3e-2)
        expect_equal(info$delta[3], 1 / 3, tolerance = 3e-2)
      }
    }
  }
})

test_that("cpp VE yields the same results as the R procedure on modular networks", {
  M <- 4L
  true_pi <- c(0.6, 0.3, 0.1)
  nbNodes <- 1000
  K <- length(true_pi)
  true_alpha <- diag(c(0.9, 0.5, 0.1))
  sbm_sampler <- sbm::sampleSimpleSBM(
    nbNodes = nbNodes,
    blockProp = true_pi,
    connectParam = list(mean = true_alpha),
    directed = TRUE
  )
  set.seed(1234)
  sbm_realisations <- lapply(seq(M), sbm_sampler[["rNetwork"]])
  effective_pim <- lapply(sbm_realisations, function(sbm_real) {
    sbm_real[["indMemberships"]] |> colMeans()
  })
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
  R_fit <- fitSimpleSBMPop$new(
    A = A,
    Q = 3L,
    Z = spectral_taus,
    directed = TRUE,
    free_mixture = FALSE,
    free_density = FALSE,
    init_method = "given",
    fit_opts = default_fit_opts_unipartite()
  )
  R_fit$optimize()

  # C++ logic
  ptr <- colsbm_create(
    A = A,
    Q = K,
    tau = spectral_taus,
    directed = TRUE,
    distribution = "bernoulli",
    free_mixture = FALSE,
    free_density = FALSE
  )
  colsbm_optimize(ptr, max_step = 100L, tol = 1e-3)

  info <- colsbm_info(ptr)
  labelize <- function(tau) {
    sapply(seq_len(nrow(tau)), function(i) which.max(tau[i, ]))
  }
  info$Z <- lapply(info$tau, labelize)

  # No incorrect output from C++
  expect_equal(
    object = rowSums(info$tau[[1]]),
    expected = rep(1, nrow(info$tau[[1]])),
    tolerance = 1e-6
  )
  expect_true(all(info$tau[[1]] >= 0))
  expect_type(info, type = "list")
  expect_identical(info$M, M)
  expect_identical(info$Q, K)
  expect_all_true(is.finite(info$vbound))

  # True quality checks
  expect_equal(R_fit$tau, info$tau, tolerance = 1e-6)
  ## Parameters checks

  ## Check that clusterings coincides
  expect_identical(
    sapply(seq_along(info$Z), function(m) {
      aricode::ARI(R_fit$Z[[m]], info$Z[[m]])
    }),
    rep(1, M)
  )
})

test_that("cpp VE yields the same results as the R procedure on core-periphery networks", {
  M <- 4
  true_pi <- c(0.6, 0.4)
  nbNodes <- 1000
  K <- length(true_pi)
  true_alpha <- matrix(
    c(
      0.9,
      0.7,
      0.35,
      0.05
    ),
    nrow = K,
    byrow = TRUE
  )
  sbm_sampler <- sbm::sampleSimpleSBM(
    nbNodes = nbNodes,
    blockProp = true_pi,
    connectParam = list(mean = true_alpha),
    directed = TRUE
  )
  set.seed(1234)
  sbm_realisations <- lapply(seq(M), sbm_sampler[["rNetwork"]])
  effective_pim <- lapply(sbm_realisations, function(sbm_real) {
    sbm_real[["indMemberships"]] |> colMeans()
  })
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
  R_fit <- fitSimpleSBMPop$new(
    A = A,
    Q = K,
    Z = spectral_taus,
    directed = TRUE,
    free_mixture = FALSE,
    free_density = FALSE,
    init_method = "given",
    fit_opts = default_fit_opts_unipartite()
  )
  R_fit$optimize()

  # C++ logic
  ptr <- colsbm_create(
    A = A,
    Q = K,
    tau = spectral_taus,
    directed = TRUE,
    distribution = "bernoulli",
    free_mixture = FALSE,
    free_density = FALSE
  )
  colsbm_optimize(ptr, max_step = 100L, tol = 1e-3)

  info <- colsbm_info(ptr)
  labelize <- function(tau) {
    sapply(seq_len(nrow(tau)), function(i) which.max(tau[i, ]))
  }
  info$Z <- lapply(info$tau, labelize)

  # No incorrect output from C++
  expect_equal(
    object = rowSums(info$tau[[1]]),
    expected = rep(1, nrow(info$tau[[1]])),
    tolerance = 1e-6
  )
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
  expect_identical(
    sapply(seq_along(info$Z), function(m) {
      aricode::ARI(R_fit[["Z"]][[m]], info$Z[[m]])
    }),
    rep(1, M)
  )
})
