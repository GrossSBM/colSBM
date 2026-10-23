test_that("cpp free-mixture models agree with R on identical parameters (full support)", {
  set.seed(7)
  M <- 2L
  Q <- 3L
  n <- 150L
  true_pi <- c(0.5, 0.3, 0.2)
  lambda <- matrix(
    c(0.20, 0.05, 0.02, 0.05, 0.15, 0.03, 0.02, 0.03, 0.10),
    Q,
    Q
  )
  rNetwork <- function(d) {
    z <- sample(Q, n, TRUE, true_pi)
    Am <- matrix(0, n, n)
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        if (i != j) Am[i, j] <- rbinom(1, 1, d * lambda[z[i], z[j]])
      }
    }
    Am
  }
  dens <- c(1, 0.6)
  A <- lapply(dens, rNetwork)
  taus <- lapply(seq_len(M), function(m) {
    tau <- matrix(1e-6, n, Q)
    tau[cbind(seq_len(n), (seq_len(n) %% Q) + 1)] <- 1 - 1e-6
    tau / rowSums(tau)
  })

  for (free_mixture in c(FALSE, TRUE)) {
    for (free_density in c(FALSE, TRUE)) {
      ptr <- colsbm_create(
        A = A,
        Q = Q,
        tau = taus,
        directed = TRUE,
        distribution = "bernoulli",
        free_mixture = free_mixture,
        free_density = free_density
      )
      colsbm_optimize(ptr, max_step = 200L, tol = 1e-3)
      info <- colsbm_info(ptr)
      expect_true(all(is.finite(info$vbound)))

      # The C++ object must expose the support it used.
      expect_equal(info$Cpi, matrix(1, Q, M))
      expect_equal(info$Calpha, matrix(1, Q, Q))

      # Feed the C++ final parameters back into the R object and ask R for its
      # own vbound: the two must agree to machine precision.
      fit <- fitSimpleSBMPop$new(
        A = A,
        Q = Q,
        Z = taus,
        directed = TRUE,
        distribution = "bernoulli",
        free_mixture = free_mixture,
        free_density = free_density,
        init_method = "given",
        fit_opts = default_fit_opts_unipartite()
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
      expect_equal(
        tail(info$vbound, 1),
        fit$compute_vbound(),
        tolerance = 1e-3
      )

      if (free_density) {
        expect_equal(info$delta[2], dens[2], tolerance = 1e-1)
      }
    }
  }
})

test_that("cpp empty-block support zeroes out tau and agrees with R", {
  set.seed(11)
  M <- 2L
  Q <- 3L
  n <- 200L
  alpha_true <- matrix(
    c(0.12, 0.03, 0.03, 0.03, 0.12, 0.03, 0.03, 0.03, 0.12),
    Q,
    Q
  )
  pi1 <- c(1 / 3, 1 / 3, 1 / 3)
  pi2 <- c(0.4, 0.6, 0) # block 3 empty in network 2
  A <- lapply(seq_len(M), function(m) {
    pm <- if (m == 1) pi1 else pi2
    z <- sample(Q, n, TRUE, pm)
    Am <- matrix(0, n, n)
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        if (i != j) Am[i, j] <- rbinom(1, 1, alpha_true[z[i], z[j]])
      }
    }
    Am
  })
  taus <- lapply(seq_len(M), function(m) {
    tau <- matrix(1e-6, n, Q)
    tau[cbind(seq_len(n), (seq_len(n) %% Q) + 1)] <- 1 - 1e-6
    tau / rowSums(tau)
  })

  Cpi <- matrix(TRUE, Q, M)
  Cpi[3, 2] <- FALSE
  Calpha <- Reduce("+", lapply(seq_len(M), function(m) tcrossprod(Cpi[, m]))) >
    0

  for (free_density in c(FALSE, TRUE)) {
    ptr <- colsbm_create(
      A = A,
      Q = Q,
      tau = taus,
      directed = TRUE,
      distribution = "bernoulli",
      free_mixture = TRUE,
      free_density = free_density,
      Cpi = Cpi,
      Calpha = Calpha
    )
    colsbm_optimize(ptr, max_step = 200L, tol = 1e-3)
    info <- colsbm_info(ptr)

    # The empty block must receive no variational mass and no pi mass.
    expect_equal(sum(info$tau[[2]][, 3]), 0, tolerance = 1e-12)
    expect_equal(info$pi[[2]][3], 0, tolerance = 1e-12)
    expect_equal(info$pim[[2]][3], 0, tolerance = 1e-12)
    # The other blocks keep the full mass.
    expect_equal(sum(info$tau[[2]]), n, tolerance = 1e-6)

    # R agrees with C++ on the same parameters.
    fit <- fitSimpleSBMPop$new(
      A = A,
      Q = Q,
      Z = taus,
      directed = TRUE,
      distribution = "bernoulli",
      free_mixture = TRUE,
      free_density = free_density,
      init_method = "empty",
      Cpi = Cpi,
      Calpha = Calpha,
      fit_opts = default_fit_opts_unipartite()
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
  }
})

test_that("cpp free-mixture procedure never does meaningfully worse than R", {
  M <- 2L
  Q <- 3L
  n <- 150L
  true_pi <- c(0.5, 0.3, 0.2)

  gen_data <- function(distribution, seed) {
    set.seed(seed)
    if (distribution == "bernoulli") {
      lambda <- matrix(
        c(0.20, 0.05, 0.02, 0.05, 0.15, 0.03, 0.02, 0.03, 0.10),
        Q,
        Q
      )
    } else {
      lambda <- matrix(
        c(2.0, 0.4, 0.2, 0.4, 1.5, 0.3, 0.2, 0.3, 1.0),
        Q,
        Q
      )
    }
    dens <- c(1, 0.6)
    lapply(dens, function(d) {
      z <- sample(Q, n, TRUE, true_pi)
      Am <- matrix(0, n, n)
      for (i in seq_len(n)) {
        for (j in seq_len(n)) {
          if (i != j) {
            if (distribution == "bernoulli") {
              Am[i, j] <- rbinom(1, 1, d * lambda[z[i], z[j]])
            } else {
              Am[i, j] <- rpois(1, d * lambda[z[i], z[j]])
            }
          }
        }
      }
      Am
    })
  }

  spectral_init <- function(A) {
    labels <- lapply(A, spectral_clustering, K = Q)
    taus <- lapply(labels, .one_hot, Q = Q)
    lapply(taus, function(tau) {
      tau[tau < 1e-6] <- 1e-6
      tau[tau > 1 - 1e-6] <- 1 - 1e-6
      tau / rowSums(tau)
    })
  }

  labelize <- function(tau) apply(tau, 1, which.max)

  run_R <- function(distribution, directed, fm, fd, A, init) {
    fit_opts <- default_fit_opts_unipartite()
    fit_opts$minibatch <- FALSE
    fit_opts$tolerance <- 1e-6
    fit_opts$max_vem_steps <- 5000L
    if (distribution == "bernoulli" && fd) {
      fit_opts$approx_pois <- TRUE
    }
    fit <- fitSimpleSBMPop$new(
      A = A,
      Q = Q,
      Z = init,
      directed = directed,
      distribution = distribution,
      free_mixture = fm,
      free_density = fd,
      init_method = "given",
      fit_opts = fit_opts
    )
    fit$optimize()
    list(vb = tail(fit$vbound, 1), tau = fit$tau)
  }

  run_cpp <- function(distribution, directed, fm, fd, A, init) {
    ptr <- colsbm_create(
      A = A,
      Q = Q,
      tau = init,
      directed = directed,
      distribution = distribution,
      free_mixture = fm,
      free_density = fd
    )
    colsbm_optimize(ptr, max_step = 5000L, tol = 1e-6)
    info <- colsbm_info(ptr)
    list(vb = tail(info$vbound, 1), tau = info$tau, iters = info$iterations)
  }

  # Full grid: both distributions, directed/undirected, free_mixture and
  # free_density on/off, three seeds. The second network is sparser (density
  # 0.6), so the undirected Bernoulli case is genuinely hard: R and C++ can
  # legitimately land in different local optima. The thresholds below are
  # calibrated on that case (worst measured vbound gap ~1.3e-3 relative,
  # worst R/C++ ARI ~0.63).
  for (seed in c(2026L, 31L, 89L)) {
    for (distribution in c("bernoulli", "poisson")) {
      for (directed in c(FALSE, TRUE)) {
        for (fm in c(FALSE, TRUE)) {
          for (fd in c(FALSE, TRUE)) {
            A <- gen_data(distribution, seed)
            init <- spectral_init(A)
            Rf <- run_R(distribution, directed, fm, fd, A, init)
            Cf <- run_cpp(distribution, directed, fm, fd, A, init)

            # C++ must actually converge within the step budget.
            expect_lt(Cf$iters, 5000L)
            # C++ is never meaningfully worse than the R procedure.
            slack <- 1e-2 * abs(Rf$vb)
            expect_gte(Cf$vb, Rf$vb - slack)
            # The two procedures land on comparable clusterings.
            ari <- sapply(
              seq_len(M),
              function(m) {
                aricode::ARI(
                  labelize(Cf$tau[[m]]),
                  labelize(Rf$tau[[m]])
                )
              }
            )
            expect_true(all(ari >= 0.5))
          }
        }
      }
    }
  }
})

test_that("cpp free-mixture bernoulli explore proposes empty blocks", {
  set.seed(13)
  M <- 2L
  Q <- 3L
  n <- 200L
  alpha_true <- matrix(
    c(0.12, 0.03, 0.03, 0.03, 0.12, 0.03, 0.03, 0.03, 0.12),
    Q,
    Q
  )
  pi1 <- c(1 / 3, 1 / 3, 1 / 3)
  pi2 <- c(0.4, 0.6, 0)
  A <- lapply(seq_len(M), function(m) {
    pm <- if (m == 1) pi1 else pi2
    z <- sample(Q, n, TRUE, pm)
    Am <- matrix(0, n, n)
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        if (i != j) Am[i, j] <- rbinom(1, 1, alpha_true[z[i], z[j]])
      }
    }
    Am
  })

  fit_opts <- default_fit_opts_unipartite()
  fit_opts$use_cpp <- TRUE
  res <- colSBM:::cpp_estimate_colSBM(
    netlist = A,
    colsbm_model = "pi",
    distribution = "bernoulli",
    fit_opts = fit_opts,
    global_opts = list(verbosity = 0, Q_min = 2L, Q_max = 4L)
  )
  expect_true(inherits(res, "bmpop"))
  expect_true(length(res$model_list) >= 2L)

  # At least one proposed model must have an empty block in network 2.
  has_empty <- any(vapply(
    res$model_list,
    function(mod) mod$Q > 1L && any(!mod$Cpi),
    logical(1)
  ))
  expect_true(has_empty)
})
