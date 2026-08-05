test_that("cpp VEM scaffold produces finite diagnostics", {
    A <- list(
        matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE),
        matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
    )

    ptr <- colsbm_create(A, Q = 2L)
    expect_true(inherits(ptr, "externalptr"))

    colsbm_optimize(ptr, max_step = 3L, tol = 1e-10)
    info <- colsbm_info(ptr)

    expect_true(is.list(info))
    expect_equal(info$M, 2L)
    expect_equal(info$Q, 2L)
    expect_true(is.finite(info$vbound))
    expect_equal(info$iterations, 3L)
})

test_that("cpp VE step keeps tau rows normalized", {
    A <- list(
        matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
    )

    ptr <- colsbm_create(A, Q = 2L)
    colsbm_optimize(ptr, max_step = 1L, tol = 1e-10)

    info <- colsbm_info(ptr)
    expect_equal(info$M, 1L)
    expect_equal(info$Q, 2L)
    expect_true(is.finite(info$vbound))
    expect_equal(rowSums(info$tau[[1]]), rep(1, nrow(info$tau[[1]])))
    expect_true(all(info$tau[[1]] >= 0))
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
    R_fit$Z

    # C++ logic
    ptr <- colsbm_create(A, Q = 3L)
    colsbm_optimize(ptr, max_step = 100L, tol = 1e-9)

    info <- colsbm_info(ptr)
})
