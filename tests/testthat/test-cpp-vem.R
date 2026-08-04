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
