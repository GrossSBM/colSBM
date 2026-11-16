test_that("Z and MAP$Z are relabelled consistently after reorder_parameters", {
  set.seed(1234)
  alpha <- matrix(c(
    0.9, 0.05, 0.3,
    0.05, 0.8, 0.1,
    0.4, 0.1, 0.6
  ), 3, 3, byrow = TRUE)
  coll <- generate_bipartite_collection(
    c(60, 70, 80), c(50, 60, 70),
    c(0.2, 0.3, 0.5), c(0.5, 0.3, 0.2),
    alpha, 3,
    return_memberships = TRUE
  )
  A <- lapply(coll, function(x) x$incidence_matrix)

  fit <- fitBipartiteSBMPop$new(
    A = A, Q = c(3, 3), distribution = "bernoulli",
    fit_opts = list(verbosity = 0)
  )
  fit$optimize()

  # Blocks are ordered by decreasing marginal connectivity
  expect_false(is.unsorted(rev(rowMeans(fit$MAP$alpha))))
  expect_false(is.unsorted(rev(colMeans(fit$MAP$alpha))))

  for (m in seq_len(fit$M)) {
    for (d in 1:2) {
      tau_map <- apply(fit$tau[[m]][[d]], 1, which.max)
      # Labels must match the reordered tau
      expect_equal(unname(fit$Z[[m]][[d]]), unname(tau_map))
      expect_equal(unname(fit$MAP$Z[[m]][[d]]), unname(tau_map))
      # Node names are kept
      expect_equal(names(fit$Z[[m]][[d]]), names(tau_map))
      # Labels must match the reordered block proportions
      expect_equal(
        as.vector(tabulate(fit$MAP$Z[[m]][[d]], fit$Q[d]) / fit$n[[d]][m]),
        as.vector(fit$MAP$pim[[m]][[d]])
      )
    }
  }
})
