test_that("plot functions fails if called on other objects", {
  expect_error(plot.fitSimpleSBMPop("a"))
  expect_error(plot.fitBipartiteSBMPop("a"))
  expect_error(plot.bmpop("a"))
})

test_that("plot.bmpop and plot.fitSimpleSBMPop plot without error", {
  set.seed(1234)
  Net <- lapply(
    list(.7, .7, .2, .2),
    function(p) {
      A <- matrix(0, 15, 15)
      A[lower.tri(A)][sample(15 * 14 / 2, size = round(p * 15 * 14 / 2))] <- 1
      A <- A + t(A)
    }
  )

  cl <- estimate_colSBM(Net,
    colsbm_model = "delta",
    directed = FALSE,
    distribution = "bernoulli",
    nb_run = 1,
    global_opts = list(
      nb_cores = 2,
      verbosity = 0,
      plot_details = 0
    )
  )

  expect_no_error(plot(cl))
  expect_no_error(plot(cl$best_fit, net_id = 1))
  expect_no_error(plot(cl$best_fit, type = "graphon"))
  expect_no_error(plot(cl$best_fit, type = "meso"))
  expect_no_error(plot(cl$best_fit, type = "meso", mixture = TRUE))
  expect_no_error(plot(cl$best_fit, type = "block", net_id = 1))

  data(foodwebs)
  fw_fit <- estimate_colSBM(netlist = foodwebs[1:3], colsbm_model = "iid")

  expect_no_error(plot(fw_fit))
  expect_no_error(plot(fw_fit, type = "graphon"))
  expect_no_error(plot(fw_fit, type = "meso"))
  expect_no_error(plot(fw_fit, type = "meso", mixture = TRUE))
  expect_no_error(plot(fw_fit, type = "block", net_id = 1))
})

test_that("plot.bisbmpop plot.fitBipartiteSBMPop plot without error", {
  set.seed(1234)
  alpha1 <- matrix(c(0.8, 0.1, 0.2, 0.7), byrow = TRUE, nrow = 2)

  first_collection <- generate_bipartite_collection(
    nr = 50,
    nc = 25,
    pi = c(0.5, 0.5),
    rho = c(0.5, 0.5),
    alpha = alpha1,
    M = 2
  )


  # A collection where joint modelisation makes sense
  cl_bipartite <- estimate_colBiSBM(
    netlist = first_collection,
    colsbm_model = "iid",
    global_opts = list(
      nb_cores = 2,
      verbosity = 0,
      plot_details = 0,
      backend = "no_mc"
    )
  )
  expect_no_error(plot(cl_bipartite))
  expect_no_error(plot(cl_bipartite, criterion = "ICL"))
  expect_no_error(plot(cl_bipartite, criterion = "vbound"))
  expect_error(plot(cl_bipartite, criterion = "not possible"))
  expect_no_error(plot(cl_bipartite$best_fit, net_id = 1, type = "graphon"))
  expect_no_error(plot(cl_bipartite$best_fit, type = "meso"))
  expect_no_error(plot(cl_bipartite$best_fit, type = "meso", values = TRUE))
  expect_no_error(plot(cl_bipartite$best_fit, type = "meso", mixture = TRUE))
  expect_no_error(plot(cl_bipartite$best_fit,
    type = "meso", values = TRUE,
    mixture = TRUE
  ))
  expect_no_error(plot(cl_bipartite$best_fit, type = "block", net_id = 1))
})
