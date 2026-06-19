common_fit_opts <- list(max_vem_steps = 3000L)

test_that("clusterize_unipartite_networks works with valid inputs", {
  fast_fit_opts <- list(max_vem_steps = 100L)
  set.seed(0)
  alpha1 <- matrix(c(0.8, 0.1, 0.2, 0.7), byrow = TRUE, nrow = 2)
  alpha2 <- matrix(c(0.8, 0.5, 0.5, 0.2), byrow = TRUE, nrow = 2)
  first_collection <- generate_unipartite_collection(
    n = 12,
    pi = c(0.5, 0.5),
    alpha = alpha1,
    M = 2
  )
  second_collection <- generate_unipartite_collection(
    n = 12,
    pi = c(0.5, 0.5),
    alpha = alpha2,
    M = 2
  )
  netlist <- append(first_collection, second_collection)

  result <- clusterize_unipartite_networks(
    netlist = netlist,
    colsbm_model = "iid",
    nb_run = 1L,
    global_opts = list(nb_cores = 1L, backend = "no_mc", verbosity = 0L, Q_max = 3L),
    fit_opts = fast_fit_opts,
    verbose = FALSE,
    temp_save_path = NULL
  )

  expect_type(result, "list")
  expect_named(result, c("partition", "cluster", "elapsed_time", "clustering_history"))
  expect_true(all(sapply(result$partition, inherits, "fitSimpleSBMPop")))
  expect_length(result$cluster, length(netlist))
  expect_true(all(colnames(result$clustering_history) == names(result$cluster)))
})

test_that("clusterize_unipartite_networks can handle splitting apart two networks", {
  fast_fit_opts <- list(max_vem_steps = 100L)
  set.seed(0)
  alpha1 <- matrix(c(0.8, 0.1, 0.2, 0.7), byrow = TRUE, nrow = 2)
  alpha2 <- matrix(c(0.8, 0.5, 0.5, 0.2), byrow = TRUE, nrow = 2)
  first_collection <- generate_unipartite_collection(
    n = 12,
    pi = c(0.5, 0.5),
    alpha = alpha1,
    M = 1
  )
  second_collection <- generate_unipartite_collection(
    n = 12,
    pi = c(0.5, 0.5),
    alpha = alpha2,
    M = 1
  )
  netlist <- append(first_collection, second_collection)

  result2split <- clusterize_unipartite_networks(
    netlist = netlist,
    colsbm_model = "iid",
    nb_run = 1L,
    global_opts = list(nb_cores = 1L, backend = "no_mc", verbosity = 0L, Q_max = 3L),
    fit_opts = fast_fit_opts,
    verbose = FALSE,
    temp_save_path = NULL
  )

  expect_type(result2split, "list")
  expect_named(result2split, c("partition", "cluster", "elapsed_time", "clustering_history"))
  expect_true(all(sapply(result2split$partition, inherits, "fitSimpleSBMPop")))
  expect_length(result2split$cluster, length(netlist))
  expect_true(all(colnames(result2split$clustering_history) == names(result2split$cluster)))
})

test_that("clusterize_unipartite_networks handles invalid colsbm_model", {
  netlist <- list(matrix(0, 10, 10), matrix(0, 10, 10))

  expect_error(
    clusterize_unipartite_networks(
      netlist = netlist,
      colsbm_model = "invalid_model",
      global_opts = list(nb_cores = 1L, backend = "no_mc", verbosity = 0L),
      fit_opts = list(max_vem_steps = 50L),
      verbose = FALSE,
      temp_save_path = NULL
    ),
    "`colsbm_model` must be one of"
  )
})

test_that("clusterize_unipartite_networks handles empty netlist", {
  expect_error(
    clusterize_unipartite_networks(
      netlist = list(),
      colsbm_model = "iid",
      global_opts = list(nb_cores = 1L, backend = "no_mc", verbosity = 0L),
      fit_opts = list(max_vem_steps = 50L),
      verbose = FALSE,
      temp_save_path = NULL
    ),
    "`netlist` must be a list of matrices"
  )
})

test_that("clusterize_unipartite_networks saves temporary results", {
  fast_fit_opts <- list(max_vem_steps = 100L)
  set.seed(1)
  alpha <- matrix(c(0.7, 0.2, 0.2, 0.7), byrow = TRUE, nrow = 2)
  netlist <- generate_unipartite_collection(
    n = 10,
    pi = c(0.5, 0.5),
    alpha = alpha,
    M = 2
  )
  save_path <- tempfile(fileext = ".rds")

  result <- clusterize_unipartite_networks(
    netlist = netlist,
    colsbm_model = "iid",
    nb_run = 1L,
    global_opts = list(nb_cores = 1L, backend = "no_mc", verbosity = 0L, Q_max = 3L),
    fit_opts = fast_fit_opts,
    verbose = FALSE,
    temp_save_path = save_path
  )

  expect_true(file.exists(save_path))
  saved <- readRDS(save_path)
  expect_named(saved, c("partition", "cluster", "elapsed_time", "clustering_history"))
  expect_identical(names(saved$cluster), names(result$cluster))
})

test_that("clusterize_bipartite_networks works with valid inputs", {
  set.seed(0)
  alpha1 <- matrix(c(0.8, 0.1, 0.2, 0.7), byrow = TRUE, nrow = 2)
  alpha2 <- matrix(c(0.8, 0.5, 0.5, 0.2), byrow = TRUE, nrow = 2)
  first_collection <- generate_bipartite_collection(
    nr = 50, nc = 25,
    pi = c(0.5, 0.5), rho = c(0.5, 0.5),
    alpha = alpha1, M = 2
  )
  second_collection <- generate_bipartite_collection(
    nr = 50, nc = 25,
    pi = c(0.5, 0.5), rho = c(0.5, 0.5),
    alpha = alpha2, M = 2
  )
  netlist <- append(first_collection, second_collection)

  result <- clusterize_bipartite_networks(
    netlist = netlist,
    colsbm_model = "iid",
    global_opts = list(nb_cores = 1, backend = "no_mc", Q1_max = 2L, Q2_max = 2L),
    fit_opts = common_fit_opts
  )

  expect_type(result, "list")
  expect_true(all(sapply(result$partition, inherits, "bisbmpop")))
})


test_that("clusterize_bipartite_networks can handle splitting apart two networks", {
  set.seed(0)
  alpha1 <- matrix(c(0.8, 0.1, 0.2, 0.7), byrow = TRUE, nrow = 2)
  alpha2 <- matrix(c(0.8, 0.5, 0.5, 0.2), byrow = TRUE, nrow = 2)
  first_collection <- generate_bipartite_collection(
    nr = 50, nc = 25,
    pi = c(0.5, 0.5), rho = c(0.5, 0.5),
    alpha = alpha1, M = 1
  )
  second_collection <- generate_bipartite_collection(
    nr = 50, nc = 25,
    pi = c(0.5, 0.5), rho = c(0.5, 0.5),
    alpha = alpha2, M = 1
  )
  netlist <- append(first_collection, second_collection)

  result2split <- clusterize_bipartite_networks(
    netlist = netlist,
    colsbm_model = "iid",
    global_opts = list(nb_cores = 1L, backend = "no_mc", verbosity = 2L, Q1_max = 2L, Q2_max = 2L),
    fit_opts = common_fit_opts
  )

  expect_type(result2split, "list")
  expect_named(result2split, c("partition", "cluster", "elapsed_time", "clustering_history"))
  expect_true(all(sapply(result2split$partition, inherits, "bisbmpop")))
  expect_length(result2split$cluster, length(netlist))
  expect_true(all(colnames(result2split$clustering_history) == names(result2split$cluster)))
})

test_that("clusterize_bipartite_networks handles invalid colsbm_model", {
  netlist <- list(matrix(0, 10, 10), matrix(0, 10, 10))

  expect_error(
    clusterize_bipartite_networks(
      netlist = netlist,
      colsbm_model = "invalid_model",
      global_opts = list(nb_cores = 1),
      fit_opts = common_fit_opts
    ),
    "`colsbm_model` must be one of"
  )
})

test_that("clusterize_bipartite_networks handles empty netlist", {
  netlist <- list()

  expect_error(
    clusterize_bipartite_networks(
      netlist = netlist,
      colsbm_model = "iid",
      global_opts = list(nb_cores = 1),
      fit_opts = common_fit_opts
    ),
    "`netlist` must be a list of matrices."
  )
})

test_that("clusterize_bipartite_networks works with different distributions", {
  set.seed(0)
  alpha1 <- matrix(c(0.8, 0.1, 0.2, 0.7), byrow = TRUE, nrow = 2)
  alpha2 <- matrix(c(0.8, 0.5, 0.5, 0.2), byrow = TRUE, nrow = 2)
  first_collection <- generate_bipartite_collection(
    nr = 50, nc = 25,
    pi = c(0.5, 0.5), rho = c(0.5, 0.5),
    alpha = alpha1, M = 2
  )
  second_collection <- generate_bipartite_collection(
    nr = 50, nc = 25,
    pi = c(0.5, 0.5), rho = c(0.5, 0.5),
    alpha = alpha2, M = 2
  )
  netlist <- append(first_collection, second_collection)

  result <- clusterize_bipartite_networks(
    netlist = netlist,
    colsbm_model = "iid",
    distribution = "poisson",
    global_opts = list(nb_cores = 1),
    fit_opts = common_fit_opts
  )

  expect_type(result, "list")
  expect_true(all(sapply(result$partition, inherits, "bisbmpop")))
})

test_that("clusterize_bipartite_networks works with full_inference = TRUE", {
  set.seed(0)
  alpha1 <- matrix(c(0.8, 0.1, 0.2, 0.7), byrow = TRUE, nrow = 2)
  alpha2 <- matrix(c(0.8, 0.5, 0.5, 0.2), byrow = TRUE, nrow = 2)
  first_collection <- generate_bipartite_collection(
    nr = 50, nc = 25,
    pi = c(0.5, 0.5), rho = c(0.5, 0.5),
    alpha = alpha1, M = 2
  )
  second_collection <- generate_bipartite_collection(
    nr = 50, nc = 25,
    pi = c(0.5, 0.5), rho = c(0.5, 0.5),
    alpha = alpha2, M = 2
  )
  netlist <- append(first_collection, second_collection)

  result <- clusterize_bipartite_networks(
    netlist = netlist,
    colsbm_model = "iid",
    full_inference = TRUE,
    global_opts = list(nb_cores = 1),
    fit_opts = common_fit_opts
  )

  expect_type(result, "list")
  expect_true(all(sapply(result$partition, inherits, "bisbmpop")))
})

test_that("clusterize_bipartite_networks validates full_collection_init class", {
  netlist <- list(matrix(0, 5, 5), matrix(0, 5, 5))

  expect_error(
    clusterize_bipartite_networks(
      netlist = netlist,
      colsbm_model = "iid",
      full_collection_init = list(),
      global_opts = list(nb_cores = 1L, backend = "no_mc", verbosity = 0L),
      fit_opts = common_fit_opts,
      verbose = FALSE,
      temp_save_path = NULL
    ),
    "full_collection_init should be a bisbmpop object"
  )
})

test_that("clusterize_bipartite_networks validates partition_init type", {
  netlist <- list(matrix(0, 5, 5), matrix(0, 5, 5))

  expect_error(
    clusterize_bipartite_networks(
      netlist = netlist,
      colsbm_model = "iid",
      partition_init = 1L,
      net_id = c("a", "b"),
      global_opts = list(nb_cores = 1L, backend = "no_mc", verbosity = 0L),
      fit_opts = common_fit_opts,
      verbose = FALSE,
      temp_save_path = NULL
    ),
    "partition_init should be a list of fitBipartite objects"
  )
})
