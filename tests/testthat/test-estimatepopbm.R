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
    global_opts = list(
      nb_cores = 1L,
      backend = "no_mc",
      verbosity = 0L,
      Q_max = 3L
    ),
    fit_opts = fast_fit_opts,
    verbose = FALSE,
    temp_save_path = NULL
  )

  expect_type(result, "list")
  expect_named(
    result,
    c("partition", "cluster", "elapsed_time", "clustering_history")
  )
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
    global_opts = list(
      nb_cores = 1L,
      backend = "no_mc",
      verbosity = 0L,
      Q_max = 3L
    ),
    fit_opts = fast_fit_opts,
    verbose = FALSE,
    temp_save_path = NULL
  )

  expect_type(result2split, "list")
  expect_named(
    result2split,
    c("partition", "cluster", "elapsed_time", "clustering_history")
  )
  expect_true(all(sapply(result2split$partition, inherits, "fitSimpleSBMPop")))
  expect_length(result2split$cluster, length(netlist))
  expect_true(all(
    colnames(result2split$clustering_history) == names(result2split$cluster)
  ))
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
    global_opts = list(
      nb_cores = 1L,
      backend = "no_mc",
      verbosity = 0L,
      Q_max = 3L
    ),
    fit_opts = fast_fit_opts,
    verbose = FALSE,
    temp_save_path = save_path
  )

  expect_true(file.exists(save_path))
  saved <- readRDS(save_path)
  expect_named(
    saved,
    c("partition", "cluster", "elapsed_time", "clustering_history")
  )
  expect_identical(names(saved$cluster), names(result$cluster))
})

test_that("clusterize_bipartite_networks works with valid inputs", {
  set.seed(0)
  alpha1 <- matrix(c(0.8, 0.1, 0.2, 0.7), byrow = TRUE, nrow = 2)
  alpha2 <- matrix(c(0.8, 0.5, 0.5, 0.2), byrow = TRUE, nrow = 2)
  first_collection <- generate_bipartite_collection(
    nr = 50,
    nc = 25,
    pi = c(0.5, 0.5),
    rho = c(0.5, 0.5),
    alpha = alpha1,
    M = 2
  )
  second_collection <- generate_bipartite_collection(
    nr = 50,
    nc = 25,
    pi = c(0.5, 0.5),
    rho = c(0.5, 0.5),
    alpha = alpha2,
    M = 2
  )
  netlist <- append(first_collection, second_collection)

  result <- clusterize_bipartite_networks(
    netlist = netlist,
    colsbm_model = "iid",
    global_opts = list(
      nb_cores = 1,
      backend = "no_mc",
      Q1_max = 2L,
      Q2_max = 2L
    ),
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
    nr = 50,
    nc = 25,
    pi = c(0.5, 0.5),
    rho = c(0.5, 0.5),
    alpha = alpha1,
    M = 1
  )
  second_collection <- generate_bipartite_collection(
    nr = 50,
    nc = 25,
    pi = c(0.5, 0.5),
    rho = c(0.5, 0.5),
    alpha = alpha2,
    M = 1
  )
  netlist <- append(first_collection, second_collection)

  result2split <- clusterize_bipartite_networks(
    netlist = netlist,
    colsbm_model = "iid",
    global_opts = list(
      nb_cores = 1L,
      backend = "no_mc",
      verbosity = 2L,
      Q1_max = 2L,
      Q2_max = 2L
    ),
    fit_opts = common_fit_opts
  )

  expect_type(result2split, "list")
  expect_named(
    result2split,
    c("partition", "cluster", "elapsed_time", "clustering_history")
  )
  expect_true(all(sapply(result2split$partition, inherits, "bisbmpop")))
  expect_length(result2split$cluster, length(netlist))
  expect_true(all(
    colnames(result2split$clustering_history) == names(result2split$cluster)
  ))
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
    nr = 50,
    nc = 25,
    pi = c(0.5, 0.5),
    rho = c(0.5, 0.5),
    alpha = alpha1,
    M = 2
  )
  second_collection <- generate_bipartite_collection(
    nr = 50,
    nc = 25,
    pi = c(0.5, 0.5),
    rho = c(0.5, 0.5),
    alpha = alpha2,
    M = 2
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
    nr = 50,
    nc = 25,
    pi = c(0.5, 0.5),
    rho = c(0.5, 0.5),
    alpha = alpha1,
    M = 2
  )
  second_collection <- generate_bipartite_collection(
    nr = 50,
    nc = 25,
    pi = c(0.5, 0.5),
    rho = c(0.5, 0.5),
    alpha = alpha2,
    M = 2
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

# --- clusterize_bipartite_networks_graphon -----------------------------------

test_that("graphon rejects non-iid colSBM models", {
  netlist <- list(matrix(0, 5, 5), matrix(0, 5, 5))

  expect_error(
    clusterize_bipartite_networks_graphon(
      netlist = netlist,
      colsbm_model = "pi",
      global_opts = list(nb_cores = 1L, backend = "no_mc", verbosity = 0L),
      fit_opts = list(max_vem_steps = 50L),
      verbose = FALSE,
      temp_save_path = NULL
    ),
    "Only iid model is implemented"
  )
})

test_that("graphon works without a temporary save file", {
  set.seed(0)
  alpha1 <- matrix(c(0.8, 0.1, 0.2, 0.7), byrow = TRUE, nrow = 2)
  alpha2 <- matrix(c(0.8, 0.5, 0.5, 0.2), byrow = TRUE, nrow = 2)
  netlist <- append(
    generate_bipartite_collection(
      nr = 25,
      nc = 15,
      pi = c(0.5, 0.5),
      rho = c(0.5, 0.5),
      alpha = alpha1,
      M = 1
    ),
    generate_bipartite_collection(
      nr = 25,
      nc = 15,
      pi = c(0.5, 0.5),
      rho = c(0.5, 0.5),
      alpha = alpha2,
      M = 1
    )
  )

  # temp_save_path = NULL should simply disable checkpointing (as for the
  # other clustering procedures), not crash the process.
  expect_no_error(
    clusterize_bipartite_networks_graphon(
      netlist = netlist,
      colsbm_model = "iid",
      nb_run = 1L,
      global_opts = list(
        nb_cores = 1L,
        backend = "no_mc",
        verbosity = 0L,
        Q1_max = 2L,
        Q2_max = 2L
      ),
      fit_opts = list(max_vem_steps = 200L),
      verbose = FALSE,
      temp_save_path = NULL
    )
  )
})

test_that("graphon returns a consistent partition with a monotone bicl_history", {
  set.seed(0)
  alpha1 <- matrix(c(0.8, 0.1, 0.2, 0.7), byrow = TRUE, nrow = 2)
  alpha2 <- matrix(c(0.8, 0.5, 0.5, 0.2), byrow = TRUE, nrow = 2)
  netlist <- append(
    generate_bipartite_collection(
      nr = 25,
      nc = 15,
      pi = c(0.5, 0.5),
      rho = c(0.5, 0.5),
      alpha = alpha1,
      M = 1
    ),
    generate_bipartite_collection(
      nr = 25,
      nc = 15,
      pi = c(0.5, 0.5),
      rho = c(0.5, 0.5),
      alpha = alpha2,
      M = 1
    )
  )

  res <- clusterize_bipartite_networks_graphon(
    netlist = netlist,
    colsbm_model = "iid",
    nb_run = 1L,
    net_id = c("n1", "n2"),
    global_opts = list(
      nb_cores = 1L,
      backend = "no_mc",
      verbosity = 0L,
      Q1_max = 2L,
      Q2_max = 2L
    ),
    fit_opts = list(max_vem_steps = 200L),
    verbose = FALSE,
    temp_save_path = tempfile(fileext = ".Rds")
  )

  expect_named(res, c("partition", "cluster", "bicl_history"))
  expect_length(res$cluster, length(netlist))
  expect_named(res$cluster, c("n1", "n2"))
  # one sub-collection per cluster label
  expect_length(res$partition, length(unique(res$cluster)))
  # cluster labels are a contiguous permutation 1..K
  expect_setequal(unique(res$cluster), seq_along(unique(res$cluster)))
  # the last recorded BIC-L is the BIC-L of the returned partition
  expect_equal(
    tail(res$bicl_history, 1),
    compute_bicl_partition(res$partition, verbose = FALSE),
    tolerance = 1e-6
  )
  # accepted merges only ever improve the total BIC-L
  expect_true(all(diff(res$bicl_history) > -1e-6))
})

test_that("graphon exposes fusion_history only with keep_history = TRUE", {
  set.seed(0)
  alpha1 <- matrix(c(0.8, 0.1, 0.2, 0.7), byrow = TRUE, nrow = 2)
  alpha2 <- matrix(c(0.8, 0.5, 0.5, 0.2), byrow = TRUE, nrow = 2)
  netlist <- append(
    generate_bipartite_collection(
      nr = 25,
      nc = 15,
      pi = c(0.5, 0.5),
      rho = c(0.5, 0.5),
      alpha = alpha1,
      M = 1
    ),
    generate_bipartite_collection(
      nr = 25,
      nc = 15,
      pi = c(0.5, 0.5),
      rho = c(0.5, 0.5),
      alpha = alpha2,
      M = 1
    )
  )
  go <- list(
    nb_cores = 1L,
    backend = "no_mc",
    verbosity = 0L,
    Q1_max = 2L,
    Q2_max = 2L
  )
  fit_opts <- list(max_vem_steps = 200L)

  res_true <- clusterize_bipartite_networks_graphon(
    netlist = netlist,
    colsbm_model = "iid",
    nb_run = 1L,
    keep_history = TRUE,
    global_opts = go,
    fit_opts = fit_opts,
    verbose = FALSE,
    temp_save_path = tempfile(fileext = ".Rds")
  )
  expect_named(
    res_true,
    c("partition", "cluster", "bicl_history", "fusion_history")
  )
  # the history starts from the separated partition and grows by one per
  # accepted merge, mirroring bicl_history
  expect_length(res_true$fusion_history, length(res_true$bicl_history))

  res_false <- clusterize_bipartite_networks_graphon(
    netlist = netlist,
    colsbm_model = "iid",
    nb_run = 1L,
    keep_history = FALSE,
    global_opts = go,
    fit_opts = fit_opts,
    verbose = FALSE,
    temp_save_path = tempfile(fileext = ".Rds")
  )
  expect_named(res_false, c("partition", "cluster", "bicl_history"))
})

test_that("graphon does not advertise temporary saves when keep_history = FALSE", {
  set.seed(0)
  alpha1 <- matrix(c(0.8, 0.1, 0.2, 0.7), byrow = TRUE, nrow = 2)
  alpha2 <- matrix(c(0.8, 0.5, 0.5, 0.2), byrow = TRUE, nrow = 2)
  netlist <- append(
    generate_bipartite_collection(
      nr = 25,
      nc = 15,
      pi = c(0.5, 0.5),
      rho = c(0.5, 0.5),
      alpha = alpha1,
      M = 1
    ),
    generate_bipartite_collection(
      nr = 25,
      nc = 15,
      pi = c(0.5, 0.5),
      rho = c(0.5, 0.5),
      alpha = alpha2,
      M = 1
    )
  )

  out <- capture.output(
    clusterize_bipartite_networks_graphon(
      netlist = netlist,
      colsbm_model = "iid",
      nb_run = 1L,
      global_opts = list(
        nb_cores = 1L,
        backend = "no_mc",
        verbosity = 0L,
        Q1_max = 2L,
        Q2_max = 2L
      ),
      fit_opts = list(max_vem_steps = 200L),
      verbose = TRUE,
      keep_history = FALSE,
      temp_save_path = tempfile(fileext = ".Rds")
    ),
    type = "message"
  )
  expect_false(any(grepl("Temporary results will be saved", out)))
})

test_that("graphon rejects a partition_init mixing bisbmpop and fitBipartiteSBMPop", {
  set.seed(0)
  alpha1 <- matrix(c(0.8, 0.1, 0.2, 0.7), byrow = TRUE, nrow = 2)
  netlist <- generate_bipartite_collection(
    nr = 25,
    nc = 15,
    pi = c(0.5, 0.5),
    rho = c(0.5, 0.5),
    alpha = alpha1,
    M = 2
  )
  go <- list(
    nb_cores = 1L,
    backend = "no_mc",
    verbosity = 0L,
    Q1_max = 2L,
    Q2_max = 2L
  )
  fit_opts <- list(max_vem_steps = 200L)
  fit1 <- estimate_colBiSBM(
    netlist = netlist[1],
    colsbm_model = "iid",
    nb_run = 1L,
    global_opts = go,
    fit_opts = fit_opts
  )
  fit2 <- estimate_colBiSBM(
    netlist = netlist[2],
    colsbm_model = "iid",
    nb_run = 1L,
    global_opts = go,
    fit_opts = fit_opts
  )

  # a partition_init must be homogeneous: mixing a bisbmpop with a
  # fitBipartiteSBMPop would make the graphon distances undefined
  expect_error(
    clusterize_bipartite_networks_graphon(
      netlist = netlist,
      colsbm_model = "iid",
      nb_run = 1L,
      partition_init = list(fit1, fit2$best_fit),
      global_opts = go,
      fit_opts = fit_opts,
      verbose = FALSE,
      temp_save_path = NULL
    ),
    "same class"
  )
})

# --- clusterize_bipartite_networks_d_a ---------------------------------------

test_that("d_a returns the last fit when max_nb_steps is exhausted", {
  set.seed(0)
  alpha1 <- matrix(c(0.8, 0.1, 0.2, 0.7), byrow = TRUE, nrow = 2)
  alpha2 <- matrix(c(0.8, 0.5, 0.5, 0.2), byrow = TRUE, nrow = 2)
  netlist <- append(
    generate_bipartite_collection(
      nr = 20,
      nc = 12,
      pi = c(0.5, 0.5),
      rho = c(0.5, 0.5),
      alpha = alpha1,
      M = 1
    ),
    generate_bipartite_collection(
      nr = 20,
      nc = 12,
      pi = c(0.5, 0.5),
      rho = c(0.5, 0.5),
      alpha = alpha2,
      M = 1
    )
  )
  save_dir <- tempfile()
  dir.create(save_dir)

  # at step 1 the descending BIC-L is always an improvement over -Inf, so the
  # loop keeps going until max_nb_steps is reached; the function must then
  # return the last ascending result instead of NULL.
  # VEM non-convergence warnings are expected here: on these tiny synthetic
  # networks some block configurations (e.g. Q = 2,1) saturate max_vem_steps.
  res <- suppressWarnings(
    clusterize_bipartite_networks_d_a(
      netlist = netlist,
      colsbm_model = "iid",
      max_nb_steps = 1L,
      nb_run = 1L,
      global_opts = list(
        nb_cores = 1L,
        backend = "no_mc",
        verbosity = 0L,
        Q1_max = 2L,
        Q2_max = 2L
      ),
      fit_opts = list(max_vem_steps = 200L),
      verbose = FALSE,
      temp_save_path = save_dir
    )
  )

  expect_type(res, "list")
  expect_named(res, c("partition", "cluster", "bicl_history"))
  expect_length(res$cluster, length(netlist))
})

test_that("d_a errors when max_nb_steps is below 1", {
  netlist <- list(matrix(0, 4, 4), matrix(0, 4, 4))

  expect_error(
    clusterize_bipartite_networks_d_a(
      netlist = netlist,
      colsbm_model = "iid",
      max_nb_steps = 0L,
      global_opts = list(nb_cores = 1L, backend = "no_mc", verbosity = 0L),
      fit_opts = list(max_vem_steps = 30L),
      verbose = FALSE,
      temp_save_path = tempdir()
    ),
    "max_nb_steps"
  )
})

# --- clusterize_bipartite_networks_low_penalty_first -------------------------

test_that("low_penalty_first returns partition, cluster and a full collection", {
  set.seed(0)
  alpha1 <- matrix(c(0.8, 0.1, 0.2, 0.7), byrow = TRUE, nrow = 2)
  alpha2 <- matrix(c(0.8, 0.5, 0.5, 0.2), byrow = TRUE, nrow = 2)
  netlist <- append(
    generate_bipartite_collection(
      nr = 25,
      nc = 15,
      pi = c(0.5, 0.5),
      rho = c(0.5, 0.5),
      alpha = alpha1,
      M = 1
    ),
    generate_bipartite_collection(
      nr = 25,
      nc = 15,
      pi = c(0.5, 0.5),
      rho = c(0.5, 0.5),
      alpha = alpha2,
      M = 1
    )
  )

  # full_inference = TRUE forces every fit down to the singleton leaves, so the
  # returned partition is always made of sub-collection fits instead of the
  # initial full-collection fit (deterministic exercise of the sub-fit path)
  res <- colSBM:::clusterize_bipartite_networks_low_penalty_first(
    netlist = netlist,
    colsbm_model = "iid",
    nb_run = 1L,
    full_inference = TRUE,
    global_opts = list(
      nb_cores = 1L,
      backend = "no_mc",
      verbosity = 0L,
      Q1_max = 2L,
      Q2_max = 2L
    ),
    fit_opts = list(max_vem_steps = 200L),
    verbose = FALSE,
    temp_save_path = NULL
  )

  expect_named(
    res,
    c(
      "partition",
      "cluster",
      "elapsed_time",
      "clustering_history",
      "full_collection"
    )
  )
  expect_true(inherits(res$full_collection, "bisbmpop"))
  expect_length(res$cluster, length(netlist))
  expect_true(all(sapply(res$partition, inherits, "fitBipartiteSBMPop")))
  # the initial full-collection fit is performed with a null penalty factor
  expect_equal(res$full_collection$fit_opts$penalty_factor, 0)
  # the low-penalty option must be propagated to every sub-collection fit
  expect_true(all(
    vapply(res$partition, function(f) f$fit_opts$penalty_factor, numeric(1)) ==
      0
  ))
})
