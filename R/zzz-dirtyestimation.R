#' Estimate a colSBM on a collection of networks using C++ and without clever
#' model selection
#'
#' @param netlist A list of matrices.
#' @param colsbm_model Which colSBM to use, one of "iid", "pi", "delta",
#' "deltapi".
#' @param net_id A vector of string, the name of the networks.
#' @param directed A boolean, are the networks directed or not.
#' @param distribution A string, the emission distribution, either "bernoulli"
#' (the default) or "poisson"
#' @param fit_sbm A list of fitted model using the \code{sbm} package.
#' Use to speed up the initialization.
#' @param nb_run An integer, the number of run the algorithm do.
#' @param global_opts Global options for the outer algorithm and the output
#' @param fit_opts Fit options for the VEM algorithm
#' @param Z_init An optional list of cluster memberships for q in 1:Q_max.
#' Default to NULL.
#' @param fit_init Do not use!
#' Optional fit init from where initializing the algorithm.
#'
#'
#' @details The list of parameters \code{global_opts} essentially tunes the
#' exploration process.
#'  \itemize{
#'  \item \code{nb_cores} integer for number of cores used for
#'  parallelization. Default is 1.
#'  \item \code{sbm_init} boolean specifying wether `{sbm}` package
#' should be used to initialize the algorithm.
#'  \item \code{spectral_init} boolean specifying wether a spectral
#' clustering should be used to initialize the algorithm.
#'  \item \code{nb_init} the number of different inits to perform
#' to start the algorithm fit. Default to 10.
#'  \item \code{verbosity} integer for verbosity (0, 1, 2, 3, 4). Default is 1.
#'   0 will disable completely the output of the function.
#'   modelisation is preferred
#'  \item \code{Q_max} integer for the max size to explore. Default
#' is computed with the following formula:
#' `floor(log(sum(sapply(netlist, function(A) nrow(A)))) + 2)`
#'  \item \code{nb_models} the number of models to keep for each values of Q.
#'  Default is 5.
#'  \item \code{depth} specifies how far around the best model the exploration will fit models. Default is 3.
#'  \item \code{plot_details} integer to control the display of the exploration process. Values are 0 or 1. Default is 1.
#'  \item \code{max_pass} the maximum number of passes that will be
#'  executed. Default is 10.
#'  \item \code{backend} the parallelization backend to use.
#' Options available are "no_mc", "future", "parallel". Default is
#' "future".
#' Note: we plan to unify everything behind future.
#' }
#'
#' @return A bmpop object listing a collection of fitted models for the
#' collection of networks
#' @export
#'
#' @seealso [colSBM::clusterize_unipartite_networks()], \code{\link[colSBM]{bmpop}},
#' \code{\link[colSBM]{fitSimpleSBMPop}}, `browseVignettes("colSBM")`
#' @import sbm
#' @examples
#' # Trivial example with Gnp networks:
#' Net <- lapply(
#'   list(.7, .7, .2, .2),
#'   function(p) {
#'     A <- matrix(0, 15, 15)
#'     A[lower.tri(A)][sample(15 * 14 / 2, size = round(p * 15 * 14 / 2))] <- 1
#'     A <- A + t(A)
#'   }
#' )
#' \dontrun{
#' cl <- estimate_colSBM(Net,
#'   colsbm_model = "delta",
#'   directed = FALSE,
#'   distribution = "bernoulli",
#'   nb_run = 1
#' )
#' }
cpp_estimate_colSBM <-
  function(netlist,
           colsbm_model,
           net_id = NULL,
           directed = NULL,
           distribution = "bernoulli",
           fit_sbm = NULL,
           nb_run = 3L,
           global_opts = list(),
           fit_opts = list(),
           Z_init = NULL,
           fit_init = NULL,
           verbose = FALSE) {
    netlist <- check_networks_list(
      networks_list = netlist,
      min_length = 1L
    )

    mixture_and_density <- check_unipartite_colsbm_models(colsbm_model)
    free_mixture <- mixture_and_density[["free_mixture"]]
    free_density <- mixture_and_density[["free_density"]]

    go <- default_global_opts_unipartite(netlist)
    global_opts <- utils::modifyList(x = go, val = global_opts)
    fo <- default_fit_opts_unipartite()
    fit_opts <- utils::modifyList(x = fo, val = fit_opts)

    nb_cores <- global_opts$nb_cores
    if (is.null(global_opts$Q_max)) {
      Q_max <- floor(log(sum(sapply(netlist, nrow))) + 2)
    } else {
      Q_max <- global_opts$Q_max
    }

    if (is.null(global_opts$Q_min)) {
      Q_min <- 2L
    } else {
      Q_min <- global_opts$Q_min
    }

    # Initialisation phase
    # If fit_init is provided, use it as the initial model
    # Otherwise, create the initial model with spectral clustering
    if (!is.null(fit_init)) {
      if (verbose) {
        say("Using provided init for the model at Q=2")
      }
      start_model <- fit_init
    } else {
      if (verbose) {
        say("Initialising spectral clustering")
      }
      spectral_inits <- lapply(netlist, spectral_clustering, K = 2) |> futurize::futurize(seed = TRUE)
      spectral_taus <- lapply(spectral_inits, .one_hot, Q = 2)
      spectral_taus <- lapply(spectral_taus, function(mat) {
        mat[mat < 1e-6] <- 1e-6
        mat[mat > 1 - 1e-6] <- 1 - 1e-6
        mat / rowSums(mat)
      })
      if (verbose) {
        say("Fitting the model at Q=2")
      }
      start_model <- fitSimpleSBMPop$new(
        A = netlist,
        Q = 2L,
        Z = spectral_inits,
        net_id = names(netlist),
        distribution = distribution,
        free_mixture = free_mixture,
        free_density = free_density,
        init_method = "spectral",
        Cpi = NULL,
        Calpha = NULL,
        logfactA = NULL,
        fit_opts = fit_opts
      )

      start_model$optimize()
      if (verbose) {
        say("Model fitted !")
      }
    }

    # Compute Q=1
    model_1 <- fitSimpleSBMPop$new(
      A = netlist,
      Q = 1L,
      Z = lapply(netlist, function(mat) matrix(1, nrow(mat), 1)),
      net_id = names(netlist),
      distribution = distribution,
      free_mixture = free_mixture,
      free_density = free_density,
      init_method = "spectral",
      Cpi = NULL,
      Calpha = NULL,
      logfactA = NULL,
      fit_opts = fit_opts
    )

    model_1$optimize()

    model_list <- list(model_1, start_model)

    best_bicl <- max(sapply(model_list, function(mod) mod$BICL))

    for (iter in seq(10)) {
      # Forward phase: estimate models for Q from 2 to Q_max
      # Pass the initial model to cpp_forward which will compute the splits
      splitted_model_list <- cpp_forward(start_model, netlist, Q_start = 2L, Q_max, distribution, free_mixture, free_density, fit_opts, verbose = verbose, current_model_list = model_list)

      model_list <- compare_model_lists(
        model_list1 = model_list,
        model_list2 = splitted_model_list,
        verbose = verbose
      )

      if (max(sapply(model_list, function(mod) mod$BICL)) > best_bicl) {
        if (verbose) {
          say("BICL criterion increased in forward phase, going for backward")
        }
        best_bicl <- max(sapply(model_list, function(mod) mod$BICL))
      } else {
        if (verbose) {
          say("BICL criterion did not increase in forward phase, stopping optimization.")
        }
        break
      }


      merged_model_list <- cpp_backward(model_list = model_list, netlist = netlist, Q_min = Q_min, distribution = distribution, free_mixture = free_mixture, free_density = free_density, fit_opts = fit_opts, verbose = verbose)

      model_list <- compare_model_lists(
        model_list1 = model_list,
        model_list2 = merged_model_list,
        verbose = verbose
      )
      if (max(sapply(model_list, function(mod) mod$BICL)) > best_bicl) {
        if (verbose) {
          say("BICL criterion increased in backward phase, going for forward")
        }
        best_bicl <- max(sapply(model_list, function(mod) mod$BICL))
      } else {
        if (verbose) {
          say("BICL criterion did not increase in backward phase, stopping optimization.")
        }
        break
      }
    }

    mybmpop <- bmpop$new(
      netlist = netlist,
      net_id = net_id,
      directed = directed,
      distribution = distribution,
      free_density = free_density,
      free_mixture = free_mixture,
      fit_sbm = fit_sbm,
      global_opts = global_opts,
      Z_init = Z_init,
      fit_opts = fit_opts
    )

    mybmpop[["model_list"]] <- model_list
    mybmpop[["BICL"]] <- vapply(model_list, function(mod) mod[["BICL"]], numeric(1))
    mybmpop[["vbound"]] <- vapply(model_list, function(mod) tail(mod[["vbound"]], 1), numeric(1))
    mybmpop[["ICL"]] <- vapply(model_list, function(mod) tail(mod[["ICL"]], 1), numeric(1))
    mybmpop[["best_fit"]] <- model_list[[which.max(mybmpop[["BICL"]])]]

    return(mybmpop)
  }

#' Forward phase estimation for colSBM
#'
#' @param start_model The starting model for the forward phase
#' @param netlist A list of matrices.
#' @param Q_max Maximum number of clusters to consider
#' @param distribution A string, the emission distribution, either "bernoulli" or "poisson"
#' @param free_mixture Boolean indicating if mixture is free
#' @param free_density Boolean indicating if density is free
#' @param fit_opts Fit options for the VEM algorithm
#' @param fit_init Initial fitted model to start the forward phase
#'
#' @return A list of fitted models for Q from 2 to Q_max
cpp_forward <-
  function(start_model, netlist, Q_start, Q_max, distribution, free_mixture, free_density, fit_opts, verbose = FALSE, iteration = NULL, current_model_list) {
    # Initialize taus
    M <- length(netlist)

    if (verbose) {
      say("Computing splits for Q_start=", Q_start, level = 1L)
    }

    # Compute splits from the initial model
    next_splits <- purrr::transpose(lapply(seq(M), function(m) split_clust(X = netlist[[m]], Z = start_model[["Z"]][[m]], Q = Q_start)) |> futurize::futurize(seed = TRUE))
    model_list <- list(start_model)
    for (Q in seq(3, Q_max)) {
      if (verbose) {
        say("Fitting models at Q=", Q, level = 1L)
      }
      models <- lapply(seq_along(next_splits), function(split_idx) {
        tmp_fit <- fitSimpleSBMPop$new(
          A = netlist,
          Q = Q,
          Z = next_splits[[split_idx]],
          net_id = names(netlist),
          distribution = distribution,
          free_mixture = free_mixture,
          free_density = free_density,
          init_method = "spectral",
          Cpi = NULL,
          Calpha = NULL,
          logfactA = NULL,
          fit_opts = fit_opts
        )
        tmp_fit$optimize()
        tmp_fit
      }) |> futurize::futurize(seed = TRUE)
      best_model_idx <- which.max(sapply(models, function(model) model[["BICL"]]))
      if (verbose) {
        bicl <- models[[best_model_idx]][["BICL"]]
        say(
          "Best splitted model index is ", best_model_idx, " with BICL = ", bicl,
          level = 1L
        )
      }
      model_list <- append(model_list, models[[best_model_idx]])
      if (verbose) {
        say("Computing next splits", level = 1L)
      }
      next_splits <- purrr::transpose(lapply(seq(M), function(m) split_clust(X = netlist[[m]], Z = model_list[[Q - 1]][["Z"]][[m]], Q = Q)))
    }

    return(model_list)
  }

#' Backward phase estimation for colSBM
#'
#' @param model_list A list of fitted models from the forward phase
#' @param netlist A list of matrices.
#' @param Q_min Minimum number of clusters to consider
#' @param distribution A string, the emission distribution, either "bernoulli" or "poisson"
#' @param free_mixture Boolean indicating if mixture is free
#' @param free_density Boolean indicating if density is free
#' @param fit_opts Fit options for the VEM algorithm
#' @param verbose Logical, should we print messages?
#'
#' @return A list of fitted models for Q from Q_max down to Q_min
cpp_backward <-
  function(model_list, netlist, Q_min, distribution, free_mixture, free_density, fit_opts, verbose = FALSE) {
    M <- length(netlist)
    # Get the maximum number of clusters from the model list
    Q_max <- length(model_list)

    if (verbose) {
      say(
        "Starting backward phase from Q=", Q_max, " down to Q=", Q_min,
        level = 1L
      )
    }

    # Initialize the result list with the models from the forward phase
    backward_model_list <- model_list

    # Process from Q_max down to Q_min + 1
    for (Q in seq(Q_max, Q_min + 1, by = -1)) {
      if (verbose) {
        say("Fitting merges for Q=", Q - 1, level = 1L)
      }

      # Compute merges using merge_clust function
      Z_merges <- purrr::transpose(lapply(seq(M), function(m) merge_clust(Z = model_list[[Q - 1]][["Z"]][[m]], Q = Q)) |> futurize::futurize(seed = TRUE))

      # Fit models for each merge
      merged_models <- lapply(seq_along(Z_merges), function(merge_idx) {
        # Create a new model with merged clusters
        tmp_fit <- fitSimpleSBMPop$new(
          A = netlist,
          Q = Q - 1,
          Z = Z_merges[[merge_idx]],
          net_id = names(netlist),
          distribution = distribution,
          free_mixture = free_mixture,
          free_density = free_density,
          init_method = "spectral",
          Cpi = NULL,
          Calpha = NULL,
          logfactA = NULL,
          fit_opts = fit_opts
        )

        # Optimize the model
        tmp_fit$optimize()
        tmp_fit
      }) |> futurize::futurize(seed = TRUE)

      # Select the best model among the merged ones
      best_model_idx <- which.max(sapply(merged_models, function(model) model[["BICL"]]))

      if (verbose) {
        bicl <- merged_models[[best_model_idx]][["BICL"]]
        say(
          "Best merged model index is ", best_model_idx, " with BICL = ", bicl,
          level = 1L
        )
      }

      # Add the best merged model to the backward list
      backward_model_list[[Q - 1]] <- merged_models[[best_model_idx]]
    }

    return(backward_model_list)
  }

#' Compare two list of models to merge them by the best BICL
#'
#' @param model_list1 a list of fit...SBM with a BICL attribute
#' @param model_list2 a second list of fit...SBM with a BICL
#' attribute
#' @param verbose a boolean that controls if the function print
#' messages. Default to TRUE
#'
#' @return a list of fit...SBM merged by selecting the best BICL
#' for each model size Q provided
compare_model_lists <- function(model_list1, model_list2, verbose = TRUE) {
  # --- Normalise each list: keep the best model per Q ---------------------
  # A list may hold several models fitted at the same Q (e.g. from different
  # splits/merges). Only the best BICL per Q matters for the comparison, so
  # we collapse each list to one model per Q first.
  best_per_Q <- function(model_list) {
    if (length(model_list) == 0L) {
      return(list())
    }
    Qs <- vapply(model_list, function(mod) mod[["Q"]], numeric(1))
    BICLs <- vapply(model_list, function(mod) mod[["BICL"]], numeric(1))
    # For each Q, keep the index of the model with the highest BICL.
    # In case of a tie, the first occurrence wins (stable order).
    best_idx <- tapply(seq_along(Qs), Qs, function(idx) idx[which.max(BICLs[idx])])
    model_list[best_idx]
  }

  list1 <- best_per_Q(model_list1)
  list2 <- best_per_Q(model_list2)

  Q1 <- vapply(list1, function(mod) mod[["Q"]], numeric(1))
  Q2 <- vapply(list2, function(mod) mod[["Q"]], numeric(1))

  # --- Align by Q value (robust to list order) ----------------------------
  common_Q <- intersect(Q1, Q2)
  only1_Q <- setdiff(Q1, Q2)
  only2_Q <- setdiff(Q2, Q1)

  # --- Compare BICL on the common Qs --------------------------------------
  merged <- vector("list", length(common_Q))
  for (i in seq_along(common_Q)) {
    k <- common_Q[i]
    mod1 <- list1[[match(k, Q1)]]
    mod2 <- list2[[match(k, Q2)]]
    # In case of a tie, keep the first list (mod1)
    if (mod1[["BICL"]] >= mod2[["BICL"]]) {
      merged[[i]] <- mod1
    } else {
      merged[[i]] <- mod2
    }
    if (verbose) {
      print_bicl <- merged[[i]][["BICL"]]
      say(
        "Best model for Q=", k, " has BICL = ", print_bicl,
        level = 1L
      )
    }
  }

  # --- Re-add the non-aligned models --------------------------------------
  if (length(only1_Q) > 0L) {
    merged <- c(merged, list1[match(only1_Q, Q1)])
  }
  if (length(only2_Q) > 0L) {
    merged <- c(merged, list2[match(only2_Q, Q2)])
  }

  # --- Sort by Q ascending ------------------------------------------------
  merged_Q <- vapply(merged, function(mod) mod[["Q"]], numeric(1))
  merged[order(merged_Q)]
}

#' Partition of a collection of unipartite networks based on their common
#' mesoscale structures
#'
#' @param netlist A list of matrices.
#' @param colsbm_model Which colSBM to use, one of "iid", "pi", "delta", "deltapi",
#' @param directed A boolean, should the networks be considered as directed or not.
#' @param net_id A vector of string, the name of the networks.
#' @param distribution A string, the emission distribution, either "bernoulli"
#' (the default) or "poisson"
#' @param nb_run An integer, the number of run the algorithm do.
#' @param global_opts Global options for the outer algorithm and
#' the output. See `estimate_colSBM` for more informations on the elements of the list.
#' @param fit_opts Fit options for the VEM algorithm
#' @param fit_init Do not use!
#' Optional fit init from where initializing the algorithm.
#' @param full_inference The default "FALSE", the algorithm stop once splitting
#' groups of networks does not improve the BICL criterion. If "TRUE", then
#' continue to split groups until a trivial classification of one network per
#' group.
#' @param verbose A boolean, should the function be verbose or not. Default to
#' TRUE.
#'
#' @param temp_save_path A string, the path where to save the temporary results.
#' Defaults to a temporary file.
#'
#' @importFrom future.apply future_lapply
#' @import cli
#' @importFrom utils modifyList
#'
#' @return A list with four elements:
#' \item{partition}{A list of models giving the best partition.}
#' \item{cluster}{A vector of integers giving the cluster of each network.}
#' \item{elapsed_time}{The total time taken by the clustering procedure.}
#' \item{clustering_history}{A matrix with M columns and has much rows as there are cuts during partitioning.}
#'
#' @details
#' This function makes call to `estimate_colSBM`.
#' @export
#'
#' @seealso [colSBM::estimate_colSBM()],
#' \code{\link[colSBM]{fitSimpleSBMPop}}, `browseVignettes("colSBM")`
#'
#' @examples
#' alpha1 <- matrix(c(0.8, 0.1, 0.2, 0.7), byrow = TRUE, nrow = 2)
#' alpha2 <- matrix(c(0.8, 0.5, 0.5, 0.2), byrow = TRUE, nrow = 2)
#' first_collection <- generate_unipartite_collection(
#'   n = 50,
#'   pi = c(0.5, 0.5),
#'   alpha = alpha1, M = 2
#' )
#' second_collection <- generate_unipartite_collection(
#'   n = 50,
#'   pi = c(0.5, 0.5),
#'   alpha = alpha2, M = 2
#' )
#'
#' netlist <- append(first_collection, second_collection)
#'
#' \dontrun{
#' cl_separated <- clusterize_networks(
#'   netlist = netlist,
#'   colsbm_model = "iid"
#' )
#' }
clusterize_unipartite_networks_cpp <- function(netlist,
                                               colsbm_model,
                                               directed = FALSE,
                                               net_id = NULL,
                                               distribution = "bernoulli",
                                               nb_run = 3L,
                                               global_opts = list(),
                                               fit_opts = list(),
                                               fit_init = NULL,
                                               full_inference = FALSE,
                                               verbose = TRUE,
                                               temp_save_path = tempfile(fileext = ".Rds")) {
  check_unipartite_colsbm_models(colsbm_model = colsbm_model)
  # Check if a netlist is provided, try to cast it if not
  netlist <- check_networks_list(networks_list = netlist)
  net_id <- check_net_id_and_initialize(net_id = net_id, networks_list = netlist)
  check_colsbm_emission_distribution(emission_distribution = distribution)
  check_networks_list_match_emission_distribution(
    networks_list = netlist,
    emission_distribution = distribution
  )
  go <- default_global_opts_unipartite(netlist = netlist)
  go$plot_details <- 0
  go <- utils::modifyList(go, global_opts)
  global_opts <- go
  fo <- default_fit_opts_unipartite()
  fo <- utils::modifyList(fo, fit_opts)
  fit_opts <- fo

  # Fit the initial model on the full collection
  if (verbose) {
    if (!is.null(temp_save_path)) {
      cli::cli_alert_info("A save file will be created at {.val {temp_save_path}} and updated after each step")
    }
    cli::cli_h1("Fitting the full collection")
  }
  start_time <- Sys.time()
  my_sbmpop <- cpp_estimate_colSBM(
    netlist = netlist,
    colsbm_model = colsbm_model,
    directed = directed,
    net_id = net_id,
    distribution = distribution,
    nb_run = nb_run,
    global_opts = global_opts,
    fit_opts = fit_opts
  )

  clustering_queue <- list(my_sbmpop)
  list_model_binary <- list()
  cluster <- rep(1, length(netlist))
  names(cluster) <- net_id

  clustering_history <- as.data.frame(matrix(cluster, nrow = 1L))


  if (verbose) {
    cli::cli_h1("Beginning clustering")
  }
  # Process the clustering queue
  while (length(clustering_queue) > 0) {
    if (!is.null(temp_save_path)) {
      saveRDS(list(
        clustering_queue = clustering_queue,
        list_model_binary = list_model_binary,
        clustering_history = clustering_history
      ), temp_save_path)
    }
    fit <- clustering_queue[[1]]
    clustering_queue <- clustering_queue[-1]

    # If the collection contains only one network, add it to the final list
    if (inherits(fit, "bmpop") && fit$best_fit$M == 1) {
      list_model_binary <- append(list_model_binary, list(fit$best_fit))
      next
    }
    if (inherits(fit, "bisbmpop") && fit$best_fit$M == 1) {
      list_model_binary <- append(list_model_binary, list(fit$best_fit))
      next
    }
    if (inherits(fit, "fitBipartiteSBMPop") && fit$M == 1) {
      list_model_binary <- append(list_model_binary, list(fit))
      next
    }

    # Compute the dissimilarity matrix
    dist_bm <- compute_dissimilarity_matrix(collection = fit)
    # Partition the networks based on the dissimilarity matrix
    cl <- partition_networks_list_from_dissimilarity(
      networks_list = fit$A,
      dissimilarity_matrix = dist_bm,
      nb_groups = 2L
    )

    if (verbose) {
      cli::cli_h2("Trying to split the collection of {.val {fit$net_id}}")
    }
    # Fit models for the sub-collections
    fits <- colsbm_lapply(
      c(1, 2),
      function(k) {
        Z_init <- lapply(
          fit$model_list,
          function(mod) mod$Z[cl == k]
        )

        if (verbose) {
          cli::cli_alert_info("Fitting a sub collection with : {.val {fit$net_id[cl == k]}}")
        }

        return(
          cpp_estimate_colSBM(
            netlist = fit$A[cl == k],
            colsbm_model = colsbm_model,
            net_id = fit$net_id[cl == k],
            distribution = distribution,
            nb_run = min(sum(cl == k), nb_run),
            Z_init = Z_init,
            global_opts = global_opts,
            fit_opts = fit_opts,
            fit_sbm = fit$fit_sbm[cl == k],
          )
        )
      },
      backend = global_opts[["backend"]],
      nb_cores = global_opts[["nb_cores"]]
    )


    bicl_increased <- (fits[[1]]$best_fit$BICL + fits[[2]]$best_fit$BICL > fit$best_fit$BICL)
    # Decide whether to continue splitting or add to final list
    if (full_inference || bicl_increased) {
      clustering_queue <- append(clustering_queue, fits)
      if (verbose && bicl_increased) {
        cli::cli_alert_success("Splitting collections improved the BIC-L criterion")
      }
      if (verbose && full_inference) {
        cli::cli_alert_info("Full inference mode enabled, continuing to split collections")
      }
      prev_cluster <- unique(cluster[fit$net_id])

      #  Making room for a new cluster
      cluster[cluster > prev_cluster] <- cluster[cluster > prev_cluster] + 1

      # Assign the new cluster to the networks
      cluster[fit$net_id[cl == 2]] <- prev_cluster + 1
      clustering_history <- rbind(clustering_history, matrix(unname(cluster), nrow = 1))
    } else {
      list_model_binary <- append(list_model_binary, list(fit$best_fit))
      if (verbose) {
        cli::cli_alert_danger("Splitting collections {.emph decreased} the BIC-L criterion")
      }
    }
  }

  # Final message indicating the end of clustering
  if (verbose) {
    cli::cli_alert_success("Finished clustering")
  }

  colnames(clustering_history) <- net_id

  output_list <- list(
    partition = list_model_binary,
    cluster = cluster,
    elapsed_time = Sys.time() - start_time,
    clustering_history = clustering_history
  )
  if (!is.null(temp_save_path)) {
    saveRDS(output_list, temp_save_path)
    if (verbose) {
      cli::cli_alert_info("The final results are saved at {.val {temp_save_path}}")
    }
  }
  return(output_list)
}
