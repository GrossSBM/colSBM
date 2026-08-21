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
        message("Using provided init for the model at Q=2")
      }
      start_model <- fit_init
    } else {
      if (verbose) {
        message("Initialising spectral clustering")
      }
      spectral_inits <- lapply(netlist, spectral_clustering, K = 2)
      spectral_taus <- lapply(spectral_inits, .one_hot, Q = 2)
      spectral_taus <- lapply(spectral_taus, function(mat) {
        mat[mat < 1e-6] <- 1e-6
        mat[mat > 1 - 1e-6] <- 1 - 1e-6
        mat / rowSums(mat)
      })
      if (verbose) {
        message("Fitting the model at Q=2")
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
        message("Model fitted !")
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

      if (max(sapply(splitted_model_list, function(mod) mod$BICL)) > best_bicl) {
        if (verbose) {
          say("BICL criterion increased in forward phase, going for backward")
        }
        best_bicl <- max(sapply(splitted_model_list, function(mod) mod$BICL))
      } else {
        if (verbose) {
          say("BICL criterion did not increase in forward phase, stopping optimization.")
        }
        break
      }
      if (length(model_list) == length(splitted_model_list)) {
        print(length(model_list))
        model_list <- compare_model_lists(model_list, splitted_model_list)
        print(length(model_list))
      } else {
        model_list <- splitted_model_list
      }

      model_list <- cpp_backward(model_list = model_list, netlist = netlist, Q_min = Q_min, distribution = distribution, free_mixture = free_mixture, free_density = free_density, fit_opts = fit_opts, verbose = verbose)

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
    model_list <- compare_model_lists(model_list, splitted_model_list)

    return(model_list)
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
      if (!is.null(iteration)) {}
      message("Computing splits for Q_start=", Q_start)
    }

    # Compute splits from the initial model
    next_splits <- purrr::transpose(lapply(seq(M), function(m) split_clust(X = netlist[[m]], Z = start_model[["Z"]][[m]], Q = Q_start)))
    model_list <- list(start_model)
    for (Q in seq(3, Q_max)) {
      if (verbose) {
        message("Fitting models at Q=", Q)
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
      })
      best_model_idx <- which.max(sapply(models, function(model) model[["BICL"]]))
      if (verbose) {
        message("Best splitted model index is ", best_model_idx, " with BICL = ", models[[best_model_idx]][["BICL"]])
      }
      model_list <- append(model_list, models[[best_model_idx]])
      if (verbose) {
        message("Computing next splits")
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
      message("Starting backward phase from Q=", Q_max, " down to Q=", Q_min)
    }

    # Initialize the result list with the models from the forward phase
    backward_model_list <- model_list

    # Process from Q_max down to Q_min + 1
    for (Q in seq(Q_max, Q_min + 1, by = -1)) {
      if (verbose) {
        message("Processing merges for Q=", Q)
      }

      # Compute merges using merge_clust function
      Z_merges <- purrr::transpose(lapply(seq(M), function(m) merge_clust(Z = model_list[[Q - 1]][["Z"]][[m]], Q = Q)))

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
      })

      # Select the best model among the merged ones
      best_model_idx <- which.max(sapply(merged_models, function(model) model[["BICL"]]))

      if (verbose) {
        message("Best merged model index is ", best_model_idx, " with BICL = ", merged_models[[best_model_idx]][["BICL"]])
      }

      # Add the best merged model to the backward list
      backward_model_list[[Q - 1]] <- merged_models[[best_model_idx]]
    }

    return(backward_model_list)
  }

compare_model_lists <- function(model_list1, model_list2) {
  bicl_1 <- sapply(model_list1, function(mod) mod[["BICL"]])
  bicl_2 <- sapply(model_list1, function(mod) mod[["BICL"]])

  idx_max_bicl <- ifelse(bicl_1 >= bicl_2, 1, 2)

  model_list_merged <- lapply(seq_along(idx_max_bicl), function(list_idx) {
    if (idx_max_bicl[list_idx] == 1) {
      return(model_list1[[list_idx]])
    }
    return(model_list2[[list_idx]])
  })
  return(model_list_merged)
}
