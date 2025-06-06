#' Generate a unipartite network
#'
#' @param n the number of nodes
#' @param pi a vector of probability to belong to the clusters
#' @param alpha the matrix of connectivity between two clusters
#' @param distribution the emission distribution, either "bernoulli" or
#' "poisson"
#' @param return_memberships Boolean, should return memberships or not.
#'
#' @return An adjacency matrix

#' @noMd
#' @noRd
generate_unipartite_network <- function(
    n, pi, alpha,
    distribution = "bernoulli",
    return_memberships = FALSE) {
  stopifnot(
    "All alpha coefficients must be positive" = all(alpha >= 0),
    "With bernoulli, the alpha must be between 0 and 1" = (
      distribution == "poisson" |
        (distribution == "bernoulli" & all(alpha >= 0L & alpha <= 1L))),
    "All pi must be between 0 and 1" = all(pi >= 0L & pi <= 1L),
    "Pi must sum to one" = all.equal(sum(pi), 1L)
  )
  cluster_memberships <- rmultinom(n, size = 1, prob = pi)
  node_node_interaction_parameter <- t(cluster_memberships) %*% alpha %*% cluster_memberships

  # Here we switch on the distributions
  adjacency_matrix <- matrix(
    switch(distribution,
      "bernoulli" = {
        rbinom(length(node_node_interaction_parameter),
          size = 1, prob = node_node_interaction_parameter
        )
      },
      "poisson" = {
        rpois(length(node_node_interaction_parameter),
          lambda = node_node_interaction_parameter
        )
      },
      stop("distribution must be one of either 'bernoulli' or 'poisson'")
    ),
    nrow = nrow(node_node_interaction_parameter)
  )

  if (return_memberships) {
    return(list(
      adjacency_matrix = adjacency_matrix,
      block_memberships = cluster_memberships
    ))
  } else {
    return(adjacency_matrix)
  }
}

#' Generate collection of unipartite
#'
#' @param n the number of nodes or a vector of the nodes per network
#' @param pi a vector of probability to belong to the clusters
#' @param alpha the matrix of connectivity between two clusters
#' @param M the number of networks to generate
#' @param distribution the emission distribution, either "bernoulli" or
#' "poisson"
#' @param return_memberships Boolean, should return memberships or not.
#'
#' @details If n is a single value, this value will be replicated for each of
#' the M networks. If it is a vector, it must be of size M, specifying the
#' number of nodes for each network.
#'
#' @return A list of M lists, which contains : $adjacency_matrix, $clustering
#'
#' @export
generate_unipartite_collection <- function(
    n, pi, alpha, M,
    distribution = "bernoulli",
    return_memberships = FALSE) {
  if (length(n) == 1) {
    n <- rep(n, M)
  }

  # Check if n is the correct length
  if (length(n) != M) {
    stop(
      "The length of n is not correct ! It should be : ",
      M, " values and it is ", length(n)
    )
  }
  # Generate the networks
  out <- lapply(seq.int(M), function(m) {
    generate_unipartite_network(
      n = n[[m]],
      pi = pi,
      alpha = alpha,
      distribution = distribution,
      return_memberships = return_memberships
    )
  })
  return(out)
}

#' Generate a bipartite network
#'
#' @param nr the number of row nodes
#' @param nc the number of col nodes
#' @param pi a vector of probability to belong to the row clusters
#' @param rho a vector of probability to belong to the columns clusters
#' @param alpha the matrix of connectivity between two clusters
#'
#' @return An incidence matrix
#'
#' @noMd
#' @noRd
generate_bipartite_network <- function(
    nr, nc, pi, rho, alpha, distribution = "bernoulli",
    return_memberships = FALSE) {
  stopifnot(
    "All alpha coefficients must be positive" = all(alpha >= 0L),
    "With bernoulli, the alpha must be between 0 and 1" = (
      distribution == "poisson" |
        (distribution == "bernoulli" & all(alpha >= 0L & alpha <= 1L))),
    "All pi must be between 0 and 1" = all(pi >= 0L & pi <= 1L),
    "Pi must sum to one" = all.equal(sum(pi), 1L),
    "All rho must be between 0 and 1" = all(rho >= 0L & rho <= 1L),
    "Rho must sum to one" = all.equal(sum(rho), 1L)
  )

  rowblocks_memberships <- rmultinom(nr, size = 1, prob = pi)
  colblocks_memberships <- rmultinom(nc, size = 1, prob = rho)
  node_node_interaction_parameter <- t(rowblocks_memberships) %*%
    alpha %*%
    colblocks_memberships
  incidence_matrix <- matrix(
    switch(distribution,
      "poisson" = {
        rpois(length(node_node_interaction_parameter),
          lambda = node_node_interaction_parameter
        )
      },
      "bernoulli" = {
        rbinom(length(node_node_interaction_parameter),
          size = 1, prob = node_node_interaction_parameter
        )
      },
      stop("distribution must be one of either 'bernoulli' or 'poisson'")
    ),
    nrow = nrow(node_node_interaction_parameter)
  )

  if (return_memberships) {
    return(list(
      incidence_matrix = incidence_matrix,
      row_blockmemberships = as.vector(
        c(seq.int(length(pi))) %*% rowblocks_memberships
      ), # We reverse the one hot encoding
      col_blockmemberships = as.vector(
        c(seq.int(length(rho))) %*% colblocks_memberships
      ) # We reverse the one hot encoding
    ))
  } else {
    return(incidence_matrix)
  }
}

#' Generate collection of bipartite networks
#'
#' @param nr the number of row nodes  or a vector specifying the
#' number of row nodes for each of the M networks
#' @param nc the number of column nodes  or a vector specifying the
#' number of column nodes for each of the M networks
#' @param pi a vector of probability to belong to the row clusters
#' @param rho a vector of probability to belong to the columns clusters
#' @param alpha the matrix of connectivity between two clusters
#' @param M the number of networks to generate
#' @param model the colBiSBM model to use. Available: "iid", "pi", "rho",
#' "pirho"
#' @param distribution the emission distribution : "bernoulli" or "poisson"
#' @param return_memberships a boolean which choose whether the function returns
#' a list containing the memberships and the incidence matrices or just the
#' incidence matrices. Defaults to FALSE, only the matrices are returned.
#'
#' @details the model parameters if set to any other than iid will shuffle the
#' provided pi and rho
#'
#' @return A list of M lists, which contains : $incidence_matrix, $row_blockmemberships, $col_blockmemberships
#'
#' @export
generate_bipartite_collection <- function(
    nr, nc, pi, rho, alpha, M,
    model = "iid",
    distribution = "bernoulli",
    return_memberships = FALSE) {
  out <- list()

  # Check if nr and nc are vectors
  if (length(nr) == 1) {
    nr <- rep(nr, M)
  }
  if (length(nc) == 1) {
    nc <- rep(nc, M)
  }

  # Check if nr and nc are the correct length
  if (length(nr) != M) {
    stop(
      "The length of nr is not correct ! It should be : ",
      M, " values and it is ", length(nr)
    )
  }
  if (length(nc) != M) {
    stop(
      "The length of nc is not correct ! It should be : ",
      M, " values and it is ", length(nc)
    )
  }

  switch(model,
    "iid" = {
      out <- lapply(seq.int(M), function(m) {
        generate_bipartite_network(
          nr = nr[[m]],
          nc = nc[[m]],
          pi = pi,
          rho = rho,
          alpha = alpha,
          distribution = distribution,
          return_memberships = return_memberships
        )
      })
    },
    "pi" = {
      out <- lapply(seq.int(M), function(m) {
        generate_bipartite_network(
          nr = nr[[m]],
          nc = nc[[m]],
          pi = sample(pi),
          rho = rho,
          alpha = alpha,
          distribution = distribution,
          return_memberships = return_memberships
        )
      })
    },
    "rho" = {
      out <- lapply(seq.int(M), function(m) {
        generate_bipartite_network(
          nr = nr[[m]],
          nc = nc[[m]],
          pi = pi,
          rho = sample(rho),
          alpha = alpha,
          distribution = distribution,
          return_memberships = return_memberships
        )
      })
    },
    "pirho" = {
      out <- lapply(seq.int(M), function(m) {
        generate_bipartite_network(
          nr = nr[[m]],
          nc = nc[[m]],
          pi = sample(pi),
          rho = sample(rho),
          alpha = alpha,
          distribution = distribution,
          return_memberships = return_memberships
        )
      })
    },
    stop("Error unknown model. Must be one of : iid, pi, rho, pirho.")
  )

  return(out)
}

#' Perform a spectral clustering
#'
#' @importFrom stats kmeans
#'
#' @param X an adjacency matrix
#' @param K the number of clusters
#' @param kmeans.nstart the number of random starts for the kmeans algorithm.
#' Defaults to 400. Ensures consistency of the results.
#' @param kmeans.iter.max the maximum number of iterations for the kmeans
#' algorithm. Defaults to 50.
#'
#' @return A vector : The clusters labels
#'
#' @keywords internal
#'
#' @return A vector : The clusters labels
spectral_clustering <- function(X, K, kmeans.nstart = 400L, kmeans.iter.max = 50L) {
  X <- as.matrix(X)
  n <- nrow(X)
  if (K == 1) {
    return(rep(1L, nrow(X)))
  }
  X[X == -1] <- NA
  isolated <- which(rowSums(X, na.rm = TRUE) == 0)
  connected <- setdiff(seq(n), isolated)
  X <- X[connected, connected]
  if (!isSymmetric.matrix(X)) {
    X <- 1 * ((X + t(X)) > 0) # .5 * (X + t(X))
  }
  if (nrow(X) < 3) {
    return(rep(1, n))
  }
  D_moins1_2 <- diag(1 / sqrt(rowSums(X, na.rm = TRUE) + 1))
  X[is.na(X)] <- mean(X, na.rm = TRUE)
  Labs <- D_moins1_2 %*% X %*% D_moins1_2
  specabs <- eigen(Labs, symmetric = TRUE)
  if (K >= nrow(X)) {
    message("Too many clusters for Spectral Clustering")
    K_old <- K
    K <- nrow(X) - 1
  }
  if (K >= 2) {
    index <- rev(order(abs(specabs$values)))[1:K]
    U <- specabs$vectors[, index]
    U <- U / rowSums(U**2)**(1 / 2)
    U[is.na(U)] <- 0
    U[is.nan(U)] <- 0
    U[is.infinite(U)] <- 0
    cl <- try(
      expr = {
        stats::kmeans(U, K,
          iter.max = kmeans.iter.max,
          nstart = kmeans.nstart
        )$cluster
      },
      silent = TRUE
    )
    while (inherits(cl, "try-error")) {
      cl <- try(
        expr = {
          stats::kmeans(U, K,
            iter.max = kmeans.iter.max,
            nstart = kmeans.nstart
          )$cluster
        },
        silent = TRUE
      )
    }
  } else {
    cl <- rep(1, nrow(X))
  }
  clustering <- rep(1, n)
  clustering[connected] <- cl
  clustering[isolated] <- which.min(rowsum(rowSums(X, na.rm = TRUE), cl))
  clustering[isolated] <- which.min(rowsum(rowSums(X, na.rm = TRUE), cl))
  return(clustering)
}



#' Perform a spectral bi-clustering, clusters by row
#' and by columns independently
#'
#' Relies on the spectral_clustering function defined above
#'
#' @param A a bipartite adjacency matrix
#' @param Q the two numbers of clusters
#' @inheritParams spectral_clustering
#'
#' @return A list of two vectors : The clusters labels. They are accessed using
#' $row_clustering and $col_clustering
#'
#' @return A list of two vectors : The clusters labels.
#' They are accessed using $row_clustering and $col_clustering
#'
#' @keywords internal
#'
spectral_biclustering <- function(A, Q, kmeans.nstart = 400L, kmeans.iter.max = 50L) {
  # Trivial clustering : everyone is part of the cluster
  if (all(Q == c(1, 1))) {
    return(list(
      row_clustering = rep(1, nrow(A)),
      col_clustering = rep(1, ncol(A))
    ))
  }

  # Extracts the number of clusters
  K1 <- Q[1] # Row clusters
  K2 <- Q[2] # Column clusters

  row_adjacency_matrix <- tcrossprod(A)
  row_clustering <- spectral_clustering(row_adjacency_matrix,
    K1,
    kmeans.nstart = kmeans.nstart,
    kmeans.iter.max = kmeans.iter.max
  )

  col_adjacency_matrix <- crossprod(A)
  col_clustering <- spectral_clustering(col_adjacency_matrix,
    K2,
    kmeans.nstart = kmeans.nstart,
    kmeans.iter.max = kmeans.iter.max
  )

  return(list(row_clustering = row_clustering, col_clustering = col_clustering))
}

# TODO : Modify the algorithm to use the rectangular matrix and its transpose
bipartite_hierarchic_clustering <- function(X, K) {
  # Trivial clustering : everyone is part of the cluster
  if (all(K == c(1, 1))) {
    return(list(
      row_clustering = rep(1, nrow(X)),
      col_clustering = rep(1, ncol(X))
    ))
  }

  # Extracts the number of clusters
  K1 <- K[1] # Row clusters
  K2 <- K[2] # Column clusters

  row_adjacency_matrix <- tcrossprod(X)
  row_clustering <- hierarClust(X, K1)

  col_adjacency_matrix <- crossprod(X)
  col_clustering <- hierarClust(t(X), K2)

  return(list(row_clustering = row_clustering, col_clustering = col_clustering))
}

#' Perform a Hierarchical Clustering
#' @importFrom stats cutree dist hclust
#' @param X An Adjacency Matrix
#' @param K the number of wanted clusters
#'
#' @noMd
#' @noRd
#'
#' @return A vector : The clusters labels
hierarClust <- function(X, K) {
  if (K == 1) {
    return(rep(1L, nrow(X)))
  }
  diss <- cluster::daisy(x = X, metric = "manhattan", warnBin = FALSE)
  if (anyNA(diss)) {
    return(rep(1L, nrow(X)))
  } else {
    clust <- cluster::agnes(x = X, metric = "manhattan", method = "ward")
  }
  # clust    <- stats::hclust(d = distance , method = "ward.D2")
  return(stats::cutree(tree = clust, k = K))
}

#' Split a list of clusters
#'
#' @param X an adjacency matrix
#' @param Z a vector of cluster memberships
#' @param Q The number of maximal clusters
#'
#' @noMd
#' @noRd
#'
#' @return A list of Q clustering of Q+1 clusters
split_clust <- function(X, Z, Q, is_bipartite = FALSE) {
  X <- matrix(X,
    nrow = nrow(X),
    ncol = ncol(X)
  )
  Z_split <- lapply(
    seq(Q),
    FUN = function(q) {
      if (sum(Z == q) <= 3) {
        return(Z)
      }
      Z_new <- Z
      mx <- ifelse(is_bipartite, mean(X[Z == q, ], na.rm = TRUE), mean(X[Z == q, Z == q], na.rm = TRUE))
      X[is.na(X)] <- mx
      if (isSymmetric.matrix(X) || is_bipartite) {
        Xsub <- X[Z == q, , drop = FALSE]
      } else {
        Xsub <- cbind(X[Z == q, , drop = FALSE], t(X[, Z == q, drop = FALSE]))
      }
      Xsub <- as.matrix(Xsub)
      if (nrow(unique(Xsub, MARGIN = 1)) <= 3) {
        return(Z)
      }
      C <- stats::kmeans(x = Xsub, centers = 2)$cluster
      # C  <-  hierarClust(Xsub, 2)# + Q#stats::kmeans(x = .5 * (X[Z==q,] + t(X[,Z==q])), centers = 2)$cluster + Q
      if (length(unique(C)) == 2) {
        if (!is_bipartite) {
          p1 <- c(
            mean(X[Z == q, , drop = FALSE][C == 1, ]),
            mean(X[, Z == q, drop = FALSE][, C == 1])
          )
          p2 <- c(
            mean(X[Z == q, , drop = FALSE][C == 2, ]),
            mean(X[, Z == q, drop = FALSE][, C == 2])
          )
        } else {
          p1 <- c(
            mean(X[Z == q, , drop = FALSE][C == 1, ])
          )
          p2 <- c(
            mean(X[Z == q, , drop = FALSE][C == 2, ])
          )
        }
        c <- which.max(abs(p1 - p2))
        md <- sample(
          x = 2, size = 2, replace = FALSE,
          prob = c(
            max(1 / ncol(Xsub), p1[c]),
            max(1 / ncol(Xsub), p2[c])
          )
        )
        Z_new[Z == q][C == which.min(md)] <- Q + 1
      }
      # Z_new[Z==q]  <-  stats::kmeans(x = Xsub, centers = 2)$cluster + Q
      # Z_new[Z_new == Q + 2]  <-  q
      return(Z_new)
    }
  )
  Z_split <- Z_split[which(sapply(X = Z_split, FUN = function(x) !is.null(x)))]
  return(Z_split)
}

#' Merge a list of clusters
#'
#' @importFrom utils combn
#'
#' @param Z a vector of cluster memberships
#' @param Q the number of original clusters
#'
#' @return A list of Q(Q-1)/2 clustering of Q-1 clusters
#'
#'  @noMd
#' @noRd
merge_clust <- function(Z, Q) {
  Z_merge <- lapply(
    X = 1:choose(Q, 2),
    FUN = function(q) {
      Z[Z == utils::combn(Q, 2)[2, q]] <- utils::combn(Q, 2)[1, q]
      if (utils::combn(Q, 2)[2, q] < Q) {
        Z[Z > utils::combn(Q, 2)[2, q]] <- Z[Z > utils::combn(Q, 2)[2, q]] - 1
      }
      return(Z)
    }
  )
  return(Z_merge)
}



#
# Fonction interne a VEM
#

F.bern <- function(X, alpha, tau) {
  return(X %*% tau %*% log(alpha) +
    (1 - X - diag(1, nrow(X))) %*% tau %*% log(1 - alpha))
}

rotate <- function(x) t(apply(x, 2, rev))


dist_param <- function(param, param_old) {
  sqrt(sum((param - param_old)**2))
}

#' A simple function that returns the default global options for the unipartite
#' @noRd
default_global_opts_unipartite <- function(netlist) {
  n <- sapply(netlist, nrow)
  list(
    Q_min = 1L,
    Q_max = floor(log(sum(n))) + 2,
    sbm_init = TRUE,
    spectral_init = TRUE,
    nb_init = 10L,
    nb_models = 5L,
    depth = 3L,
    plot_details = 1L,
    max_pass = 10L,
    verbosity = 0L,
    nb_cores = 1L
  )
}

#' A simple function that returns default fit options for unipartite
#' @noRd
default_fit_opts_unipartite <- function() {
  list(
    algo_ve = "fp",
    approx_pois = FALSE,
    minibatch = TRUE,
    verbosity = 0L
  )
}

#' A simple function that returns the default global options for the bipartite
#' @noRd
default_global_opts_bipartite <- function(netlist) {
  list(
    Q1_min = 1L,
    Q2_min = 1L,
    Q1_max = floor(log(sum(sapply(netlist, function(A) nrow(A)))) + 2),
    Q2_max = floor(log(sum(sapply(netlist, function(A) ncol(A)))) + 2),
    nb_init = 10L,
    nb_models = 5L,
    backend = "future",
    depth = 1L,
    plot_details = 1L,
    max_pass = 10L,
    verbosity = 1L,
    nb_cores = 1L,
    compare_stored = TRUE
  )
}

#' A simple function that returns default fit options for bipartite
#' @noRd
default_fit_opts_bipartite <- function() {
  list(
    algo_ve = "fp",
    minibatch = TRUE,
    verbosity = 0,
    tolerance = 1e-6,
    max_vem_steps = 3000L,
    greedy_exploration_max_steps = 50,
    greedy_exploration_max_steps_without_improvement = 5,
    kmeans_nstart = 400L,
    kmeans_iter_max = 100L,
    penalty_factor = 0.5
  )
}

#' Title
#'
#' @param X An adjacency matrix
#' @param K An integer, the number of folds
#'
#' @return A matrix of the same size than X with class integer as coefficient
#'
#' @noMd
#' @noRd
build_fold_matrix <- function(X, K) {
  n <- ncol(X)
  arrange <- sample.int(n)
  labels <- cut(seq(n), breaks = K, labels = FALSE)
  fold_matrix <- diag(n)
  for (i in seq(n)) {
    fold_matrix[i, ] <- (labels + labels[i]) %% K
  }
  fold_matrix <- fold_matrix[arrange, arrange] + 1
  diag(fold_matrix) <- 0
  return(fold_matrix)
}

#' Compute BIC-L of a provided partition (list of bmpop or bisbmpop)
#'
#' @param partition A list of (or a single) bmpop, bisbmpop, fitBipartiteSBMPop,
#' fitSimpleSBMPop objects
#' @param penalty_factor A numeric value, the penalty factor to use for the
#' computation of the BIC-L. Defaults to 0.5.
#' @param verbose A boolean, should the function print additional information.
#' Defaults to TRUE.
#'
#' @return A numeric value, the BIC-L of the partition
#' @noRd
compute_bicl_partition <- function(partition, penalty_factor = 0.5, verbose = TRUE) {
  if (inherits(partition, "bmpop") | inherits(partition, "bisbmpop")) {
    if (verbose) {
      cli::cli_alert_info(
        "A {.type {partition}} object was provided. The BIC-L is computed from the best fit."
      )
    }
    return(partition$best_fit$compute_BICL(
      penalty_factor = penalty_factor,
      store = FALSE
    ))
  }
  if (inherits(partition, "fitBipartiteSBMPop") | inherits(partition, "fitSimpleSBMPop")) {
    return(partition$compute_BICL(
      penalty_factor = penalty_factor,
      store = FALSE
    ))
  }
  if (inherits(partition, "list")) {
    if (all(sapply(partition, inherits, "bmpop") | sapply(partition, inherits, "bisbmpop"))) {
      if (verbose) {
        cli::cli_alert_info(
          "A list of {.type {partition[[1]]}} objects was provided. The BIC-L is computed from the best fit."
        )
      }
      return(sum(sapply(partition, function(col) {
        col$best_fit$compute_BICL(
          penalty_factor = penalty_factor,
          store = FALSE
        )
      })))
    }
    if (all(sapply(partition, inherits, "fitBipartiteSBMPop") | sapply(partition, inherits, "fitSimpleSBMPop"))) {
      return(sum(sapply(partition, function(col) {
        col$compute_BICL(
          penalty_factor = penalty_factor,
          store = FALSE
        )
      })))
    }
  }
  stop("The provided partition is not a valid object for BIC-L computation.")
}

.xlogx <- function(x) {
  ifelse(x < 2 * .Machine$double.eps, 0, x * log(x))
}
.xlogy <- function(x, y, eps = NULL) {
  ifelse(x < 2 * .Machine$double.eps, 0, x * .log(y, eps = eps))
}
.quadform <- function(x, y) tcrossprod(x %*% y, x)
.tquadform <- function(x, y) crossprod(x, y %*% x)
logistic <- function(x) 1 / (1 + exp(-x))
logit <- function(x) log(x / (1 - x))
.logit <- function(x, eps = NULL) {
  if (is.null(eps)) {
    res <- log(x / (1 - x))
  } else {
    res <- log(pmax(pmin(x, 1 - eps), eps) / pmax(pmin(1 - x, 1 - eps), eps))
  }
  return(res)
}


.threshold <- function(x, eps = 1e-9) {
  #  x <- .softmax(x)
  x[x < eps] <- eps
  x[x > 1 - eps] <- 1 - eps
  x <- x / .rowSums(x, nrow(x), ncol(x))
  x
}

.softmax <- function(x) {
  # x_max <- apply(x, 1, max)
  # OPTIM
  x_max <- matrixStats::rowMaxs(x)
  x <- exp(x - x_max)
  x <- x / .rowSums(x, nrow(x), ncol(x))
  x
}
.log <- function(x, eps = NULL) {
  if (is.null(eps)) {
    res <- log(x)
  } else {
    # OPTIM
    # x[x >= 1 - eps] <- 1 - eps
    x[x <= eps] <- eps
    res <- log(x)
  }
  return(res)
}
.one_hot <- function(x, Q) {
  O <- matrix(0, length(x), Q)
  O[cbind(seq.int(length(x)), x)] <- 1
  return(O)
}

.rev_one_hot <- function(X) {
  return(as.vector(max.col(X)))
}
