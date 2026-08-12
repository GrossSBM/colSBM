library(microbenchmark)
M <- 4
true_pi <- c(0.1, 0.3, 0.6)
nbNodes <- 5000
K <- length(true_pi)
true_alpha <- diag(c(0.9, 0.5, 0.1))
sbm_sampler <- sbm::sampleSimpleSBM(nbNodes = nbNodes, blockProp = true_pi, connectParam = list(mean = true_alpha), directed = TRUE)
set.seed(1234)
sbm_realisations <- lapply(seq(M), sbm_sampler[["rNetwork"]])
A <- lapply(seq(M), function(m) sbm_realisations[[m]][["networkData"]])

# A spectral clustering for the start
spectral_labels <- lapply(A, spectral_clustering, K = 3)

spectral_taus <- lapply(spectral_labels, .one_hot, Q = 3)

spectral_taus <- lapply(spectral_taus, function(tau) {
    tau[tau < 1e-6] <- 1e-6
    tau[tau > 1 - 1e-6] <- 1 - 1e-6
    tau <- tau / rowSums(tau)
})

bench <- microbenchmark(
    "R" = { # colSBM R logic
        R_fit <- fitSimpleSBMPop$new(A = A, Q = 3L, Z = spectral_taus, directed = TRUE, free_mixture = FALSE, free_density = FALSE, init_method = "spectral", fit_opts = default_fit_opts_unipartite())
        R_fit$optimize()
    },
    "C++" = {
        # C++ logic
        ptr <- colsbm_create(A = A, Q = 3L, tau = spectral_taus, directed = TRUE, distribution = "bernoulli", free_mixture = FALSE, free_density = FALSE)
        colsbm_optimize(ptr, max_step = 100L, tol = 1e-3)
    }, times = 10L
)
